#!/usr/bin/env bash
# =============================================================================
# Reproduce ALL per-unit spatial viewers  (ZERO LLM tokens per patient)
# =============================================================================
# Hybrid pipeline: cheap parts LOCAL, heavy morphology stitch on HPC sbatch.
# The LLM built the generator ONCE; every viewer below is produced by scripts.
#
#   stage manifest : LOCAL  build viewer_units.tsv + units.txt from the object
#   stage ship     : LOCAL  scp scripts + manifest + units.txt -> HPC build dir
#   stage stitch   : HPC    sbatch ARRAY (1 task/unit) -> rawimg_<unit>.json
#   stage pull     : LOCAL  scp rawimg_*.json back down
#   stage viewers  : LOCAL  build_sample_viewer.py --all -> <prefix><unit>_*.html
#   stage verify   : LOCAL  verify_overlay.py on every viewer (alignment QC)
#
# Usage:
#   ./run_all.sh manifest ship stitch pull viewers verify     # step by step
#   ./run_all.sh all                                          # everything
# One (patient,slide) unit == one array task == one viewer. Multi-slide patients
# A unit with cores on both slides yields one viewer per slide.
# =============================================================================
set -euo pipefail

PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYLOCAL=<env>/bin/python  # e.g. a conda env named scvi
HPC=<hpc_login>
DEST=<hpc_project_dir>
BUILD=$DEST/rawimg_build
VIEWERS="$PIPE/viewers"

stage_manifest() {
  echo "== [manifest] =="
  "$PYLOCAL" "$PIPE/build_sample_viewer.py" --dump-manifest
  tail -n +2 "$PIPE/viewer_units.tsv" | cut -f1 > "$PIPE/units.txt"
  echo "units.txt: $(wc -l < "$PIPE/units.txt") units"
}

stage_ship() {
  echo "== [ship] scp -> $BUILD =="
  ssh "$HPC" "mkdir -p $BUILD"
  scp "$PIPE/stitch_one.py" "$PIPE/stitch_morphology.sbatch" \
      "$PIPE/viewer_units.tsv" "$PIPE/units.txt" "$HPC:$BUILD/"
}

stage_stitch() {
  local N JID
  N=$(ssh "$HPC" "wc -l < $BUILD/units.txt")
  echo "== [stitch] sbatch --array=1-$N  (2 CPU / 8 GB / 20 min each) =="
  echo "   total est. compute: $N tasks x ~15-60 s = a few CPU-minutes"
  JID=$(ssh "$HPC" "cd $BUILD && BUILD_DIR=$BUILD sbatch --parsable --array=1-$N stitch_morphology.sbatch")
  echo "   submitted array job $JID ; waiting for completion..."
  # quoted heredoc -> nothing expands locally; JID/BUILD passed as remote env vars
  ssh "$HPC" "JID=$JID BUILD=$BUILD bash -s" <<'REMOTE'
    for i in $(seq 1 120); do
      R=$(squeue -j "$JID" -h -o '%T' 2>/dev/null | wc -l)
      [ "$R" -eq 0 ] && { echo 'array drained'; break; }
      echo "  still running: $R task(s)"; sleep 15
    done
    echo "rawimg files: $(ls $BUILD/rawimg_*.json 2>/dev/null | wc -l)"
REMOTE
}

stage_pull() {
  echo "== [pull] scp rawimg_*.json -> $PIPE =="
  scp "$HPC:$BUILD/rawimg_*.json" "$PIPE/"
  echo "local rawimg files: $(ls "$PIPE"/rawimg_*.json 2>/dev/null | wc -l)"
}

stage_viewers() {
  echo "== [viewers] build ALL (loads matrix once) =="
  "$PYLOCAL" "$PIPE/build_sample_viewer.py" --all --rawimg-dir "$PIPE/rawimg_all" --outdir "$VIEWERS"
  echo "viewers: $(ls "$VIEWERS"/<prefix>*_sample_viewer.html 2>/dev/null | wc -l)"
}

stage_verify() {
  echo "== [verify] alignment QC on every viewer =="
  for v in "$VIEWERS"/<prefix>*_sample_viewer.html; do
    echo "--- $(basename "$v") ---"
    "$PYLOCAL" "$PIPE/verify_overlay.py" "$v" "${v%.html}_align.png" | sed -n '1,4p' || echo "  (no morphology / skip)"
  done
}

stage_postprocess() {
  echo "== [postprocess] degap → regrid → finalize region =="
  for v in "$VIEWERS"/<prefix>*_sample_viewer.html; do
    "$PYLOCAL" "$PIPE/degap_viewer.py"          "$v" 2>&1 | tail -1
    "$PYLOCAL" "$PIPE/regrid_viewer.py"          "$v" --ncol 3 2>&1 | tail -1
    "$PYLOCAL" "$PIPE/finalize_region_polygon.py"   "$v" 2>&1 | tail -1
  done
}

stage_deploy() {
  local DBX=<deploy_dir>
  local TS=$(date +%Y%m%d_%H%M%S)
  echo "== [deploy] backup dropbox viewers to _bak_${TS}_pre_deploy + copy new =="
  mkdir -p "$DBX/_bak_${TS}_pre_deploy"
  cp -p "$DBX"/<prefix>*_sample_viewer.html "$DBX/_bak_${TS}_pre_deploy/" 2>/dev/null || true
  cp -p "$VIEWERS"/<prefix>*_sample_viewer.html "$DBX/"
  echo "  deployed: $(ls "$DBX"/<prefix>*_sample_viewer.html | wc -l) files → $DBX"
}

[[ $# -eq 0 ]] && { grep -E '^#( |=)' "$0" | sed 's/^# \{0,1\}//'; exit 0; }
for s in "$@"; do
  case "$s" in
    all) stage_manifest; stage_ship; stage_stitch; stage_pull; stage_viewers; stage_postprocess; stage_verify; stage_deploy;;
    manifest|ship|stitch|pull|viewers|postprocess|verify|deploy) "stage_$s";;
    *) echo "unknown stage: $s"; exit 1;;
  esac
done
