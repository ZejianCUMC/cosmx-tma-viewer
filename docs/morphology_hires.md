# CGC viewer — increase morphology (rawImage) resolution for a sample

Reusable recipe to sharpen the "RAW morphology — CellComposite image" layer of any
per-patient viewer. Applied successfully to two wide wide cores **@ 0.25 µm/px**
(2026-09-14). Use on-demand for specific samples; NOT applied cohort-wide (would add ~146 MB).

## Root cause of low-res morphology
`stitch_one.py` downscales each mosaic to a fixed **`target_long = 2500 px`** over the FULL
FOV bounding box. For cores whose FOVs are spread across the slide (wide µm extent) the 2500 px
is stretched thin → ~5 µm/px (blurry). Compact cores (FOVs clustered) are already ~0.6 µm/px.
`degap`/`regrid` run AFTER the downscale and only remove empty space (they preserve px/µm), so
they cannot recover resolution — and every viewer rebuild reuses the cached low-res stitch.

## The fix (extent-preserving → ONLY morphology changes)
Set `target_long` **adaptively** from the micron span so density hits a target µm/px:
`TARGET = min(CAP, max(2500, round(max(span_x,span_y) * umpx / desired_umpx)))`.
Because `rawExtent` is the micron bbox (independent of target_long), the DOWNSTREAM
build+degap+regrid geometry is bit-identical — only `rawImage` pixel density increases.
Provided scripts: `stitch_one_hires_1umpx.py` (~1 µm/px, CAP 14000) and
`stitch_one_hires_025umpx.py` (~0.25 µm/px, CAP 55000).

## One-time pipeline fix (already applied durably)
`degap_viewer.py` and `regrid_viewer.py` now set `Image.MAX_IMAGE_PIXELS=None` after the PIL
import. Required: a 0.25 µm/px wide-core mosaic (e.g. 50856x15406 = 783 M px) trips
PIL's decompression-bomb limit (~179 M) and degap/regrid silently fail (gap not removed → bloated
viewer). With the limit raised, degap removes the empty gap and compacts to tissue-only.

## Step-by-step (example: unit UID at desired µm/px)
Vars: DEST=<hpc_project_dir>,
VP=result/viewer_pipeline (local), V6=result/majortype_clean/subtype/cgc_hierarchical_labels_v6.tsv.gz
(the LOCKED v6 labels — always use this; verify md5=2a47ac98… vs FINAL_ANNOTATION_LOCK_v6.json).

1. Pick script: 0.25 µm/px = stitch_one_hires_025umpx.py (or 1 µm/px). Ship to a dedicated HPC
   build dir with viewer_units.tsv + a units.txt (target UIDs, one per line) + the sbatch:
     ssh <hpc> "mkdir -p $DEST/rawimg_build_hiresX"
     scp stitch_one_hires_025umpx.py <hpc>:$DEST/rawimg_build_hiresX/stitch_one.py
     scp viewer_units.tsv stitch_morphology_hires.sbatch <hpc>:$DEST/rawimg_build_hiresX/
     printf 'UID\n' | ssh <hpc> "cat > $DEST/rawimg_build_hiresX/units.txt"
2. Submit (resources in stitch_morphology_hires.sbatch = 2 CPU / 32 GB / 40 min; 32 GB needed for
   wide cores at 0.25 µm/px, 8 GB is enough at 1 µm/px or for compact cores):
     ssh <hpc> "cd $DEST/rawimg_build_hiresX && BUILD_DIR=\$PWD sbatch stitch_morphology.sbatch"
   HPC jobs pre-approved (small/<10h). Check: sacct -j <JID> ; expect rawimg_<UID>.json.
3. Pull the hi-res json to a local rawimg dir:
     scp <hpc>:$DEST/rawimg_build_hiresX/rawimg_UID.json  $VP/rawimg_hiresX/
4. Rebuild the viewer with LOCKED v6 labels + the hi-res morphology, then postprocess:
     python build_sample_viewer.py UID --labels $V6 --rawimg-dir $VP/rawimg_hiresX --outdir $VP/viewer_testX
     python degap_viewer.py  $VP/viewer_testX/CGC_UID_sample_viewer.html
     python regrid_viewer.py $VP/viewer_testX/CGC_UID_sample_viewer.html --ncol 3
     python finalize_region_polygon.py $VP/viewer_testX/CGC_UID_sample_viewer.html
5. VERIFY only morphology changed (must print ALL identical = True, rawExtent equal = True):
     python morphology_hires_method/verify_only_morphology_changed.py \
       $VP/viewer_testX/CGC_UID_sample_viewer.html  <DBX>/CGC_UID_sample_viewer.html
6. Deploy (the deploy directory = share/CGC/viewer/): local-backup the current the deploy directory file FIRST (NO _bak in
   the deploy directory), cp test → the deploy directory + local viewer/, and update the cache so future rebuilds stay hi-res:
     cp -p rawimg_hiresX/rawimg_UID.json rawimg_all/rawimg_UID.json   # backup old locally first

## Density / size guidance (calibrated)
- Morphology JPEG ≈ **0.55 MB per FOV** at 0.25 µm/px (q82); ~0.15 MB/FOV at 1 µm/px.
- Final viewer size = existing all-genes size (cell-count-driven) + morphology. Cohort-wide 0.25 µm/px
  ≈ 667 → 814 MB (+146 MB); a few big cores (a 36-FOV core → ~97 MB).
- desired_umpx choices: 1.0 (safe/small, 8 GB), 0.5 (2× sharper, 16 GB), 0.25 (4× sharper, 32 GB),
  native 0.12 µm/px (max; large files + needs gap-aware stitch for wide cores).
- WIDE cores waste the pixel budget on the empty inter-FOV gap (huge intermediate mosaic). For
  near-native res on those, the cleaner long-term fix is a gap-aware compact stitch (stitch only
  tissue FOVs); not implemented here — the extent-preserving approach above is what is validated.

## Guarantee
Only `DATA.rawImage` changes. `rawExtent`, coordinates, cell-type labels (v6/IFN-response), all
6,175 genes, IF channels, important-region overlays and UI stay byte-identical (verified by MD5, step 5).
