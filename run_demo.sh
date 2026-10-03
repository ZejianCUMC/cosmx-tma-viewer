#!/usr/bin/env bash
# Synthetic end-to-end demo: data -> viewers -> landing page -> browser checks. No real data involved.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python}
$PY input/make_demo_data.py demo/data
$PY input/build_viewer.py --cells demo/data/cells.tsv.gz --cores demo/data/cores.tsv --units demo/data/units.tsv \
    --fov-positions demo/data/fov_positions.tsv --config input/demo_config.json --out demo/out \
    --counts-dir demo/data/counts --fov-images demo/data/fov_images --fov-images-hires demo/data/fov_images_hires \
    --he-images demo/data/he --regions demo/data/regions.tsv
$PY input/build_index.py demo/out --flags demo/data/flags.tsv --config input/demo_config.json
if [ -n "${CHROME:-}" ] || command -v google-chrome >/dev/null || command -v chromium >/dev/null || command -v chromium-browser >/dev/null; then
  $PY input/verify_viewer.py demo/out/viewer_DEMO_*.html ${CHROME:+--chrome "$CHROME"} --report demo/verify.json
else
  echo "Chrome/Chromium not found: skipped browser checks (set CHROME=/path/to/chrome)"
fi
echo "open demo/out/index.html"
