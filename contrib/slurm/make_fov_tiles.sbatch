#!/bin/bash
# Full-resolution FOV tiles for the viewer zoom-in layer (input/make_fov_tiles.py), run next to the raw images.
# Fill the <...> placeholders. ~20 s for 50 FOVs at 2048 px; 2 CPU / 4 GB is ample.
#SBATCH --job-name=viewer_fov_tiles
#SBATCH --partition=<PARTITION>
#SBATCH --account=<ACCOUNT>
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=viewer_fov_tiles_%j.log
set -euo pipefail
PY=<PATH_TO_PYTHON_WITH_PILLOW>
REPO=<PATH_TO_cosmx-tma-viewer>
# Either a list (slide, fov, path) ...
$PY $REPO/input/make_fov_tiles.py tiles --list fovs_with_paths.tsv --tile-px 2048
# ... or one slide folder of CellComposite images:
# $PY $REPO/input/make_fov_tiles.py tiles --source <CellComposite_dir> --slide <SLIDE> --fovs 1,2,3 --tile-px 2048
