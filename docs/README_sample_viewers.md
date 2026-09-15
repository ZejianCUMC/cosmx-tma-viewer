# CosMx per-unit spatial viewers — reproducible generator

Rebuilds **all ~23 patient viewers** (each with the *RAW morphology — CellComposite image*
backdrop) by running scripts. **Zero LLM tokens per patient** — the LLM wrote the generator
once; every viewer is produced deterministically by `build_sample_viewer.py` + the HPC stitch.

## One-command reproduction
```bash
cd <this dir>
./run_all.sh all          # manifest → ship → HPC stitch array → pull → build → verify
```
Or step by step (recommended first time, so you can eyeball the QC):
```bash
./run_all.sh manifest ship stitch pull viewers verify
```

## What each deliverable does
| File | Where | Role |
|------|-------|------|
| `build_sample_viewer.py` | LOCAL (`scvi` env) | Reads the cell object `.h5ad`, builds per-patient DATA, injects the golden template, embeds `rawimg_<unit>.json` if present → `viewers/<prefix><unit>_sample_viewer.html`. |
| `stitch_one.py` | HPC compute node (`GenomicTools` env) | Stitches one unit's FOV CellComposites into a global-µm mosaic (per-FOV vertical flip + **whole-mosaic flipud fix** + downscale 2500px + JPEG q82) → `rawimg_<unit>.json` `{rawImage, rawExtent}`. |
| `stitch_morphology.sbatch` | HPC | SLURM **array**, one task per unit (2 CPU / 8 GB / 20 min). |
| `run_all.sh` | LOCAL | Orchestrates the full reproduction (stages above). |
| `verify_overlay.py` | LOCAL (`scvi` env) | Alignment QC: overlays cells on the embedded mosaic using the **viewer-exact Y-flip** and prints a blob-proof PanCK-vs-tissue check. |
| `sample_viewer_template.html` | — | Golden viewer shell (`const DATA=__DATA__`), carries all evolved features (rawimg gate, virtual-H&E fallback, L2/L3 subtypes, IF channels, genes). |

## Unit model & patient list
The unit is **(patient, slide)**. `build_sample_viewer.py --dump-manifest` writes `viewer_units.tsv`
(the authoritative list: `unit_id, patient, slide, n_fovs, n_cells, fovs`) from the object's
`patient_id`/`slide`/`fov_num`. **one unit per (patient, slide)** (a patient with cores on
**both** TMA slides → one viewer per slide, `<prefix><UNIT>_<SLIDE>_sample_viewer.html`; single-slide
patients → `<prefix><UNIT>_sample_viewer.html`). The same unique `patient_id` list also appears in
`result/fov_meta_qc.tsv`.

## Coordinate model (verified 2026-08-26)
- world µm = `global_px × UMPX (0.12028)`; `x_global_px` = FOV left edge, `y_global_px` = FOV top edge.
- per cell: `gx = x0 + x_FOV_px`, `gy = y0 − y_FOV_px` (CosMx local_y grows as global_y shrinks) →
  each FOV is pasted **vertically flipped**. Confirmed: `x_slide_mm×1000 == (x0+x_FOV_px)×UMPX`,
  `y_slide_mm×1000 == (y0−y_FOV_px)×UMPX` (corr 1.0).
- `DATA.x = x_slide_mm×1000`, `DATA.y = y_slide_mm×1000`;
  `rawExtent = [xmin·UMPX, (ymin0−FOVW)·UMPX, (xmax0+FOVW)·UMPX, ymax0·UMPX]` (µm, unchanged by the flip).

## ⚠ The Y-flip fix (why `stitch_one.py` does a final flipud)
The viewer JS flips Y for the **cell scatter** when `DATA.global` is true, but **not** for the
`drawImage` morphology layer. So the stitched mosaic must be flipped top-to-bottom **once more**
(`FLIP_TOP_BOTTOM`) before base64-embedding — `stitch_one.py` does this (`flipud_applied: true`),
`rawExtent` kept in original µm. Blob-like cores can *look* aligned even when mirrored, so QC is
**quantitative**: per-cell `Mean.PanCK` vs local mosaic green must be **positive** and epithelial
cells must sit on brighter tissue. (as-is mosaic = −0.20 wrong; flipud = +0.23 correct.)

Verification uses the viewer-exact mapping (do **not** verify against the mosaic's own frame):
```
a=min(DATA.y); b=max(DATA.y); e=rawExtent
col = (x-e[0])/(e[2]-e[0]) * W
row = ((a+b-y)-e[1])/(e[3]-e[1]) * H     # includes the Y-flip
```

## Test status
One unit built end-to-end (HPC stitch + local assembly):
DATA parses, `rawImage` present, 6.02 MB, `if(DATA.rawImage)` truthy, and
`verify_overlay.py` → **corr(PanCK,green)=+0.25, epi/non-epi brightness=1.17, 94% of cells on tissue,
VERDICT PASS**.

## Notes / robustness
- Units with 0 valid FOVs are skipped and logged (`missing_fovs` recorded in each `rawimg_<unit>.json`).
- If a `rawimg_<unit>.json` is absent, the viewer still builds and the JS gate falls back to the
  per-cell virtual-H&E panel (`__rgbpath__`) — no crash.
- Total HPC cost for the full array: 27 tasks × ~15–60 s ≈ a few CPU-minutes (2 CPU / 8 GB / 20 min cap each).
- Paths (object, labels, HPC build dir, envs) are `DEFAULTS` at the top of each script — edit there to relocate.
- Outputs land in `viewers/`; nothing is deployed. Copy to the study `viewer/` dir only after review.

## CORRECTNESS FIX (2026-08-26) — authoritative patient assignment
The object's `obs.patient_id` is OVER-ASSIGNED by the `tma_location` mapping
(a unit can be credited with an order of magnitude more FOVs and cells than it
cells, TMA_1-only). `build_sample_viewer.py` now derives the per-cell patient
from the AUTHORITATIVE `fov_meta_qc.tsv` (slide,fov)->patient join
(`H5adLite.patient_auth`), used for BOTH the manifest and the per-unit cell
selection; `stitch_one.py` inherits the corrected FOV lists via the manifest.
Validation: unit / patient / cell totals match `cohort_config.json` `expected`, and every
per-patient count matches `oncore_cells_per_patient` in cohort_cellcounts.json
(0 mismatches). The spurious extra per-slide unit is gone.
