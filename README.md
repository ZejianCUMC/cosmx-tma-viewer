# CosMx sample viewer

Scripts that turn a CosMx spatial-transcriptomics AnnData object into a
self-contained, offline HTML viewer — one file per (patient, slide) unit, with
no server and no network dependency.

Each built viewer carries, inline:

- every cell as a point in slide microns, coloured by major type (L1), subtype
  (L2 / L3), putative-tumour flag, any of ~6,000 panel genes, or an IF channel;
- the real CosMx `CellComposite` morphology mosaic as a synced second panel;
- two linked pan/zoom canvases, a searchable feature picker, click-to-recolour
  legends, and per-cell hover;
- a frame per TMA core carrying that core's own pathology and treatment
  timepoint (see `docs/core_annotation.md`).

**This repository is scripts only.** Built viewers embed per-patient clinical
metadata and run 5–80 MB each; they, and every source matrix / label table /
clinical workbook, are excluded by `.gitignore`. Keep it that way.

## Layout

| Path | What it is |
|---|---|
| `input/` | Every script, plus `cohort_config.example.json` and the viewer shell `sample_viewer_template.html`. |
| `docs/` | Feature checklist (the source of truth the verifier tests against), per-core annotation notes, morphology-resolution recipe, pipeline notes. |

Inside `input/`, the scripts fall into four groups:

| Group | Files |
|---|---|
| Generator | `build_sample_viewer.py`, `sample_viewer_template.html` |
| Geometry post-processing | `degap_viewer.py`, `regrid_viewer.py`, `regrid_by_core.py`, `assign_subcores.py`, `finalize_region_polygon.py` |
| Per-core annotation | `inject_core_clinical.py`, `make_template_cores.py` |
| Morphology (HPC) | `stitch_one.py`, `stitch_one_hires_1umpx.py`, `stitch_one_hires_025umpx.py`, `stitch_morphology*.sbatch`, `verify_only_morphology_changed.py` |
| Verification / driver | `verify_viewer_features.py`, `verify_overlay.py`, `build_all_viewer_index.py`, `run_all.sh` |

## How a viewer is built

```
build_sample_viewer.py     object + labels + FOV-to-unit join -> DATA -> template -> HTML
  |
assign_subcores.py         authoritative core_id, then physical sub-core split   (PRE-degap)
  |
degap_viewer.py            compact empty FOV-layout bands, cells + mosaic together
  |
regrid_by_core.py          one tidy grid tile per sub-core   (or regrid_viewer.py, legacy)
  |
inject_core_clinical.py    attach per-core pathology / timepoint / grade / stage
  |
finalize_region_polygon.py convex hull for configured important regions, final coordinates
  |
verify_viewer_features.py  ~30 mechanical checks; non-zero exit gates a deploy
```

`run_all.sh` drives a cohort end to end (manifest -> HPC stitch -> pull -> build ->
post-process -> verify -> deploy). The post-processing scripts only ever rewrite
the `const DATA={...}` span, never the HTML/CSS/JS shell, so the shell stays
byte-identical to `input/sample_viewer_template.html`. That contract is what
lets a UI change be applied to already-built viewers by swapping the shell
instead of rebuilding from the expression matrix.

## Cohort configuration

Anything specific to a study — unit ids, per-unit pathology labels, curated important region
region coordinates, expected cohort counts — is **study data, not code**, and is
read at runtime from `cohort_config.json` in the repository root. That file is
gitignored. Only `cohort_config.example.json`, which contains placeholders, is
tracked.

```bash
cp cohort_config.example.json cohort_config.json   # then fill in for your cohort
```

Every key is optional: a missing file or key disables the feature or check that
uses it rather than failing. With no config at all the pipeline still builds
viewers — important-region overlays are simply absent and the cohort-integrity checks report
`skipped` instead of comparing against expected values.

## Conventions

- Every script opens with a header block giving description, inputs, outputs,
  environment, key dependencies and validation date, then a `DEFAULTS` dict
  holding every tunable value.
- Data paths in `DEFAULTS` are `<project_root>/…` placeholders, not real
  locations: pass them explicitly (`--h5ad`, `--labels`, `--fov-meta-qc`,
  `--template`, `--rawimg-dir`, `--outdir`, `--src`, …) or edit `_ROOT` once for
  your own layout. Paths that resolve relative to the script itself (template,
  output dir) work as-is.
- `verify_viewer_features.py` greps the built HTML for literal UI markers
  (`id="psz"`, `.li[data-cat]`, `function drawRegion(p)`, …). Renaming one of those
  in the template fails the gate even if the UI still works — update both.

## Requirements

Python 3.9+ with `numpy`, `pandas`, `scipy`, `h5py`, `Pillow`, and `openpyxl`
for the clinical-workbook step. The viewers themselves need only a browser.
The morphology stitchers additionally expect the CosMx `flatFiles` /
`DecodedFiles` trees and are written to run as a SLURM array, one task per unit.
