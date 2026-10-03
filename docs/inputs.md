# Input formats (`build_viewer.py`)

All tables are tab-separated (`.tsv`, optionally gzipped). Column names are case-sensitive. Extra columns are allowed.

## Required

### `--cells` (one row per cell)
| column | meaning |
|---|---|
| `cell` | unique cell id; must match `cells.txt` of the count export if `--counts-dir` is given |
| `slide` | slide label (FOV ids are per slide) |
| `fov` | FOV id (integer) |
| `core` | core id (must appear in `--cores`) |
| `x_local_px`, `y_local_px` | cell centroid in FOV-local pixels (CosMx `CenterX_local_px`, `CenterY_local_px`) |
| `L1`, `L2`, `L3` | three-level cell-type labels (L2 groups L3 in the legend) |
| any channel column | listed in config `channels` (display name → column), e.g. mean IF intensities |

### `--cores` (one row per core)
| column | meaning |
|---|---|
| `unit` | viewer unit the core belongs to (one viewer per unit; one slide per unit) |
| `core` | core id |
| `slide` | slide label |
| `group` | optional core group key (frame colour, see config `core_groups`) |
| `fovs` | optional comma-separated FOV list; FOVs listed here are imaged even if they have no cells |
| other columns | annotation fields referenced by config (`core_short_fields`, `core_tooltip_fields`, `core_summary_fields`) |

### `--units` (one row per unit)
`unit` plus any fields referenced by config `patient_fields` (shown as header chips).

### `--fov-positions`
`slide`, `fov`, `x_global_px`, `y_global_px` (CosMx `fov_positions_file`: FOV origin in global pixels).

## Optional

| option | content |
|---|---|
| `--counts-dir` | CosMx-style export: `counts.mtx` (MatrixMarket, genes × cells), `genes.txt`, `cells.txt`. Stored per unit as normalised (1e4) log1p values quantised to 8 bit per gene. |
| `--fov-images` | `<slide>_F<fov:05d>.jpg` (or .png), one square image per FOV, row 0 = local y 0. Any size; ~256 px gives the default 0.5 px/µm overview. Build with `make_fov_tiles.py --tile-px 256`. |
| `--fov-images-hires` | same naming, larger tiles (2048 px ≈ 0.25 µm/px) for the zoom-in layer; units without tiles simply have no full-resolution layer. |
| `--he-images` | `<core>.jpg` per core; shown fitted to the core frame in the H&E panel (display alignment, not registration). |
| `--regions` | `unit`, `core`, `id`, `tier` (1–3), `cx_um`, `cy_um`, `r_um` in **native slide microns** (same frame as the cells: `X = px·(x_global + x_local)`, `Y = px·(−y_global + y_local)`), plus any fields named in config `region.info_fields`. |

## Index flags (`build_index.py --flags`)
`unit`, `core`, `label` (badge text), `detail` (note under the core), `tooltip` (badge hover text).

## Config
See `input/cohort_config.example.json`. Every key is optional; defaults live in `DEFAULTS` at the top of
`build_viewer.py`. Keep your filled-in `cohort_config.json` out of version control.
