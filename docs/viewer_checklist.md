# CosMx per-unit viewer — feature checklist

**Purpose.** Every time the viewer pipeline is rebuilt or a change is deployed, a reviewing agent
MUST verify each row of this checklist against the deployed set. If any row fails, the deploy is
BLOCKED until the root cause is fixed.

Run the mechanical check with:
```bash
<env>/bin/python  # e.g. a conda env named scvi \
  <data_dir>/result/viewer_pipeline/verify_viewer_features.py \
  <deploy_dir>
```
It exits non-zero on any failure and writes a per-viewer matrix. Then eyeball 2–3 viewers with
`verify_overlay.py` for cells-on-tissue.

---

## Feature checklist (each item = one row in the verifier)

### DATA integrity
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| D1 | DATA JSON parseable | `const DATA=…{…}` in HTML | brace-balanced extract → `json.loads` | template edit corrupted the block or `__DATA__` placeholder leaked in |
| D2 | Cell count matches per-unit manifest | `viewer_units.tsv` → `n_cells` | `len(DATA.x) == manifest[unit].n_cells` | wrong (patient,slide) join in `patient_auth` (over-assign bug re-emerged) |
| D3 | All parallel arrays same length | build_sample_viewer.py | `len(x)==len(y)==len(fov)==len(type)==len(subtype)==len(tumor)==len(genes[each])` | subset filter drifted between arrays |
| D4 | Expected viewer files present in the deploy directory | manifest | `len(glob(<prefix>*_sample_viewer.html)) == expected.units` | one unit failed silently in the build stage |
| D5 | Patient metadata populated | `PATIENT_FIELDS` in build | `DATA.patient.Patient` non-empty; `Outcome` non-NA | patient join broken; obs field renamed upstream |

### Cell typing / labels
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| L1 | L2/L3 labels present | `subtype/cgc_hierarchical_labels_v2.tsv.gz` | `len(unique(DATA.subtype))>=3` and no "Unknown" > 20% | labels file path stale; join key mismatch |
| L2 | `major_celltype` categorical | h5ad `obs.major_celltype` | `len(set(DATA.type))>=3` | dropped column in a schema refresh |
| L3 | Tumor mask populated | h5ad tumor_putative | `sum(DATA.tumor)>0` for tumor-carrying units (some units legitimately have none) | tumor call regressed |

### Coordinate / geometry
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| G1 | **Degap applied** (empty bands compacted) | `degap_viewer.py` output | `DATA.global == false` | postprocess stage skipped |
| G2 | **Regrid applied** (tidy grid) | `regrid_viewer.py` output | `DATA.rawExtent[0]==0 and DATA.rawExtent[1]==0` | postprocess stage skipped |
| G3 | No residual empty band ≥2000 µm (after degap+regrid) | derived from DATA.x/y | `max_gap(DATA.x, bin=50) < 2000` — inter-core regrid padding gets up to ~1500 µm; a real un-degap regression shows 5,000–11,000 µm | degap step skipped |
| G4 | rawExtent aspect matches embedded image | `DATA.rawExtent` + PIL of `DATA.rawImage` | `abs((W/H) - ((e[2]-e[0])/(e[3]-e[1])))/ratio < 0.05` | regrid re-tiled cells but not mosaic (or vice-versa) |
| G5 | Cells-on-tissue overlay corr > 0.10 | `verify_overlay.py` | corr(PanCK, green mosaic) ≥ 0.10, % cells on tissue ≥ 60% | Y-flip regressed or wrong rawExtent |

### RAW morphology track
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| M1 | `rawImage` present (data URI) | build injects from `rawimg_all/rawimg_<unit>.json` | HTML contains `"rawImage":"data:image/jpeg;base64,` | `--rawimg-dir` pointed at pipeline root (where only a stray `rawimg_<unit>.json` may sit) instead of `rawimg_all/` |
| M2 | `rawImage` magic bytes valid | inline base64 | PIL open + `verify()` succeeds | truncation / corrupt regrid crop |
| M3 | Legend "CellComposite morphology" swatch entry visible | template legend branch | HTML contains `CellComposite morphology` | template regression |

### Important-region highlight
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| R1 | Checkbox UI present on all viewers | template | HTML contains `id="regionchk"` and `Important region` | template regressed |
| R2 | `drawRegion()` JS present on all viewers | template | HTML contains `function drawRegion(p)` and `showRegion` | template regressed |
| R3 | `DATA.importantRegion` key exists on ALL viewers (null or object) | build + finalize | HTML contains `"importantRegion":null` OR `"importantRegion":{` | finalize step skipped |
| R4 | Units listed in `cohort_config.json` `expected.region_units` have `importantRegion` **object** with valid polygon | the curated important region table | each configured unit → `DATA.importantRegion` is a dict with `poly`(≥ `min_pts`), `cx`, `cy`, `pdx` non-null | `important_regions` in `cohort_config.json` missing rows OR `finalize_region_polygon.py` skipped |
| R5 | Non-configured units have `importantRegion:null` (checkbox greys out) | template handler | remaining 24 viewers → `"importantRegion":null` | postprocess forgot to seed the key |
| R6 | Cell count inside important region matches the curated table | curated region table `n` | within 25% of `expected.region_units[<unit>].approx_cells` | `important_regions` coord drift vs the curated table |

### Interactive UI
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| U1 | Size + opacity sliders (`#psz`, `#pop`) | template | `id="psz"` and `id="pop"` present, both have `oninput` handlers | header block edit dropped controls |
| U2 | Reset-view button | template | `id="reset"` present with `onclick` binding `fit()` | header edit |
| U3 | Feature combo box (both panels) | template | `id="ci0"` `id="ci1"` `id="cp0"` `id="cp1"` present | panel block edit |
| U4 | Legend cell-type click-toggle (highlight/hide) | template (post-legend-toggle version) | `.li[data-cat]` selector + `p.state` map in JS | reverted to `.bak_prelegendtoggle` accidentally |
| U5 | Right-panel default = RAW morphology when available | template `panels[]` init | `panels=[{feat:"__type__"},{feat:DATA.rawImage?"__rawimg__":"__rgbpath__"}]` present | template edit changed defaults |

### Manifest / cohort integrity
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| C1 | Authoritative patient join (no over-assign) | `fov_meta_qc.tsv` join in `H5adLite.patient_auth` | sum(DATA.n across all viewers) == `expected.cell_total`; no unit splits per slide unless it is listed in `expected.multi_slide_units` | build fell back to `obs.patient_id` (fov_meta_qc missing or path broken) |
| C2 | Multi-slide units split per slide | manifest | one `<unit>_<slide>` file per slide for every unit in `expected.multi_slide_units`, and no others | manifest split logic regressed |
| C3 | Every unit has a valid outcome | h5ad `patient_outcome` | `DATA.patient.Outcome` ∈ `expected.outcome_allowed` | clinical join broken |

### Deploy hygiene
| # | Feature | Source of truth | Test | Fails when… |
|---|---|---|---|---|
| P1 | Pre-deploy backup exists | `stage_deploy` in `run_all.sh` | Directory `_bak_<TS>_pre_deploy/` under the deploy directory viewer/ with 26 files | deploy was manual `cp` bypassing `run_all.sh` |
| P2 | File sizes reasonable | expected 0.4–9 MB post-degap-regrid | all `.size > 100 KB` and `< 20 MB` | 09-08 regression: files <500KB and morphology missing OR file grew past 20 MB from unnecessary duplication |
| P3 | Index page linked | the configured `naming.index_file` in the same dir | file exists and links to every unit | index rebuild forgotten |

---

## Fast triage (when a check fails)

| Check that failed | Likely cause | Fix |
|---|---|---|
| M1 (morphology) | `--rawimg-dir` wrong | rerun `./run_all.sh viewers postprocess deploy` — `run_all.sh` uses `$PIPE/rawimg_all` |
| G1, G2 (degap/regrid) | postprocess stage skipped | rerun `./run_all.sh postprocess deploy` |
| T3, T5 (importantRegion key missing) | finalize step skipped | rerun `./run_all.sh postprocess deploy` (finalize is inside postprocess) |
| T4 (region polygon missing/wrong count) | `important_regions` in `cohort_config.json` drifted from the curated region table | re-derive the coords, update the config, rerun `viewers postprocess deploy` |
| C1 (over-assign; a unit gains an unexpected per-slide split) | `fov_meta_qc.tsv` path stale | check `DEFAULTS['fov_meta_qc']` in `build_sample_viewer.py`, then rerun |
| U1–U5 (UI features) | template edited | diff against `sample_viewer_template.html.bak_pretls_20260908_132737` |

---

## Never-do list (regressions that have happened before)

1. **Never rebuild with `--rawimg-dir $PIPE`** — only a stray `rawimg_<unit>.json` may sit at the pipeline root; the full set lives in `rawimg_all/`. Silent morphology-loss for 25/26 units.
2. **Never copy `viewers/*.html` to the deploy directory without running `postprocess` first.** That reintroduces the 09-04 and 09-08 regressions (no degap, no regrid, no important region).
3. **Never inject important-region overlays via ad-hoc post-scripts to individual viewers.** The 09-03 build did that and the code was lost. All important-region overlays are now driven by `important_regions` in `cohort_config.json`.
4. **Never overwrite `sample_viewer_template.html` without saving a `.bak_<what>_<ts>` first.**
5. **Never trust `obs.patient_id` for the per-unit join.** It's over-assigned by the `tma_location` map. Always go through `H5adLite.patient_auth(fov_meta_qc)`.


---

## 12 September 2026 verification run (baseline PASS)

Command:
```bash
<env>/bin/python verify_viewer_features.py \
  <deploy_dir> --json <out>.json
```

Result: **all 26 per-viewer checks PASS × all 7 cohort checks PASS** (0 failures).
See `verify_report_20260908_134738.json`. This is the reference baseline for future
regressions — any deviation should trace to the "Never-do list" or the "Fast triage" table above.
