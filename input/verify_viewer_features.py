#!/usr/bin/env python
"""
verify_viewer_features.py
=========================
Mechanical checker for the per-unit viewer feature checklist.
Runs every row of viewer_checklist.md against every deployed viewer HTML in a
directory and prints a pass/fail matrix. Exits non-zero on ANY failure so it can
gate deploys in CI/agents. See viewer_checklist.md for rationale.

Usage:
    verify_viewer_features.py <viewer_dir>                # e.g. the deploy directory share
    verify_viewer_features.py <viewer_dir> --json out.json # machine-readable report
"""
import argparse, base64, glob, io, json, os, re, sys

# -- cohort config (study data lives OUTSIDE the code) ------------------------
# Unit ids, per-unit pathology labels and cohort counts are study data, not
# code, so they are read from cohort_config.json instead of being hardcoded.
# See cohort_config.example.json. A missing file or key disables the feature
# that uses it; nothing else changes.
def _cohort_cfg():
    import json as _json
    for _d in (os.path.dirname(os.path.abspath(__file__)),
               os.path.dirname(os.path.dirname(os.path.abspath(__file__))), os.getcwd()):
        _p = os.path.join(_d, "cohort_config.json")
        if os.path.exists(_p):
            try:
                return _json.load(open(_p, encoding="utf-8"))
            except Exception as _e:
                print(f"[warn] bad cohort_config.json ({_e}); continuing without it")
                return {}
    return {}

_COHORT = _cohort_cfg()
_PREFIX = ((_COHORT.get("naming") or {}).get("file_prefix") or "sample_")
_INDEX = ((_COHORT.get("naming") or {}).get("index_file") or "index_all_viewers.html")
from collections import defaultdict
try:
    import numpy as np
except ImportError:
    np = None
try:
    from PIL import Image
    _HAVE_PIL = True
except ImportError:
    _HAVE_PIL = False

# ---- expected constants (source of truth) ---------------------------------
_E = (_COHORT.get("expected") or {})
EXPECTED_UNITS = _E.get("units")                       # None -> skip the count check
EXPECTED_CELL_TOTAL = _E.get("cell_total")             # None -> skip the total check
MULTI_SLIDE_UNITS = set(_E.get("multi_slide_units") or ())
REGION_UNITS = {k: v for k, v in (_E.get("region_units") or {}).items() if not k.startswith("_")}
OUTCOME_ALLOWED = set(_E.get("outcome_allowed") or ())  # empty -> skip the outcome check
GENE_MIN = 30
SIZE_MIN_MB, SIZE_MAX_MB = 0.1, 20.0
PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))

# ---- DATA extractor -------------------------------------------------------
def extract_data(txt):
    m = re.search(r"const DATA\s*=\s*", txt)
    if not m: return None
    s = m.end(); d = 0; i = s
    while i < len(txt):
        c = txt[i]
        if c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: i += 1; break
        i += 1
    try: return json.loads(txt[s:i])
    except Exception: return None

def unit_from_filename(fn):
    base = os.path.basename(fn)
    m = re.match(re.escape(_PREFIX) + r"(.+)_sample_viewer\.html$", base)
    return m.group(1) if m else None

# ---- per-viewer feature tests ---------------------------------------------
def scan_viewer(path):
    """Return dict of {check_id: (bool_pass, msg)} for one viewer."""
    r = {}
    txt = open(path).read()
    size_mb = len(txt)/1e6
    unit = unit_from_filename(path)
    patient = unit.split("_")[0] if unit else None

    D = extract_data(txt)
    r["D1_data_parseable"] = (D is not None, "extract_data returned None" if D is None else f"n={D.get('n')}")
    if D is None: return r, size_mb, unit, patient

    # ---- DATA integrity ----------------------------------------------------
    n = int(D.get("n", 0))
    arrays = ["x","y","fov","type","subtype","tumor"]
    lens = {k: len(D.get(k,[])) for k in arrays}
    same_len = all(v == n for v in lens.values())
    r["D3_array_lengths_match"] = (same_len, "OK" if same_len else f"n={n} but {lens}")

    patient_field = (D.get("patient") or {}).get("Patient","")
    r["D5_patient_metadata"] = (bool(patient_field and patient_field != "NA"),
                                f"Patient='{patient_field}'")

    # ---- Cell typing labels ------------------------------------------------
    subtype_vals = set(D.get("subtype", []))
    r["L1_subtype_labels"] = (len(subtype_vals) >= 3,
                              f"n_unique={len(subtype_vals)} sample={list(subtype_vals)[:4]}")
    unknown = sum(1 for v in D.get("subtype",[]) if str(v).lower()=="unknown")
    r["L1b_unknown_pct"] = (unknown < 0.20*n if n else True,
                             f"unknown={unknown}/{n} ({100*unknown/max(n,1):.1f}%)")
    r["L2_major_type"] = (len(set(D.get("type",[]))) >= 3,
                          f"n_unique={len(set(D.get('type',[])))}")

    # ---- Coordinate / geometry --------------------------------------------
    gl = D.get("global")
    r["G1_degap_applied"] = (gl is False, f"DATA.global={gl!r} (expect False)")
    ext = D.get("rawExtent")
    if ext and len(ext)==4:
        r["G2_regrid_applied"] = (float(ext[0])==0.0 and float(ext[1])==0.0,
                                   f"rawExtent starts=({ext[0]:.1f},{ext[1]:.1f}) expect (0,0)")
    else:
        r["G2_regrid_applied"] = (False, f"rawExtent missing/malformed: {ext}")

    # gap detection (post-degap): no empty band >500µm in x or y
    if np is not None and n > 0:
        def biggest_empty_band(arr, binsize=50):
            arr = np.asarray(arr, dtype=float)
            if len(arr) < 5: return 0.0
            lo, hi = float(arr.min()), float(arr.max())
            edges = np.arange(lo, hi+binsize, binsize)
            h, _ = np.histogram(arr, bins=edges)
            empty = h == 0
            best = 0; i = 0
            while i < len(empty):
                if empty[i]:
                    j = i
                    while j < len(empty) and empty[j]: j += 1
                    best = max(best, (j-i)*binsize)
                    i = j
                else: i += 1
            return best
        gx = biggest_empty_band(D.get("x",[]))
        gy = biggest_empty_band(D.get("y",[]))
        r["G3_no_residual_gaps"] = (max(gx,gy) < 2000,
                                     f"max_x_gap={gx:.0f}µm max_y_gap={gy:.0f}µm (expect <2000; larger = degap skipped)")
    else:
        r["G3_no_residual_gaps"] = (True, "skipped (numpy missing or empty)")

    # ---- Morphology --------------------------------------------------------
    has_raw = ('"rawImage":"data:image' in txt) and (D.get("rawImage") is not None)
    r["M1_morphology_present"] = (has_raw, "OK" if has_raw else "rawImage missing")

    # rawExtent aspect vs image aspect
    if has_raw and _HAVE_PIL and ext and len(ext)==4:
        try:
            ru = D["rawImage"].split(",",1)[1]
            img = Image.open(io.BytesIO(base64.b64decode(ru)))
            img.verify()
            img = Image.open(io.BytesIO(base64.b64decode(ru)))
            W, H = img.size
            r["M2_image_bytes_ok"] = (True, f"{W}x{H}")
            iw, ih = ext[2]-ext[0], ext[3]-ext[1]
            ratio_img = W/max(1,H); ratio_ext = iw/max(1,ih)
            drift = abs(ratio_img-ratio_ext)/max(ratio_img,ratio_ext)
            r["G4_aspect_match"] = (drift < 0.08,
                                     f"img={ratio_img:.2f} ext={ratio_ext:.2f} drift={100*drift:.1f}%")
        except Exception as e:
            r["M2_image_bytes_ok"] = (False, f"PIL failed: {type(e).__name__}")
            r["G4_aspect_match"] = (False, "skipped (image parse failed)")
    else:
        r["M2_image_bytes_ok"] = (has_raw, "skipped (no PIL or no morphology)")
        r["G4_aspect_match"] = (True, "skipped (no PIL or no morphology)")

    r["M3_morphology_legend"] = ("CellComposite morphology" in txt,
                                  "legend swatch text present" if "CellComposite morphology" in txt else "swatch text missing")

    # ---- Important-region highlight -----------------------------------------------------
    r["R1_region_checkbox_ui"] = ('id="regionchk"' in txt and 'Important region' in txt,
                                "UI present" if 'id="regionchk"' in txt else "checkbox missing")
    r["T2_drawRegion_js"] = ("function drawRegion(p)" in txt and "showRegion" in txt,
                          "OK" if "function drawRegion(p)" in txt else "JS missing")
    region = D.get("importantRegion", "MISSING")
    r["R3_region_key_present"] = (region is None or isinstance(region, dict),
                                f"importantRegion={type(region).__name__ if region is not None else 'null'}")

    expected_region = patient in REGION_UNITS and unit == patient  # configured units are single-slide
    if expected_region:
        region_meta = REGION_UNITS[patient]
        ok = isinstance(region, dict) and region.get("poly") and region.get("pdx") and \
             len(region["poly"]) >= region_meta["min_pts"]
        r["R4_region_polygon"] = (ok, (f"pdx={region.get('pdx')} poly_pts={len(region.get('poly',[]))} "
                                        f"expected_pdx={region_meta['pdx']}") if isinstance(region,dict)
                                        else f"importantRegion is {type(region).__name__}")
        # cell-count sanity
        nin = int(region.get("n_cells_inside", 0)) if isinstance(region,dict) else 0
        rel = abs(nin - region_meta["approx_cells"])/region_meta["approx_cells"]
        r["R6_region_cell_count"] = (rel < 0.25,
                                   f"n_inside={nin} expected≈{region_meta['approx_cells']} drift={100*rel:.1f}%")
    else:
        r["R4_region_polygon"] = (True, "skipped (not a important-region unit)")
        r["R5_region_null_when_unset"] = (region is None or (isinstance(region,dict) and not region.get("poly")),
                                        f"importantRegion={region}")
        r["R6_region_cell_count"] = (True, "skipped (not a important-region unit)")

    # ---- Interactive UI ----------------------------------------------------
    r["U1_size_opacity_sliders"] = ('id="psz"' in txt and 'id="pop"' in txt,
                                     "OK" if 'id="psz"' in txt and 'id="pop"' in txt else "sliders missing")
    r["U2_reset_button"]   = ('id="reset"' in txt, "OK" if 'id="reset"' in txt else "reset missing")
    r["U3_combo_boxes"]    = (all(f'id="ci{i}"' in txt for i in (0,1)) and all(f'id="cp{i}"' in txt for i in (0,1)),
                              "OK" if 'id="ci0"' in txt and 'id="ci1"' in txt else "combos missing")
    r["U4_legend_toggle"]  = ('.li[data-cat]' in txt and 'p.state' in txt,
                              "OK" if '.li[data-cat]' in txt else "legend-toggle CSS/JS missing")
    r["U5_rawimg_default_right"] = ('DATA.rawImage?"__rawimg__"' in txt,
                                    "right panel defaults to morphology" if 'DATA.rawImage?"__rawimg__"' in txt else "default init missing")

    # ---- Cohort integrity --------------------------------------------------
    _oc = (D.get("patient") or {}).get("Outcome", "NA")
    r["C3_outcome_valid"] = (not OUTCOME_ALLOWED or _oc in OUTCOME_ALLOWED,
                             f"Outcome={_oc!r}" if OUTCOME_ALLOWED else "skipped (no outcome_allowed configured)")

    # ---- Deploy hygiene ----------------------------------------------------
    r["P2_size_reasonable"] = (SIZE_MIN_MB <= size_mb <= SIZE_MAX_MB,
                                f"{size_mb:.2f} MB (expect {SIZE_MIN_MB}-{SIZE_MAX_MB})")
    r["_meta_n_cells"] = (True, n)

    return r, size_mb, unit, patient


# ---- cohort-level tests ----------------------------------------------------
def cohort_checks(viewer_dir, per_viewer):
    """Checks that need the full 26-viewer set."""
    checks = {}
    units_present = sorted([u for _,_,u,_ in per_viewer if u])
    checks["D4_expected_viewers"] = (EXPECTED_UNITS is None or len(units_present) == EXPECTED_UNITS,
                                f"found {len(units_present)}/{EXPECTED_UNITS}" if EXPECTED_UNITS is not None
                                else "skipped (no expected.units configured)")
    # C2: multi-slide split
    multi_files = [u for u in units_present if u.endswith("_TMA_1") or u.endswith("_TMA_2")]
    multi_patients = set(u.split("_TMA_")[0] for u in multi_files)
    checks["C2_multi_slide_split"] = (not MULTI_SLIDE_UNITS or multi_patients == MULTI_SLIDE_UNITS,
                                       f"multi-slide units found={sorted(multi_patients)} expect={sorted(MULTI_SLIDE_UNITS)}"
                                       if MULTI_SLIDE_UNITS else "skipped (no multi_slide_units configured)")
    # No spurious <unit>
    _spurious = sorted({u for u in units_present if "_TMA_" in u
                        and u.split("_TMA_")[0] not in MULTI_SLIDE_UNITS} ) if MULTI_SLIDE_UNITS else []
    checks["C1_no_over_assign"] = (not _spurious,
                                    "no unexpected per-slide split (would signal an over-assigned FOV join)"
                                    if not _spurious else f"REGRESSION: unexpected split units {_spurious}")
    total = sum(v for r,_,_,_ in per_viewer for k,(_,v) in r.items() if k == "_meta_n_cells")
    checks["C1b_total_cells"] = (EXPECTED_CELL_TOTAL is None or total == EXPECTED_CELL_TOTAL,
                                  f"sum(n)={total} expect={EXPECTED_CELL_TOTAL}" if EXPECTED_CELL_TOTAL is not None
                                  else f"sum(n)={total} (no expected.cell_total configured)")
    # P1: pre-deploy backup exists
    bak_dirs = sorted(glob.glob(os.path.join(viewer_dir, "_bak_*_pre_*")))
    checks["P1_backup_present"] = (len(bak_dirs) >= 1,
                                    f"backup dirs found: {[os.path.basename(d) for d in bak_dirs]}")
    # P3: index page
    idx = os.path.join(viewer_dir, _INDEX)
    checks["P3_index_page"] = (os.path.exists(idx),
                                f"{'present' if os.path.exists(idx) else 'MISSING'} at {idx}")
    # important region coverage
    region_found = {u for _,_,u,p in per_viewer if p in REGION_UNITS and u == p}
    checks["R7_region_units_present"] = (not REGION_UNITS or region_found == set(REGION_UNITS.keys()),
                                       f"important-region units present={sorted(region_found)} expect={sorted(REGION_UNITS.keys())}"
                                       if REGION_UNITS else "skipped (no region_units configured)")
    return checks


# ---- main -----------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument("viewer_dir", help="dir containing <prefix>*_sample_viewer.html")
    ap.add_argument("--json", help="write machine-readable report")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a.viewer_dir, _PREFIX + "*_sample_viewer.html")))
    if not files:
        print(f"[FAIL] no {_PREFIX}*_sample_viewer.html in {a.viewer_dir}")
        sys.exit(2)

    per_viewer = [scan_viewer(f) for f in files]
    cohort = cohort_checks(a.viewer_dir, per_viewer)

    # collect check ids in a consistent order
    all_ids = []
    seen = set()
    for r,_,_,_ in per_viewer:
        for k in r:
            if k.startswith("_"): continue
            if k not in seen: seen.add(k); all_ids.append(k)

    n_pass = defaultdict(int); n_fail = defaultdict(int); fail_examples = defaultdict(list)
    for f,(r,sz,u,p) in zip(files, per_viewer):
        for k in all_ids:
            ok,msg = r.get(k,(True,"n/a"))
            (n_pass if ok else n_fail)[k] += 1
            if not ok and len(fail_examples[k]) < 3:
                fail_examples[k].append(f"{u}: {msg}")

    # ---- print per-check summary ----
    total_v = len(files)
    print(f"\n=== per-viewer checks (n={total_v} viewers) ===")
    print(f"{'check':32s}  pass / fail  first failure")
    for k in all_ids:
        p, f_ = n_pass[k], n_fail[k]
        mark = "✅" if f_ == 0 else "❌"
        first = fail_examples[k][0] if fail_examples[k] else ""
        print(f"  {mark} {k:29s}  {p:2d} / {f_:2d}  {first}")

    print(f"\n=== cohort checks ===")
    for k,(ok,msg) in cohort.items():
        print(f"  {'✅' if ok else '❌'} {k:29s}  {msg}")

    # verbose per-viewer detail on any failure
    fails = [(f,r,u,p) for f,(r,_,u,p) in zip(files,per_viewer) if any(not v[0] for k,v in r.items() if not k.startswith("_"))]
    if fails and a.verbose:
        print(f"\n=== per-viewer failure detail ({len(fails)} viewers) ===")
        for f,r,u,p in fails:
            bad = [(k,v[1]) for k,v in r.items() if not v[0] and not k.startswith("_")]
            print(f"\n  {u} ({p}):")
            for k,m in bad: print(f"    ✗ {k}: {m}")

    # overall verdict
    fail_check = any(f_ > 0 for f_ in n_fail.values()) or any(not ok for ok,_ in cohort.values())
    print(f"\n=== OVERALL: {'FAIL' if fail_check else 'PASS'} ({sum(n_fail.values())} per-viewer failures, {sum(1 for ok,_ in cohort.values() if not ok)} cohort failures) ===")

    if a.json:
        rep = {"viewer_dir": a.viewer_dir, "n_viewers": total_v,
               "per_check": {k: {"pass": n_pass[k], "fail": n_fail[k], "examples": fail_examples[k]} for k in all_ids},
               "cohort": {k: {"pass": bool(ok), "msg": msg} for k,(ok,msg) in cohort.items()},
               "overall_pass": not fail_check}
        json.dump(rep, open(a.json,"w"), indent=2)
        print(f"[json] written {a.json}")

    sys.exit(1 if fail_check else 0)

if __name__ == "__main__":
    main()
