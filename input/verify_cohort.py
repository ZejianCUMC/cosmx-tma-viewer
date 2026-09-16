#!/usr/bin/env python
"""Per-viewer verification of the core-annotated cohort: DATA integrity, frame
geometry (no dot inside/outside the wrong frame, no frame or slot overlap),
core->clinical coverage, and morphology preservation vs the deployed set."""
import argparse, base64, glob, io, json, os, re, sys
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

def _cohort_cfg():
    import json as _json
    for _d in (os.path.dirname(os.path.abspath(__file__)),
               os.path.dirname(os.path.dirname(os.path.abspath(__file__))), os.getcwd()):
        _p = os.path.join(_d, "cohort_config.json")
        if os.path.exists(_p):
            try:
                return _json.load(open(_p, encoding="utf-8"))
            except Exception:
                return {}
    return {}


_PREFIX = ((_cohort_cfg().get("naming") or {}).get("file_prefix") or "sample_")

# Optional reference set to compare morphology resolution against
# (e.g. the currently deployed viewers). Pass --reference-dir to enable.
REFERENCE_DIR = None

def data_of(p):
    t = io.open(p, encoding="utf-8").read()
    m = re.search(r"const DATA\s*=\s*", t); s = m.end(); d = 0; i = s
    while i < len(t):
        c = t[i]
        if c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: i += 1; break
        i += 1
    return json.loads(t[s:i]), t

def img_dims(u):
    raw = base64.b64decode(u.split(",", 1)[1])
    return Image.open(io.BytesIO(raw)).size, len(raw)

rows, bad = [], 0
ap = argparse.ArgumentParser(description="Verify a core-annotated viewer set.")
ap.add_argument("viewer_dir")
ap.add_argument("--prefix", default=_PREFIX)
ap.add_argument("--reference-dir", default=None,
                help="compare morphology resolution against this set (e.g. the deployed one)")
args = ap.parse_args()
REFERENCE_DIR = args.reference_dir

for p in sorted(glob.glob(os.path.join(args.viewer_dir, args.prefix + "*_sample_viewer.html"))):
    unit = os.path.basename(p)[len(args.prefix):-len("_sample_viewer.html")]
    D, txt = data_of(p)
    x = np.array(D["x"]); y = np.array(D["y"]); fov = np.array(D["fov"])
    ci = np.array([D["coreFov"][str(f)] for f in fov])
    errs = []
    # every core's cells inside its own frame, no foreign cells, no overlaps
    clear = []
    for k, c in enumerate(D["cores"]):
        b = c["box"]; own = ci == k
        inx = (x >= b[0]) & (x <= b[2]) & (y >= b[1]) & (y <= b[3])
        if int((~inx & own).sum()): errs.append(f"{c['id']}:cells outside own frame")
        if int((inx & ~own).sum()): errs.append(f"{c['id']}:foreign cells in frame")
        clear.append(min(x[own].min()-b[0], b[2]-x[own].max(), y[own].min()-b[1], b[3]-y[own].max()))
        if not c.get("clin"): errs.append(f"{c['id']}:no clinical string")
    for a in range(len(D["cores"])):
        for b2 in range(a+1, len(D["cores"])):
            A, B = D["cores"][a]["box"], D["cores"][b2]["box"]
            if min(A[2],B[2])-max(A[0],B[0]) > 0 and min(A[3],B[3])-max(A[1],B[1]) > 0:
                errs.append(f"frame overlap {D['cores'][a]['key']}/{D['cores'][b2]['key']}")
            S, T = D["cores"][a]["slot"], D["cores"][b2]["slot"]
            if min(S[2],T[2])-max(S[0],T[0]) > 1 and min(S[3],T[3])-max(S[1],T[1]) > 1:
                errs.append(f"slot overlap {D['cores'][a]['key']}/{D['cores'][b2]['key']}")
    # header contract + UI markers
    if "Grade" in D["patient"] or "Stage" in D["patient"] or "Timepoint" in D["patient"]:
        errs.append("core-level field still in the patient chip row")
    if "coresUnmapped" in D: errs.append("cores with no cells still shipped")
    for lit in ['id="psz"', 'id="ci0"', '.li[data-cat]', 'function drawRegion(p)', 'id="regionchk"',
                'CellComposite morphology', 'DATA.rawImage?"__rawimg__"', 'renderCoreTags']:
        if lit not in txt: errs.append(f"missing UI literal {lit}")
    # Morphology must not lose RESOLUTION vs the deployed viewer. Compare um/px,
    # not pixel count: the old regrid always laid out ncol=3 columns even for a
    # single core, so its canvas carried 2 columns of black padding. A tighter
    # canvas at the same density is an improvement, not a regression.
    dep = os.path.join(REFERENCE_DIR, os.path.basename(p)) if REFERENCE_DIR else ""
    mor = "n/a"
    if D.get("rawImage") and dep and os.path.exists(dep):
        (nw, nh), nb = img_dims(D["rawImage"])
        ne = D["rawExtent"]; n_umpx = (ne[2]-ne[0])/nw
        dD, _ = data_of(dep)
        if dD.get("rawImage"):
            (dw, dh), db = img_dims(dD["rawImage"])
            de = dD["rawExtent"]; d_umpx = (de[2]-de[0])/dw
            mor = f"{n_umpx:.2f} um/px (dep {d_umpx:.2f})"
            if n_umpx > d_umpx * 1.05:
                errs.append(f"MORPHOLOGY COARSER: {n_umpx:.3f} um/px vs deployed {d_umpx:.3f}")
    ncore = len({c["id"] for c in D["cores"]})
    rows.append((unit, D["n"], D["fovs"], ncore, len(D["cores"]), min(clear), os.path.getsize(p)/1e6, mor, errs))
    bad += len(errs)

print(f"{'unit':14s} {'cells':>7s} {'FOV':>4s} {'cores':>5s} {'frames':>6s} {'clear':>6s} {'MB':>6s}  morphology")
tot_c = tot_f = 0
for u, n, f, nc, nf, cl, mb, mor, errs in rows:
    tot_c += n; tot_f += f
    print(f"{u:14s} {n:7d} {f:4d} {nc:5d} {nf:6d} {cl:6.1f} {mb:6.1f}  {mor}")
    for e in errs: print(f"    !! {e}")
print(f"\n{len(rows)} viewers | {tot_c:,} cells | {tot_f} FOVs | "
      f"{sum(r[3] for r in rows)} cores -> {sum(r[4] for r in rows)} frames")
print(f"VERDICT: {'PASS - no geometry, clinical or UI errors' if bad==0 else f'FAIL - {bad} errors'}")
sys.exit(1 if bad else 0)
