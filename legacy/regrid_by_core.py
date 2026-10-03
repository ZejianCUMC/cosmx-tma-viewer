#!/usr/bin/env python
# =============================================================================
# viewer — CORE-AWARE regrid  (replaces regrid_viewer.py's spatial guess)
# =============================================================================
# Description : Re-arrange a sample viewer's tissue pieces into a tidy N-column
#               grid, one tile per SUB-CORE. Grouping is NOT regrid_viewer.py's
#               union-find guess on FOV centroids (which has no clinical
#               grounding and can group FOVs across two different core_ids); it
#               comes from the sub-core keys assign_subcores.py wrote PRE-degap
#               = authoritative core_id, then physical connected components
#               inside it. One tile is therefore exactly one physical tissue
#               piece and never spans two core_ids; sub-cores of one punch keep
#               the shared parent core_id as their label.
#               Emits DATA.cores (frame box + grid slot per sub-core) and
#               DATA.coreFov (fov -> index into DATA.cores).
# Input       : viewer HTML (post-degap; DATA.global false; DATA.coreFov and
#               DATA.coreSub present from assign_subcores.py)
# Output      : same HTML rewritten in place (DATA span only)
# Conda env   : <env>/bin/python  # e.g. a conda env named scvi
# Key deps    : numpy, Pillow
# Validated   : 2026-09-15 (<unit>: 6 core_id / 7 sub-cores, no tile overlap)
# Parameters  : See DEFAULTS dict below
# =============================================================================

import argparse, base64, json, re
from io import BytesIO

import numpy as np

try:
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None          # hi-res mosaics trip the bomb limit
except Exception:
    Image = None

# -- DEFAULTS (all tunable parameters) ----------------------------------------
DEFAULTS = {
    "ncol":         3,      # grid columns
    "frame_pad_um": 45.0,   # clearance between the cell bbox and the frame line
    "gutter_um":    170.0,  # empty space between neighbouring grid slots
    "jpeg_quality": 84,     # mosaic re-encode quality (matches regrid_viewer.py)
}


def extract_data(txt):
    """Locate and parse the `const DATA={...}` object; returns (start, end, dict)."""
    m = re.search(r"const DATA\s*=\s*", txt)
    s = m.end(); d = 0; i = s
    while i < len(txt):
        c = txt[i]
        if c == '{':
            d += 1
        elif c == '}':
            d -= 1
            if d == 0:
                i += 1
                break
        i += 1
    return s, i, json.loads(txt[s:i])


def regrid(path, out=None, **kw):
    P = {**DEFAULTS, **{k: v for k, v in kw.items() if v is not None}}
    txt = open(path, encoding="utf-8").read()
    s, e, D = extract_data(txt)

    x = np.array(D["x"], float)
    y = np.array(D["y"], float)
    fov = np.array(D.get("fov", [0] * len(x)), int)

    if D.get("global"):                                # bake render-time y-flip
        a, b = y.min(), y.max()
        y = b - (y - a)
        D["global"] = False

    f2k = D.get("coreFov") or {}
    sub = D.get("coreSub") or {}
    if not f2k or not sub:
        raise SystemExit("[err] DATA.coreFov/coreSub missing — run assign_subcores.py "
                         "on the PRE-degap build first")
    missing = sorted({int(f) for f in fov if str(int(f)) not in f2k})
    if missing:
        raise SystemExit(f"[err] FOVs without a sub-core key: {missing[:12]}")
    core = np.array([f2k[str(int(f))] for f in fov])   # sub-core key per cell
    # tile order: parent core_id, then part number (sub-cores stay adjacent)
    ids = sorted(sub, key=lambda k: (sub[k]["core"], sub[k]["part"]))

    # -- per-core bbox in the current (degapped) space -------------------------
    fp = float(P["frame_pad_um"])
    info = {}
    for cid in ids:
        k = core == cid
        x0, x1, y0, y1 = x[k].min(), x[k].max(), y[k].min(), y[k].max()
        info[cid] = dict(x0=x0, x1=x1, y0=y0, y1=y1, cx=(x0 + x1) / 2, cy=(y0 + y1) / 2,
                         w=x1 - x0, h=y1 - y0, n=int(k.sum()),
                         fovs=sorted(int(v) for v in set(fov[k])))

    # reject overlapping source bboxes: the mosaic crop would pull in a neighbour
    for a_ in range(len(ids)):
        for b_ in range(a_ + 1, len(ids)):
            A, B = info[ids[a_]], info[ids[b_]]
            if min(A["x1"], B["x1"]) - max(A["x0"], B["x0"]) > 0 and \
               min(A["y1"], B["y1"]) - max(A["y0"], B["y0"]) > 0:
                raise SystemExit(f"[err] source bboxes overlap: {ids[a_]} x {ids[b_]}; "
                                 "mosaic crop would include a neighbouring core")

    # -- tidy grid: per-column width and per-row height (a single tall core must
    #    not inflate every row, which a uniform cell size would do) ------------
    ncol = int(P["ncol"])
    gut = float(P["gutter_um"])
    nrow = int(np.ceil(len(ids) / ncol))
    pos = {cid: (k // ncol, k % ncol) for k, cid in enumerate(ids)}
    colw = [max([info[c]["w"] + 2 * fp for c in ids if pos[c][1] == j] or [0]) for j in range(ncol)]
    rowh = [max([info[c]["h"] + 2 * fp for c in ids if pos[c][0] == i] or [0]) for i in range(nrow)]
    colx = np.concatenate([[0.0], np.cumsum([w + gut for w in colw])])
    rowy = np.concatenate([[0.0], np.cumsum([h + gut for h in rowh])])
    totW = float(colx[ncol] - gut)
    totH = float(rowy[nrow] - gut)

    # -- move cells; record frame box + slot rect in FINAL coords --------------
    nx, ny = x.copy(), y.copy()
    cores_out = []
    for cid in ids:
        r, c = pos[cid]
        v = info[cid]
        sx0, sy0 = float(colx[c]), float(rowy[r])
        tx = sx0 + colw[c] / 2                          # slot centre
        ty = sy0 + rowh[r] / 2
        k = core == cid
        nx[k] += tx - v["cx"]
        ny[k] += ty - v["cy"]
        box = [tx - v["w"] / 2 - fp, ty - v["h"] / 2 - fp,
               tx + v["w"] / 2 + fp, ty + v["h"] / 2 + fp]
        parent = sub[cid]["core"]                       # sub-cores SHARE the core_id
        cores_out.append({
            "key": cid, "id": parent, "part": sub[cid]["part"], "nparts": sub[cid]["nparts"],
            "loc": parent.split("_", 1)[1] if "_" in parent else parent,
            "n": v["n"], "nfov": len(v["fovs"]), "fovs": v["fovs"],
            "box": [round(float(t), 1) for t in box],
            "slot": [round(sx0, 1), round(sy0, 1),
                     round(sx0 + colw[c], 1), round(sy0 + rowh[r], 1)],
            "src": [v["x0"], v["y0"], v["x1"], v["y1"]],   # dropped below; mosaic only
        })

    # -- re-tile the mosaic: crop each core's padded box, paste at its new box --
    if Image is not None and D.get("rawImage"):
        e0, e1, e2, e3 = D["rawExtent"]
        im = Image.open(BytesIO(base64.b64decode(D["rawImage"].split(",", 1)[1]))).convert("RGB")
        IW, IH = im.size
        ppx, ppy = IW / (e2 - e0), IH / (e3 - e1)
        canvas = Image.new("RGB", (max(1, int(round(totW * ppx))),
                                   max(1, int(round(totH * ppy)))), (0, 0, 0))
        for co in cores_out:
            vx0, vy0, vx1, vy1 = co["src"]
            sx0 = max(0, int(round((vx0 - fp - e0) * ppx)))
            sy0 = max(0, int(round((vy0 - fp - e1) * ppy)))
            sx1 = min(IW, int(round((vx1 + fp - e0) * ppx)))
            sy1 = min(IH, int(round((vy1 + fp - e1) * ppy)))
            if sx1 <= sx0 or sy1 <= sy0:
                continue
            tile = im.crop((sx0, sy0, sx1, sy1))
            canvas.paste(tile, (int(round(co["box"][0] * ppx)), int(round(co["box"][1] * ppy))))
        buf = BytesIO()
        canvas.save(buf, format="JPEG", quality=int(P["jpeg_quality"]))
        D["rawImage"] = "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()

    for co in cores_out:
        co.pop("src", None)

    D["rawExtent"] = [0.0, 0.0, round(totW, 1), round(totH, 1)]
    D["x"] = [round(float(v), 1) for v in nx]
    D["y"] = [round(float(v), 1) for v in ny]
    D["cores"] = cores_out
    D["coreFov"] = {str(f): i for i, co in enumerate(cores_out) for f in co["fovs"]}

    out = out or path
    open(out, "w", encoding="utf-8").write(txt[:s] + json.dumps(D, separators=(",", ":")) + txt[e:])
    ncore = len({co["id"] for co in cores_out})
    print(f"[regrid-core] {path.split('/')[-1]}: {ncore} core_id / {len(ids)} sub-cores -> "
          f"{nrow}x{ncol} grid ({totW:.0f}x{totH:.0f}um)  frame_pad {fp:.0f}um  gutter {gut:.0f}um")
    for co in cores_out:
        tag = co["id"] + (f" ({co['part']}/{co['nparts']})" if co["nparts"] > 1 else "")
        print(f"    {tag:18s} n={co['n']:6d} fovs={co['nfov']:2d} "
              f"box=({co['box'][0]:.0f},{co['box'][1]:.0f})-({co['box'][2]:.0f},{co['box'][3]:.0f})")
    return cores_out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Core-aware regrid using authoritative core_id.")
    ap.add_argument("path")
    ap.add_argument("--out")
    ap.add_argument("--ncol", type=int)
    ap.add_argument("--frame-pad-um", type=float, dest="frame_pad_um")
    ap.add_argument("--gutter-um", type=float, dest="gutter_um")
    a = ap.parse_args()
    regrid(a.path, a.out, ncol=a.ncol,
           frame_pad_um=a.frame_pad_um, gutter_um=a.gutter_um)
