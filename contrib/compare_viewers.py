#!/usr/bin/env python
# =============================================================================
# cosmx-tma-viewer : compare two builds of the same viewer (regression check)
# =============================================================================
# Description : For every core present in both viewers: same number of cells, cells matched by nearest neighbour after
#               fitting a per-axis scale + shift (so a re-layout or a corrected placement is not reported as a change),
#               identical L1/L2/L3 labels and FOV ids, identical decoded gene values for a set of genes, same FOV lists,
#               and the number of cells inside each region polygon. Reports the fitted scale per core, so a geometry
#               change (e.g. a core that was stretched in one build) is visible.
# Usage       : python contrib/compare_viewers.py OLD.html NEW.html [--genes G1 G2 ...] [--tol 0.05]
# Key deps    : numpy, scipy, matplotlib (Path only)
# =============================================================================
import argparse, base64, json, re, sys
import numpy as np
from scipy.spatial import cKDTree
from matplotlib.path import Path as MP

def load(f):
    t = open(f, encoding="utf-8").read(); m = re.search(r"const DATA\s*=\s*", t); return json.JSONDecoder().raw_decode(t, m.end())[0]

def genes_of(V, names):
    g = V.get("geneSparse")
    if not g: return np.zeros((V["n"], 0))
    ip = np.frombuffer(base64.b64decode(g["indptr"]), "<i4"); ix = np.frombuffer(base64.b64decode(g["indices"]), "<u2")
    va = np.frombuffer(base64.b64decode(g["data"]), "u1"); rows = np.repeat(np.arange(g["n"]), np.diff(ip)); gm = np.array(g["geneMax"])
    out = np.zeros((g["n"], len(names)))
    for j, n in enumerate(names):
        if n in g["names"]:
            x = g["names"].index(n); m = ix == x; out[rows[m], j] = va[m] / 255 * gm[x]
    return out

def run_pipeline(old, new, genes=None, tol=0.05):
    A, B = load(new), load(old); ok = True; rows = []
    for V in (A, B): V["_core"] = np.array([V["cores"][V["coreFov"][str(f)]]["id"] for f in V["fov"]])
    if genes is None and A.get("geneSparse"): genes = A["geneSparse"]["names"][:6]
    ga, gb = genes_of(A, genes or []), genes_of(B, genes or [])
    for c in A["cores"]:
        ia = np.where(A["_core"] == c["id"])[0]; ib = np.where(B["_core"] == c["id"])[0]
        if not len(ib): rows.append(f"{c['id']}: only in NEW"); ok = False; continue
        if len(ia) != len(ib): rows.append(f"{c['id']}: n {len(ib)} -> {len(ia)}"); ok = False; continue
        xa, ya = np.array(A["x"])[ia], np.array(A["y"])[ia]; xb, yb = np.array(B["x"])[ib], np.array(B["y"])[ib]
        sx, sy = xb.std() / xa.std(), yb.std() / ya.std()
        d, j = cKDTree(np.c_[xb - xb.mean(), yb - yb.mean()]).query(np.c_[(xa - xa.mean()) * sx, (ya - ya.mean()) * sy])
        for _ in range(2):
            kx = np.polyfit(xa, xb[j], 1); ky = np.polyfit(ya, yb[j], 1)
            d, j = cKDTree(np.c_[xb, yb]).query(np.c_[np.polyval(kx, xa), np.polyval(ky, ya)])
        same = all((np.array(A[k])[ia] == np.array(B[k])[ib][j]).all() for k in ("type", "subL2", "subtype", "fov"))
        gdiff = float(np.abs(ga[ia] - gb[ib][j]).max()) if ga.shape[1] else 0.0
        fovs = sorted(c["fovs"]) == sorted(next(x for x in B["cores"] if x["id"] == c["id"])["fovs"])
        good = d.max() < tol * max(1, kx[0], ky[0]) and len(set(j)) == len(j) and same and gdiff < 0.05 and fovs; ok &= good
        rows.append(f"{c['id']}: n={len(ia)} scale(old/new x,y)={kx[0]:.4f}/{ky[0]:.4f} resid={d.max():.3f} labels={'same' if same else 'DIFF'} "
                    f"gene|diff|max={gdiff:.3f} fovs={'same' if fovs else 'DIFF'}")
    for R in A.get("tlsRegions", []):
        rb = next((x for x in B.get("tlsRegions", []) if x["id"] == R["id"]), None)
        if rb is None: rows.append(f"region {R['id']}: only in NEW"); continue
        n = lambda V, Q: int(MP(np.array(Q["poly"])).contains_points(np.c_[np.array(V["x"]), np.array(V["y"])][V["_core"] == Q["core"]]).sum())
        a, b = n(A, R), n(B, rb); ok &= abs(a - b) <= 1; rows.append(f"region {R['id']}: cells inside {b} -> {a}")
    print(("SAME " if ok else "DIFF ") + f"{old} vs {new}"); [print("   " + r) for r in rows]
    return {"ok": ok, "summary": rows}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="compare two builds of the same viewer"); ap.add_argument("old"); ap.add_argument("new")
    ap.add_argument("--genes", nargs="*"); ap.add_argument("--tol", type=float, default=0.05)
    a = ap.parse_args(); sys.exit(0 if run_pipeline(a.old, a.new, a.genes, a.tol)["ok"] else 1)
