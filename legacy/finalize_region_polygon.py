#!/usr/bin/env python
# Finalize DATA.importantRegion.poly for a per-patient viewer HTML.
# Reads DATA._importantRegion.cellIdx (indexes into DATA.x/y, preserved through degap+regrid),
# computes their convex hull in the FINAL coord space, injects DATA.importantRegion.
# Idempotent: if no _importantRegion or importantRegion already populated, does nothing.
import json, re, sys, argparse
import numpy as np
try:
    from scipy.spatial import ConvexHull
    _HAVE_HULL = True
except ImportError:
    _HAVE_HULL = False

def extract_data(txt):
    m = re.search(r"const DATA\s*=\s*", txt); s = m.end(); d = 0; i = s
    while i < len(txt):
        c = txt[i]
        if c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: i += 1; break
        i += 1
    return s, i, json.loads(txt[s:i])

def convex_hull_2d(pts):
    """Simple 2D convex hull (fallback if scipy missing). pts is Nx2."""
    if len(pts) < 3:
        # not enough points — return a small triangle around centroid
        c = pts.mean(0) if len(pts) else np.array([0.0, 0.0])
        return np.array([c + [ 5, 0], c + [-3, 4], c + [-3, -4]])
    if _HAVE_HULL:
        h = ConvexHull(pts)
        return pts[h.vertices]
    # Andrew's monotone chain
    pts = pts[np.lexsort(pts.T[::-1])]
    def build(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and (h[-1][0] - h[-2][0])*(p[1] - h[-2][1]) - (h[-1][1] - h[-2][1])*(p[0] - h[-2][0]) <= 0:
                h.pop()
            h.append(tuple(p))
        return h
    lo = build(pts); up = build(pts[::-1])
    return np.array(lo[:-1] + up[:-1])

def process(path, out=None, pad_um=8.0):
    txt = open(path).read()
    if 'const DATA=' not in txt and 'const DATA =' not in txt:
        print(f"[skip] {path}: no DATA"); return False
    s, e, D = extract_data(txt)

    # If already finalized, no-op
    if D.get('importantRegion') and isinstance(D['importantRegion'], dict) and D['importantRegion'].get('poly'):
        print(f"[skip] {path.split('/')[-1]}: importantRegion already populated"); return False

    region = D.get('_importantRegion')
    if not region or not region.get('cellIdx'):
        # ensure key exists so template gate is well-defined
        D['importantRegion'] = None
        newjson = json.dumps(D, separators=(',', ':'))
        out = out or path
        open(out, 'w').write(txt[:s] + newjson + txt[e:])
        print(f"[ok] {path.split('/')[-1]}: no _importantRegion; importantRegion=null"); return False

    idx = np.array(region['cellIdx'], dtype=int)
    x = np.array(D['x'], dtype=float)
    y = np.array(D['y'], dtype=float)
    valid = (idx >= 0) & (idx < len(x))
    idx = idx[valid]
    if len(idx) < 3:
        print(f"[skip] {path.split('/')[-1]}: only {len(idx)} inside cells"); return False

    pts = np.column_stack([x[idx], y[idx]])
    # Slight padding for visibility
    c = pts.mean(0)
    pts_padded = pts + (pts - c) / np.maximum(np.linalg.norm(pts - c, axis=1, keepdims=True), 1e-6) * pad_um
    hull = convex_hull_2d(pts_padded)

    region_entry = {
        "pdx": region.get('pdx', 'important region'),
        "cx": float(round(c[0], 1)),
        "cy": float(round(c[1], 1)),
        "n_cells_inside": int(len(idx)),
        "poly": [[float(round(v[0], 1)), float(round(v[1], 1))] for v in hull],
    }
    D['importantRegion'] = region_entry
    D.pop('_importantRegion', None)  # drop intermediate

    newjson = json.dumps(D, separators=(',', ':'))
    out = out or path
    open(out, 'w').write(txt[:s] + newjson + txt[e:])
    print(f"[ok] {path.split('/')[-1]}: importantRegion set  pdx={region_entry['pdx']}  n_inside={len(idx)}  poly_pts={len(hull)}  center=({region_entry['cx']:.0f},{region_entry['cy']:.0f})")
    return True

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('paths', nargs='+', help='viewer HTML file(s)')
    ap.add_argument('--pad-um', type=float, default=8.0)
    a = ap.parse_args()
    for p in a.paths:
        try: process(p, pad_um=a.pad_um)
        except Exception as e: print(f"[err] {p}: {e}")
