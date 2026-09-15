#!/usr/bin/env python
# =============================================================================
# CGC viewer — split each authoritative core_id into physical SUB-CORES
# =============================================================================
# Description : One TMA core_id can cover two physically separate tissue pieces
#               (verified: <unit> <core_id> = FOVs 52-59 and 60-64, two blobs ~3 mm
#               apart). Framing them as one box would draw a frame across empty
#               slide; framing by the old regrid heuristic alone had no clinical
#               grounding and could group FOVs across different core_ids.
#               This script does both correctly: partition FOVs by AUTHORITATIVE
#               core_id first, then split each core_id into spatially connected
#               components. A component therefore never spans two core_ids and
#               never merges two distant tissue pieces.
# Why pre-degap: degap_viewer.py compacts every empty band to ~one FOV pitch, so
#               a 3 mm inter-piece gap and a single tissue-free FOV hole both end
#               up ~200 um apart and become indistinguishable. Component
#               detection must therefore run on the ORIGINAL build coordinates;
#               the resulting labels ride through degap untouched.
# Input       : viewer HTML straight from build_sample_viewer.py (PRE-degap)
#               fov_meta_qc.tsv  (slide,fov) -> core_id   AUTHORITATIVE
# Output      : same HTML rewritten in place; adds
#               DATA.coreFov  {fov: "<core_id>#1"}   (sub-core key per FOV)
#               DATA.coreSub  {"<core_id>#1": {"core":"<core_id>","part":1,"nparts":2}}
#               The "#k" suffix is an INTERNAL key only. Every sub-core of one
#               punch keeps the parent core_id as its displayed label and shares
#               one clinical record + one header tag + one show/hide checkbox.
# Conda env   : <env>/bin/python  # e.g. a conda env named scvi
# Key deps    : numpy, pandas, scipy
# Validated   : 2026-09-15 (<unit>)
# Parameters  : See DEFAULTS dict below
# =============================================================================

import argparse, json, re

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# -- DEFAULTS (all tunable parameters) ----------------------------------------
DEFAULTS = {
    # Adjacency reach as a multiple of the unit's median FOV pitch. 1.8x is the
    # value regrid_viewer.py used for the dense core interior.
    "near_pitch_mult": 1.8,
    # ...plus a floor in um so a unit with an unusually small pitch still joins
    # touching FOVs.
    "near_floor_um": 200.0,
    # Raster-consecutive FOVs up to this many apart in NUMBER bridge a
    # tissue-free FOV hole (CosMx numbers a core's FOVs consecutively), at the
    # wider reach below. Keeps a hole from faking a split.
    "bridge_max_fov_gap": 3,
    "bridge_pitch_mult": 2.75,
    "bridge_floor_um": 350.0,
    "min_cells_per_part": 30,   # a component smaller than this folds into its nearest sibling
}


def extract_data(txt):
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


def core_map(fov_meta_qc, patient, slide):
    m = pd.read_csv(fov_meta_qc, sep="\t")
    m = m[m["patient_id"].astype(str) == str(patient)]
    m = m[m["tma"].astype(str) == slide.replace("_", "")]
    return {int(f): str(c) for f, c in zip(m["fov"], m["core_id"])}


def components(cen, fids, pitch, P):
    """Union-find over one core_id's FOV centroids -> list of fov-number sets."""
    if len(cen) == 1:
        return [{fids[0]}]
    near = max(pitch * P["near_pitch_mult"], pitch + P["near_floor_um"])
    bridge = max(pitch * P["bridge_pitch_mult"], pitch + P["bridge_floor_um"])
    parent = list(range(len(cen)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        parent[find(a)] = find(b)

    tree = cKDTree(cen)
    for a, b in tree.query_pairs(near):            # touching FOVs
        union(a, b)
    fs = np.asarray(fids)
    for a in range(len(fs)):                       # raster-consecutive hole bridge
        for b in range(a + 1, len(fs)):
            if fs[b] - fs[a] > P["bridge_max_fov_gap"]:
                break
            if float(np.hypot(*(cen[a] - cen[b]))) <= bridge:
                union(a, b)
    comp = {}
    for i in range(len(cen)):
        comp.setdefault(find(i), []).append(i)
    return [{fids[i] for i in members} for members in comp.values()]


def process(path, fov_meta_qc, out=None, **kw):
    P = {**DEFAULTS, **{k: v for k, v in kw.items() if v is not None}}
    txt = open(path, encoding="utf-8").read()
    s, e, D = extract_data(txt)
    if not D.get("global", False):
        raise SystemExit("[err] DATA.global is already false — run this BEFORE degap_viewer.py")

    x = np.array(D["x"], float)
    y = np.array(D["y"], float)
    fov = np.array(D["fov"], int)
    patient = str(D["patient"]["Patient"])
    slide = str(D["sample"]).split()[0]
    f2c = core_map(fov_meta_qc, patient, slide)
    missing = sorted({int(f) for f in fov if int(f) not in f2c})
    if missing:
        raise SystemExit(f"[err] FOVs absent from fov_meta_qc for {patient}/{slide}: {missing[:12]}")

    cen_all = pd.DataFrame({"x": x, "y": y, "fov": fov}).groupby("fov")[["x", "y"]].mean()
    C = cen_all[["x", "y"]].values
    dnn, _ = cKDTree(C).query(C, k=min(2, len(C)))
    pitch = float(np.median(dnn[:, 1])) if len(C) > 1 else 500.0

    core = np.array([f2c[int(f)] for f in fov])
    coreFov, coreSub = {}, {}
    print(f"[subcore] {patient}/{slide}: {len(C)} FOVs, median pitch {pitch:.0f}um")
    for cid in sorted(set(core)):
        fids = sorted({int(f) for f in fov[core == cid]})
        sub = cen_all.loc[fids][["x", "y"]].values
        parts = components(sub, fids, pitch, P)
        # order parts by cell count desc, fold away slivers
        parts = sorted(parts, key=lambda s: -int(np.isin(fov, list(s)).sum()))
        keep = [p for p in parts if int(np.isin(fov, list(p)).sum()) >= P["min_cells_per_part"]]
        if not keep:
            keep = parts[:1]
        for p in parts:
            if p not in keep:
                keep[0] |= p
        # re-order kept parts top-to-bottom then left-to-right for stable naming
        keep = sorted(keep, key=lambda s: (round(float(cen_all.loc[sorted(s)]["y"].mean()) / 500),
                                           float(cen_all.loc[sorted(s)]["x"].mean())))
        for k, pset in enumerate(keep, 1):
            # Sub-cores of one punch SHARE the core_id: the unique key exists only
            # so each physical piece gets its own frame; the displayed label and
            # the clinical record stay the parent core_id.
            key = cid if len(keep) == 1 else f"{cid}#{k}"
            n = int(np.isin(fov, list(pset)).sum())
            coreSub[key] = {"core": cid, "part": k, "nparts": len(keep)}
            for f in sorted(pset):
                coreFov[str(f)] = key
            span = (float(x[np.isin(fov, list(pset))].ptp()), float(y[np.isin(fov, list(pset))].ptp()))
            print(f"    {key:14s} fovs={len(pset):2d} cells={n:6d} span={span[0]:.0f}x{span[1]:.0f}um "
                  f"{'(SPLIT)' if len(keep) > 1 else ''}")
    D["coreFov"] = coreFov
    D["coreSub"] = coreSub
    open(out or path, "w", encoding="utf-8").write(
        txt[:s] + json.dumps(D, separators=(",", ":")) + txt[e:])
    nsplit = sum(1 for v in coreSub.values() if v["nparts"] > 1)
    print(f"[subcore] {len(set(core))} core_id -> {len(coreSub)} sub-cores ({nsplit} from split cores)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Split authoritative core_ids into physical sub-cores.")
    ap.add_argument("path")
    ap.add_argument("--fov-meta-qc", required=True)
    ap.add_argument("--out")
    for k, v in DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), type=type(v), default=None)
    a = ap.parse_args()
    process(a.path, a.fov_meta_qc, a.out,
            **{k: getattr(a, k) for k in DEFAULTS if getattr(a, k) is not None})
