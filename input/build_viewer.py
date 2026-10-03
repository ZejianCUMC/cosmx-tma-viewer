#!/usr/bin/env python
# =============================================================================
# cosmx-tma-viewer : build self-contained HTML viewers (one per unit)
# =============================================================================
# Description : Builds one offline HTML viewer per unit (e.g. patient x slide) from plain tables:
#               cells (FOV-local positions + 3-level labels + optional channels), cores, units, FOV positions,
#               optional CosMx count export, per-FOV morphology images (overview and optional full resolution),
#               optional per-core H&E images and optional circular regions. Every core is placed in its own
#               layout slot by ONE translation of its native slide coordinates, so cells, morphology tiles and
#               regions share exactly the same transform. Study-specific wording, colours and fields come from
#               the cohort config (see cohort_config.example.json); nothing study-specific is in the code.
# Input       : see docs/inputs.md
# Output      : <out_dir>/<file_prefix><unit>.html (+ build_summary.tsv)
# Key deps    : numpy, pandas, scipy, Pillow
# Parameters  : DEFAULTS below; every key can be overridden in the config JSON
# =============================================================================
import argparse, base64, io, json, math, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

# -- DEFAULTS (all tunable; config keys of the same name override) ------------------------------------------
DEFAULTS = {
    "title": "CosMx TMA viewer",
    "file_prefix": "viewer_",
    "pixel_size_um": 0.12028,          # um per global pixel (CosMx)
    "fov_size_px": 4256,               # FOV side in global pixels (CosMx)
    "layout": {"ncol": 3,              # core slots per row
               "slot_um": 1760,        # minimum slot side; grown to fit the largest core + margins
               "cell_margin_um": 20,   # frame margin around the outermost cells
               "align_rows": []},      # lists of core ids per unit that share one frame top/bottom
    "morphology": {"overview_px_per_um": 0.5, "overview_jpeg_quality": 88,
                   "hires_px_per_um": 4.0, "hires_jpeg_quality": 85},
    "levels": ["Cell type · L1", "Cell type · L2", "Cell type · L3"],
    "hidden_default": [],              # category values hidden in every view by default (e.g. a QC-fail class)
    "palette": {},                     # category value -> colour; others get a 20-colour cycle
    "channels": {},                    # display name -> column in cells table (e.g. "DAPI": "mean_dapi")
    "core_groups": {},                 # group key -> {label, short, color}; frame colour per core group
    "core_group_title": "",            # header chip title for the group summary ("" = no chip)
    "core_short_fields": [],           # core fields joined into the short text on tags / hover
    "core_tooltip_fields": [],         # [{label, field}] lines in the core-label tooltip
    "core_summary_fields": [],         # [{label, field}] header chips aggregated over cores
    "patient_fields": [],              # [{label, field, format}] header chips from the units table
    "patient_warn": {},                # chip label -> regex; matching chips are highlighted
    "patient_summary_after": "",       # chip label after which the core summary chips are inserted
    "storage_prefix": "tmaviewer",     # browser localStorage namespace
    "he": {"max_px": 850, "default_pose": {"angle": 0, "flip": 1, "fy": 1}, "pose_rules": []},
    "region": {"label": "Region", "toggle": "Important region", "note": "", "tiers": {},
               "info_fields": [], "n_poly": 96},
    "left_default": "__rawimg__", "right_default": "__he__",
}

def deep_update(a, b):
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(a.get(k), dict): deep_update(a[k], v)
        elif not k.startswith("_"): a[k] = v
    return a

def uri(img, fmt="JPEG", q=88):
    buf = io.BytesIO(); img.save(buf, format=fmt, **({"quality": q} if fmt == "JPEG" else {}))
    return f"data:image/{fmt.lower()};base64," + base64.b64encode(buf.getvalue()).decode()

def find_image(d, stem):
    for ext in (".jpg", ".jpeg", ".png", ".tif", ".tiff"):
        p = Path(d) / f"{stem}{ext}"
        if p.exists(): return p
    return None

# -- HELPERS --------------------------------------------------------------------------------------------------
def load_counts(cdir, cell_ids):
    """CosMx-style export: counts.mtx (genes x cells, MatrixMarket), genes.txt, cells.txt -> CSR cells x genes."""
    from scipy.io import mmread
    from scipy.sparse import csr_matrix
    genes = [l.strip() for l in open(Path(cdir) / "genes.txt")]
    cells = [l.strip() for l in open(Path(cdir) / "cells.txt")]
    pos = pd.Series(np.arange(len(cells)), index=cells)
    miss = [c for c in cell_ids if c not in pos.index]
    if miss: sys.exit(f"{len(miss)} cells missing from the count export, e.g. {miss[:3]}")
    X = csr_matrix(mmread(Path(cdir) / "counts.mtx").T.tocsr()).astype(np.float32)
    return X, genes, pos

def gene_sparse(X, rows, genes):
    """normalise to 1e4 + log1p, then quantise per gene to uint8 codes (1..255 of the gene's max in this unit)."""
    sub = X[rows].tocsr().copy(); tot = np.asarray(sub.sum(1)).ravel(); tot[tot == 0] = 1
    sub = sub.multiply(1e4 / tot[:, None]).tocsr(); sub.data = np.log1p(sub.data)
    gmax = np.asarray(sub.max(axis=0).todense()).ravel().astype(np.float32); gmax[gmax == 0] = 1.0
    codes = np.clip(np.rint(sub.data * (255.0 / gmax[sub.indices])), 1, 255).astype(np.uint8)
    b = lambda a: base64.b64encode(a.tobytes()).decode("ascii")
    assert len(genes) < 65536
    gm = [round(float(x), 4) for x in gmax]
    return {"names": genes, "geneMax": gm, "n": int(sub.shape[0]), "indptr": b(sub.indptr.astype(np.int32)),
            "indices": b(sub.indices.astype(np.uint16)), "data": b(codes)}, gm

def paint(canvas, origin, R, tiles, clip):
    """paste tiles [(img, x0,y0,x1,y1 in viewer um)] into canvas whose pixel (0,0) is viewer um `origin`, clipped to `clip` box."""
    W, H = canvas.size; ox, oy = origin
    cx0, cy0, cx1, cy1 = [int(round((v - o) * R)) for v, o in zip(clip, (ox, oy, ox, oy))]
    for img, x0, y0, x1, y1 in tiles:
        X0, Y0, X1, Y1 = [int(round((v - o) * R)) for v, o in zip((x0, y0, x1, y1), (ox, oy, ox, oy))]
        bx0, by0, bx1, by1 = max(X0, cx0, 0), max(Y0, cy0, 0), min(X1, cx1, W), min(Y1, cy1, H)
        if bx1 <= bx0 or by1 <= by0 or X1 <= X0 or Y1 <= Y0: continue
        t = img.resize((X1 - X0, Y1 - Y0), Image.LANCZOS)
        canvas.paste(t.crop((bx0 - X0, by0 - Y0, bx1 - X0, by1 - Y0)), (bx0, by0))

def fmt_value(v, f):
    if v is None or (isinstance(v, float) and math.isnan(v)) or str(v) in ("", "nan", "NA", "None"): return None
    if f:
        try: return f.format(float(v))
        except (ValueError, TypeError): return str(v)
    return str(v)

# -- MAIN PIPELINE --------------------------------------------------------------------------------------------
def run_pipeline(cells_tsv, cores_tsv, units_tsv, fov_tsv, out_dir, config=None, counts_dir=None, fov_images=None,
                 fov_images_hires=None, he_images=None, regions_tsv=None, units=None, template=None):
    P = deep_update(json.loads(json.dumps(DEFAULTS)), config or {})
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    tpl = Path(template or Path(__file__).with_name("viewer_template.html")).read_text(encoding="utf-8")
    C = pd.read_csv(cells_tsv, sep="\t", low_memory=False, dtype={"cell": str, "core": str, "slide": str})
    K = pd.read_csv(cores_tsv, sep="\t", dtype=str).fillna("")
    U = pd.read_csv(units_tsv, sep="\t", dtype=str).fillna("").set_index("unit")
    F = pd.read_csv(fov_tsv, sep="\t", dtype={"slide": str}).set_index(["slide", "fov"])
    G = pd.read_csv(regions_tsv, sep="\t", dtype={"core": str}) if regions_tsv else pd.DataFrame(columns=["unit", "core"])
    need = {"cell", "slide", "fov", "core", "x_local_px", "y_local_px", "L1", "L2", "L3"}
    assert need <= set(C.columns), f"cells table needs {sorted(need)}"
    assert set(C.core) <= set(K.core), "every cell core must be listed in the cores table"
    PX, FW, M = P["pixel_size_um"], P["fov_size_px"], P["layout"]["cell_margin_um"]
    xg = F.loc[list(zip(C.slide, C.fov)), "x_global_px"].values; yg = F.loc[list(zip(C.slide, C.fov)), "y_global_px"].values
    C["gx"] = PX * (xg + C.x_local_px.values); C["gy"] = PX * (-yg + C.y_local_px.values)     # native um, y down
    X = genes = cpos = None
    if counts_dir: X, genes, cpos = load_counts(counts_dir, C.cell.tolist())
    CFG = {"hiddenDefault": P["hidden_default"], "levelLabels": P["levels"], "storagePrefix": P["storage_prefix"],
           "palette": P["palette"], "coreGroups": P["core_groups"], "coreGroupTitle": P["core_group_title"],
           "coreSummaryFields": P["core_summary_fields"], "patientWarn": P["patient_warn"],
           "patientSummaryAfter": P["patient_summary_after"], "heDefaultPose": P["he"]["default_pose"],
           "hePoseRules": P["he"]["pose_rules"], "leftDefault": P["left_default"], "rightDefault": P["right_default"],
           "region": {"label": P["region"]["label"], "note": P["region"]["note"], "tiers": P["region"]["tiers"]}}
    summary = []
    for unit in (units or list(U.index)):
        Ku = K[K.unit == unit].reset_index(drop=True); Cu = C[C.core.isin(Ku.core)]
        if Cu.empty: print(f"[skip] {unit}: no cells"); continue
        slides = sorted(set(Ku.slide)); assert len(slides) == 1, f"{unit}: one slide per unit (fov ids are per slide)"
        slide = slides[0]
        # --- layout: content box per core, centred in a grid slot -------------------------------------------
        boxes = {}
        for k in Ku.core:
            c = Cu[Cu.core == k]; boxes[k] = [c.gx.min() - M, c.gy.min() - M, c.gx.max() + M, c.gy.max() + M]
        side = max([P["layout"]["slot_um"]] + [max(b[2] - b[0], b[3] - b[1]) + 2 * M for b in boxes.values()])
        nc = min(P["layout"]["ncol"], len(Ku)); nr = math.ceil(len(Ku) / nc)
        cores, shift, parts = [], {}, []
        for i, k in enumerate(Ku.core):
            s = [(i % nc) * side, (i // nc) * side, (i % nc + 1) * side, (i // nc + 1) * side]; b = boxes[k]
            dx = (s[0] + s[2]) / 2 - (b[0] + b[2]) / 2; dy = (s[1] + s[3]) / 2 - (b[1] + b[3]) / 2; shift[k] = (dx, dy)
            row = Ku[Ku.core == k].iloc[0]; fields = {f: row[f] for f in Ku.columns if f not in ("unit", "core", "slide", "group")}
            cu = Cu[Cu.core == k].sort_values(["fov"], kind="stable")
            parts.append(cu.assign(vx=cu.gx + dx, vy=cu.gy + dy, ci=i))
            fovs = sorted({int(v) for v in str(row.get("fovs", "")).replace(";", ",").split(",") if v.strip()} | set(cu.fov.astype(int)))
            short = "; ".join(str(row[f]) for f in P["core_short_fields"] if row.get(f, ""))
            grp = P["core_groups"].get(row.get("group", ""), {})
            lines = [f'{t["label"]}: {row[t["field"]]}' for t in P["core_tooltip_fields"] if row.get(t["field"], "")]
            if grp: lines.append(grp.get("label", ""))
            cores.append({"key": k, "id": k, "part": 1, "nparts": 1, "n": int(len(cu)), "nfov": len(fovs), "fovs": fovs,
                          "pending_fovs": [], "box": [round(v, 2) for v in (b[0] + dx, b[1] + dy, b[2] + dx, b[3] + dy)],
                          "slot": s, "group": row.get("group", ""), "short": short, "lines": lines, "fields": fields})
        for row_ids in [r for r in P["layout"]["align_rows"] if set(r) <= set(Ku.core)]:
            sel = [c for c in cores if c["id"] in row_ids]; top = min(c["box"][1] for c in sel); bot = max(c["box"][3] for c in sel)
            for c in sel: c["box"][1], c["box"][3] = top, bot
        A = pd.concat(parts); W_, H_ = nc * side, nr * side
        # --- morphology: every FOV of every core, same transform as its cells, clipped to the core frame ----
        def fov_rect(fv, k):
            x0 = PX * F.loc[(slide, fv), "x_global_px"]; y0 = -PX * F.loc[(slide, fv), "y_global_px"]
            dx, dy = shift[k]; return x0 + dx, y0 + dy, x0 + PX * FW + dx, y0 + PX * FW + dy
        R = P["morphology"]["overview_px_per_um"]; raw = Image.new("RGB", (math.ceil(W_ * R), math.ceil(H_ * R)), (0, 0, 0))
        n_img = 0
        if fov_images:
            for c in cores:
                tiles = []
                for fv in c["fovs"]:
                    p = find_image(fov_images, f"{slide}_F{fv:05d}")
                    if p is None or (slide, fv) not in F.index: continue
                    tiles.append((Image.open(p).convert("RGB"),) + fov_rect(fv, c["id"])); n_img += 1
                paint(raw, (0, 0), R, tiles, c["box"])
        rawHi = []
        if fov_images_hires:
            Rh = P["morphology"]["hires_px_per_um"]
            for c in cores:
                b = c["box"]; img = Image.new("RGB", (math.ceil((b[2] - b[0]) * Rh), math.ceil((b[3] - b[1]) * Rh)), (0, 0, 0)); tiles = []
                for fv in c["fovs"]:
                    p = find_image(fov_images_hires, f"{slide}_F{fv:05d}")
                    if p is not None and (slide, fv) in F.index: tiles.append((Image.open(p).convert("RGB"),) + fov_rect(fv, c["id"]))
                if not tiles: continue
                paint(img, (b[0], b[1]), Rh, tiles, b)
                rawHi.append({"core": c["id"], "ext": [b[0], b[1], round(b[0] + img.size[0] / Rh, 3), round(b[1] + img.size[1] / Rh, 3)],
                              "src": uri(img, "JPEG", P["morphology"]["hires_jpeg_quality"])})
        he = {}
        if he_images:
            for c in cores:
                p = find_image(he_images, c["id"])
                if p is None: continue
                im = Image.open(p).convert("RGB"); im.thumbnail((P["he"]["max_px"],) * 2); he[c["id"]] = uri(im, "JPEG", 90)
        # --- regions (circles in native um -> polygons through the core transform) ---------------------------
        regs = []
        for _, r in G[G.unit == unit].iterrows():
            k = r.core; dx, dy = shift[k]; n = int(P["region"]["n_poly"]); rr = float(r.r_um) / math.cos(math.pi / n)
            th = np.linspace(0, 2 * math.pi, n, endpoint=False)
            poly = [[round(float(r.cx_um) + rr * math.cos(a) + dx, 2), round(float(r.cy_um) + rr * math.sin(a) + dy, 2)] for a in th]
            lines = [f'{t["label"]}: {r[t["field"]]}' for t in P["region"]["info_fields"] if str(r.get(t["field"], "")) not in ("", "nan")]
            regs.append({"id": str(r.get("id", "")), "tier": int(r.get("tier", 1)), "core": k, "poly": poly,
                         "lx": round(float(r.cx_um) + dx, 2), "ly": round(float(r.cy_um) - rr + dy, 2), "lines": lines})
        # --- expression ----------------------------------------------------------------------------------------
        gs, gmax = (None, None)
        if X is not None: gs, gmax = gene_sparse(X, cpos[A.cell].values, genes)
        ch = {name: np.round(A[col].astype(float).values, 1).tolist() for name, col in P["channels"].items() if col in A.columns}
        urow = U.loc[unit] if unit in U.index else pd.Series(dtype=str)
        pat = {}
        for f in P["patient_fields"]:
            v = fmt_value(urow.get(f["field"], None), f.get("format"))
            if v is not None: pat[f["label"]] = v
        N = len(A)
        DATA = {"sample": f"{unit} · {slide}", "global": False, "patient": pat, "cores": cores,
                "coreFov": {str(fv): i for i, c in enumerate(cores) for fv in c["fovs"]}, "coreSub": {}, "heImages": he,
                "geneSource": [0] * N, "geneMaxBySource": [gmax] if gmax else [],
                "x": np.round(A.vx.values, 2).tolist(), "y": np.round(A.vy.values, 2).tolist(), "fov": A.fov.astype(int).tolist(),
                "type": A.L1.astype(str).tolist(), "subtype": A.L3.astype(str).tolist(), "subL2": A.L2.astype(str).tolist(),
                "channels": ch, "n": N, "fovs": int(A.fov.nunique()), "rawExtent": [0, 0, W_, H_],
                "rawImage": uri(raw, "JPEG", P["morphology"]["overview_jpeg_quality"]), "tlsRegions": regs, "rawHi": rawHi}
        if gs: DATA["geneSparse"] = gs
        html = (tpl.replace("__TITLE__", P["title"]).replace("__REGION_TOGGLE__", P["region"]["toggle"])
                   .replace("/*__CONFIG__*/", json.dumps(CFG, ensure_ascii=False))
                   .replace("/*__DATA__*/", json.dumps(DATA, separators=(",", ":"), ensure_ascii=False)))
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", unit); fn = out / f'{P["file_prefix"]}{safe}.html'; fn.write_text(html, encoding="utf-8")
        summary.append({"unit": unit, "slide": slide, "cores": len(cores), "cells": N, "fov_images": n_img, "hires_cores": len(rawHi),
                        "he_cores": len(he), "regions": len(regs), "MB": round(len(html) / 1e6, 1), "file": fn.name})
        print(f"[ok] {fn.name}: {N:,} cells, {len(cores)} cores, {n_img} FOV images, {len(rawHi)} full-res cores, {len(regs)} regions, {len(html)/1e6:.1f} MB", flush=True)
    pd.DataFrame(summary).to_csv(out / "build_summary.tsv", sep="\t", index=False)
    return {"output_path": str(out), "params_used": P, "summary": summary}

# -- CLI ENTRY POINT --------------------------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__ or "Build cosmx-tma-viewer HTML files")
    ap.add_argument("--cells", required=True); ap.add_argument("--cores", required=True); ap.add_argument("--units", required=True)
    ap.add_argument("--fov-positions", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--config", help="cohort config JSON (see cohort_config.example.json)")
    ap.add_argument("--counts-dir", help="CosMx export dir with counts.mtx (genes x cells), genes.txt, cells.txt")
    ap.add_argument("--fov-images", help="dir of <slide>_F<fov:05d>.jpg overview morphology images")
    ap.add_argument("--fov-images-hires", help="dir of full-resolution FOV images for the zoom-in layer")
    ap.add_argument("--he-images", help="dir of <core>.jpg H&E images")
    ap.add_argument("--regions", help="regions TSV (unit, core, id, tier, cx_um, cy_um, r_um, info fields)")
    ap.add_argument("--only", nargs="*", help="build only these units"); ap.add_argument("--template")
    a = ap.parse_args()
    cfg = json.load(open(a.config)) if a.config else {}
    run_pipeline(a.cells, a.cores, a.units, a.fov_positions, a.out, cfg, a.counts_dir, a.fov_images, a.fov_images_hires,
                 a.he_images, a.regions, a.only, a.template)
