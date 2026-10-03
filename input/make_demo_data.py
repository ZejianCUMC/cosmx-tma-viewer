#!/usr/bin/env python
# =============================================================================
# cosmx-tma-viewer : synthetic demo cohort
# =============================================================================
# Description : Writes a small, fully synthetic input set (no real data, no real metadata) in the formats read by
#               build_viewer.py: cells, cores, units, FOV positions, a CosMx-style count export, overview and
#               full-resolution FOV images, H&E-like core images, regions and index flags.
# Output      : <out_dir>/{cells.tsv.gz, cores.tsv, units.tsv, fov_positions.tsv, regions.tsv, flags.tsv,
#               counts/{counts.mtx,genes.txt,cells.txt}, fov_images/, fov_images_hires/, he/}
# Key deps    : numpy, pandas, Pillow
# =============================================================================
import argparse, math
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image, ImageDraw, ImageFilter

DEFAULTS = {"seed": 7, "pixel_size_um": 0.12028, "fov_size_px": 4256, "cells_per_core": 2500, "n_genes": 60,
            "overview_tile_px": 256, "hires_tile_px": 1024}
TYPES = {  # L1 -> {L2 -> [L3]}; generic tissue classes
    "Epithelial": {"Luminal-like": ["Epi_A", "Epi_B"], "Basal-like": ["Epi_C"]},
    "Stromal": {"Fibroblast": ["Fib_A", "Fib_B"], "Pericyte": ["Pericyte"]},
    "Immune": {"Lymphoid": ["T_like", "B_like"], "Myeloid": ["Mac_like"]},
    "LowQuality": {"LowQuality": ["LowQuality"]}}
UNITS = [  # unit, slide, cores [(core, grid origin col,row, n_fov_x, n_fov_y, group)]
    ("DEMO_A", "S1", [("S1_C1", 0, 0, 3, 3, "group_1"), ("S1_C2", 5, 0, 3, 2, "group_2")]),
    ("DEMO_B", "S2", [("S2_C1", 0, 0, 3, 3, "group_1")])]

def run_pipeline(out_dir, **kw):
    P = {**DEFAULTS, **kw}; rng = np.random.default_rng(P["seed"]); out = Path(out_dir)
    for d in ("counts", "fov_images", "fov_images_hires", "he"): (out / d).mkdir(parents=True, exist_ok=True)
    PX, FW = P["pixel_size_um"], P["fov_size_px"]; fov_um = PX * FW
    l3 = [(a, b, c) for a in TYPES for b in TYPES[a] for c in TYPES[a][b]]
    genes = [f"GENE{i:02d}" for i in range(1, P["n_genes"] + 1)]
    markers = {c: rng.choice(len(genes), 4, replace=False) for _, _, c in l3}
    cells, cores, fovs, regions, flags, triplets = [], [], [], [], [], []
    for unit, slide, clist in UNITS:
        fov_id = 1
        for core, gx0, gy0, nx, ny, grp in clist:
            fl = []
            for j in range(ny):
                for i in range(nx):
                    fovs.append(dict(slide=slide, fov=fov_id, x_global_px=(gx0 + i) * FW, y_global_px=-(gy0 + j) * FW)); fl.append(fov_id); fov_id += 1
            cx, cy = (gx0 + nx / 2) * fov_um, (gy0 + ny / 2) * fov_um; rad = 0.45 * min(nx, ny) * fov_um
            pts = []
            while len(pts) < P["cells_per_core"]:
                x, y = rng.uniform(cx - rad, cx + rad), rng.uniform(cy - rad, cy + rad)
                if (x - cx) ** 2 + (y - cy) ** 2 < rad ** 2: pts.append((x, y))
            pts = np.array(pts); ang = np.arctan2(pts[:, 1] - cy, pts[:, 0] - cx)
            for k, (x, y) in enumerate(pts):
                i, j = int((x / fov_um) - gx0), int((y / fov_um) - gy0); i, j = min(max(i, 0), nx - 1), min(max(j, 0), ny - 1)
                f = fl[j * nx + i]; fr = [r for r in fovs if r["slide"] == slide and r["fov"] == f][0]
                xl = x / PX - fr["x_global_px"]; yl = y / PX + fr["y_global_px"]      # inverse of X=PX*(xg+xl), Y=PX*(-yg+yl)
                sector = int((ang[k] + math.pi) / (2 * math.pi) * 6) % 6
                a, b, c = l3[min(len(l3) - 1, sector + (rng.random() < 0.3) * 3)] if rng.random() > 0.08 else ("LowQuality", "LowQuality", "LowQuality")
                cid = f"{slide}::c_{f}_{k}"; n = int(rng.integers(40, 900)) if a != "LowQuality" else int(rng.integers(10, 60))
                cells.append(dict(cell=cid, slide=slide, fov=f, core=core, x_local_px=round(xl, 1), y_local_px=round(yl, 1), L1=a, L2=b, L3=c,
                                  mean_dapi=round(float(rng.gamma(4, 200)), 1), mean_panck=round(float(rng.gamma(2, 300 if a == "Epithelial" else 60)), 1),
                                  mean_cd45=round(float(rng.gamma(2, 300 if a == "Immune" else 50)), 1), mean_membrane=round(float(rng.gamma(3, 150)), 1)))
                ci = len(cells)
                lam = np.full(len(genes), n / len(genes) * 0.3); lam[markers[c]] += n / 12
                for g, v in enumerate(rng.poisson(lam)):
                    if v: triplets.append((g + 1, ci, int(v)))
            cores.append(dict(unit=unit, core=core, slide=slide, group=grp, annotation=rng.choice(["Type I", "Type II"]), fovs=",".join(map(str, fl))))
            # synthetic morphology: dark field + coloured blobs at cell positions, per FOV
            for f in fl:
                fr = [r for r in fovs if r["slide"] == slide and r["fov"] == f][0]
                for d, size in (("fov_images", P["overview_tile_px"]), ("fov_images_hires", P["hires_tile_px"])):
                    im = Image.new("RGB", (size, size), (6, 10, 8)); dr = ImageDraw.Draw(im); s = size / FW
                    for r in (cc for cc in cells if cc["core"] == core and cc["fov"] == f):
                        u, v = r["x_local_px"] * s, r["y_local_px"] * s; rr = max(1, 55 * s)
                        col = (40, 200, 60) if r["L1"] == "Epithelial" else ((60, 90, 230) if r["L1"] == "Immune" else (150, 80, 60))
                        dr.ellipse([u - rr, v - rr, u + rr, v + rr], fill=col)
                    im.filter(ImageFilter.GaussianBlur(max(0.5, 8 * s))).save(out / d / f"{slide}_F{f:05d}.jpg", quality=88)
            he = Image.new("RGB", (600, 600), (245, 240, 242)); dr = ImageDraw.Draw(he)
            for _ in range(900):
                u, v = rng.uniform(80, 520, 2)
                if (u - 300) ** 2 + (v - 300) ** 2 < 220 ** 2: dr.ellipse([u - 4, v - 4, u + 4, v + 4], fill=(120, 60, 150))
            he.filter(ImageFilter.GaussianBlur(1.2)).save(out / "he" / f"{core}.jpg", quality=90)
            if core == "S1_C1":
                regions.append(dict(unit=unit, core=core, id="R1", tier=1, cx_um=cx + 0.3 * rad, cy_um=cy, r_um=0.18 * rad, summary="synthetic aggregate"))
                flags.append(dict(unit=unit, core=core, label="Example flag", detail="synthetic flag for the demo", tooltip="Shown on the landing page"))
    pd.DataFrame(cells).to_csv(out / "cells.tsv.gz", sep="\t", index=False)
    pd.DataFrame(cores).to_csv(out / "cores.tsv", sep="\t", index=False)
    pd.DataFrame([dict(unit=u, cohort_arm=a, age_band=b) for (u, _, _), a, b in zip(UNITS, ["X", "Y"], ["60-69", "70-79"])]).to_csv(out / "units.tsv", sep="\t", index=False)
    pd.DataFrame(fovs).to_csv(out / "fov_positions.tsv", sep="\t", index=False)
    pd.DataFrame(regions).to_csv(out / "regions.tsv", sep="\t", index=False); pd.DataFrame(flags).to_csv(out / "flags.tsv", sep="\t", index=False)
    with open(out / "counts" / "counts.mtx", "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n"); fh.write(f"{len(genes)} {len(cells)} {len(triplets)}\n")
        fh.writelines(f"{g} {c} {v}\n" for g, c, v in triplets)
    (out / "counts" / "genes.txt").write_text("\n".join(genes) + "\n"); (out / "counts" / "cells.txt").write_text("\n".join(c["cell"] for c in cells) + "\n")
    print(f"demo data: {len(cells):,} cells, {len(cores)} cores, {len(UNITS)} units, {len(fovs)} FOVs -> {out}")
    return {"output_path": str(out), "params_used": P}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="write a synthetic demo input set"); ap.add_argument("out_dir")
    run_pipeline(ap.parse_args().out_dir)
