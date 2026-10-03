#!/usr/bin/env python
# =============================================================================
# CosMx morphology mosaic stitcher — ONE unit  (HPC compute node)
# =============================================================================
# Description : Stitch a (patient, slide) unit's per-FOV CellComposite JPGs into
#               a global-micron mosaic, then base64-embed it as {rawImage,
#               rawExtent} for the spatial viewer. Run as one array task per unit
#               by stitch_morphology.sbatch (NEVER on the login node).
# Input       : viewer_units.tsv  (unit_id, slide, fovs)  -- shipped from local
#               flatFiles/<run_prefix><slide>/*_fov_positions_file.csv.gz
#               DecodedFiles/<run_prefix><slide>/*/CellStatsDir/CellComposite/
#                            CellComposite_F<fov:05d>.jpg
# Output      : rawimg_build/rawimg_<unit_id>.json = {rawImage, rawExtent, meta}
# Conda env   : <hpc_project_root>/mamba/envs/GenomicTools/bin/python
# Key deps    : Pillow, numpy, pandas
# ----------------------------------------------------------------------------
# COORDINATE MODEL (verified 2026-08-26):
#   world micron = global_px * UMPX (0.12028)
#   x_global_px = FOV left edge ; y_global_px = FOV top edge
#   per-cell:  gx = x0 + x_FOV_px ,  gy = y0 - y_FOV_px   (CosMx local_y grows
#              as global_y shrinks) -> each FOV is pasted VERTICALLY FLIPPED.
#   rawExtent = [xmin*UMPX, (ymin0-FOVW)*UMPX, (xmax0+FOVW)*UMPX, ymax0*UMPX]
#               (original micron; UNCHANGED by the flip below).
#
# Y-FLIP FIX (critical, 2026-08-26): the viewer JS flips Y for the CELL scatter
#   when DATA.global is true, but NOT for the drawImage morphology layer. So the
#   assembled mosaic must be flipped top-to-bottom (np.flipud) ONCE MORE before
#   embedding, or the morphology renders vertically mirrored vs the cells.
#   Verified by per-cell PanCK vs local-green correlation: as-is = -0.20 (wrong),
#   flipud = +0.23 (correct), 93% of cells land on tissue.
# =============================================================================

import os, sys, glob, json, base64, argparse
import numpy as np, pandas as pd
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

# AtoMx run-folder prefix is acquisition-specific; see cohort_config.json.
_RUN_PREFIX = ((_cohort_cfg().get("naming") or {}).get("run_prefix") or "")

DEFAULTS = {
    "dest":  "<hpc_project_dir>",
    "units": None,          # viewer_units.tsv (defaults to alongside this script)
    "out":   None,          # output dir (defaults to <dest>/rawimg_build)
    "umpx":        0.12028, # microns per global pixel
    "fovw":        4256,    # CellComposite side length (px)
    "target_long": 2500,    # mosaic long-side (px) after downscale
    "quality":     82,      # output JPEG quality
}


def slide_run(slide):
    """obs slide label 'TMA_1' -> AtoMx run folder '<run_prefix>TMA_1'.
    The run prefix is acquisition-specific: set naming.run_prefix in cohort_config.json."""
    return f"{_RUN_PREFIX}{slide}"


def find_composite_dir(dest, slide):
    run = slide_run(slide)
    hits = glob.glob(f"{dest}/DecodedFiles/{run}/*/CellStatsDir/CellComposite")
    if not hits:
        raise FileNotFoundError(f"no CellComposite dir for {slide} under {dest}/DecodedFiles/{run}")
    return sorted(hits)[0]


def load_positions(dest, slide):
    run = slide_run(slide)
    p = f"{dest}/flatFiles/{run}/{run}_fov_positions_file.csv.gz"
    df = pd.read_csv(p)
    return df.set_index("FOV")[["x_global_px", "y_global_px"]].astype(float)


def stitch(unit_id, slide, fovs, dest, out_dir, P):
    umpx, FOVW, TARGET, Q = P["umpx"], P["fovw"], P["target_long"], P["quality"]
    pos = load_positions(dest, slide)
    ccdir = find_composite_dir(dest, slide)

    # keep only FOVs that have BOTH a position and an existing composite image
    have, missing = {}, []
    for fov in fovs:
        fn = os.path.join(ccdir, f"CellComposite_F{fov:05d}.jpg")
        if fov in pos.index and os.path.exists(fn):
            have[fov] = (float(pos.loc[fov, "x_global_px"]), float(pos.loc[fov, "y_global_px"]), fn)
        else:
            missing.append(fov)
    if not have:
        print(f"[skip] {unit_id}: 0 valid FOVs (missing {missing})")
        return None

    x0s = np.array([v[0] for v in have.values()])
    y0s = np.array([v[1] for v in have.values()])
    xmin = x0s.min(); xmax = x0s.max() + FOVW
    ymax = y0s.max(); ymin = y0s.min() - FOVW           # top / bottom edges (global px)
    span_x, span_y = xmax - xmin, ymax - ymin
    f = TARGET / max(span_x, span_y)
    W, H = int(round(span_x * f)), int(round(span_y * f))
    tw = int(round(FOVW * f))
    rawExtent = [xmin * umpx, ymin * umpx, xmax * umpx, ymax * umpx]  # UNCHANGED micron extent

    mos = Image.new("RGB", (W, H), (0, 0, 0))
    for fov, (x0, y0, fn) in sorted(have.items()):
        im = Image.open(fn).convert("RGB")
        if im.size != (FOVW, FOVW):
            print(f"[warn] {unit_id} FOV {fov} size {im.size} != {(FOVW, FOVW)}")
        imr = im.resize((tw, tw), Image.LANCZOS); im.close()
        px = int(round((x0 - xmin) * f))
        py = int(round((y0 - FOVW - ymin) * f))          # top-left of this FOV tile
        mos.paste(imr.transpose(Image.FLIP_TOP_BOTTOM), (px, py))  # per-FOV vertical flip

    # --- Y-FLIP FIX: flip the whole assembled mosaic to match the viewer cell flip ---
    mos = mos.transpose(Image.FLIP_TOP_BOTTOM)

    from io import BytesIO
    buf = BytesIO(); mos.save(buf, format="JPEG", quality=Q)
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    rawImage = "data:image/jpeg;base64," + b64

    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"rawimg_{unit_id}.json")
    json.dump({"rawExtent": rawExtent, "rawImage": rawImage,
               "meta": {"unit_id": unit_id, "slide": slide, "W": W, "H": H,
                        "n_fovs_used": len(have), "missing_fovs": missing,
                        "jpeg_bytes": buf.tell(), "flipud_applied": True}},
              open(out, "w"))
    print(f"[ok] {out}  W,H {W},{H}  fovs {len(have)}/{len(fovs)}  "
          f"jpeg {buf.tell()/1e3:.0f}KB  rawExtent {[round(v,1) for v in rawExtent]}"
          + (f"  missing {missing}" if missing else ""))
    return out


def main():
    ap = argparse.ArgumentParser(description="Stitch morphology mosaic for one (patient,slide) unit.")
    ap.add_argument("unit_id", help="unit_id from viewer_units.tsv (e.g. <unit> or <unit>)")
    ap.add_argument("--dest", default=DEFAULTS["dest"])
    ap.add_argument("--units", default=DEFAULTS["units"])
    ap.add_argument("--out", default=DEFAULTS["out"])
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    units_tsv = args.units or os.path.join(here, "viewer_units.tsv")
    out_dir = args.out or os.path.join(args.dest, "rawimg_build")

    man = pd.read_csv(units_tsv, sep="\t", dtype={"unit_id": str})
    row = man[man["unit_id"] == args.unit_id]
    if row.empty:
        sys.exit(f"unit_id '{args.unit_id}' not in {units_tsv}")
    row = row.iloc[0]
    fovs = [int(v) for v in str(row["fovs"]).split(",") if v != ""]
    stitch(args.unit_id, str(row["slide"]), fovs, args.dest, out_dir, DEFAULTS)


if __name__ == "__main__":
    main()
