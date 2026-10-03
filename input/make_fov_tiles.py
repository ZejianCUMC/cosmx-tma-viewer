#!/usr/bin/env python
# =============================================================================
# cosmx-tma-viewer : FOV morphology tiles
# =============================================================================
# Description : Resize instrument FOV images (e.g. CosMx CellStatsDir/CellComposite/CellComposite_F00001.jpg, one
#               square image per FOV) into the tiles read by build_viewer.py: <slide>_F<fov:05d>.jpg.
#               Overview tiles: ~256 px (0.5 px/um for a 512-um FOV). Full-resolution zoom tiles: ~2048 px (0.25 um/px).
#               Run where the raw images live (e.g. a cluster node); copy only the tiles.
# Input       : --source DIR with images whose name contains F<5-digit fov>; --slide label; or --list TSV (slide, fov, path)
# Output      : <out>/<slide>_F<fov:05d>.jpg + manifest.tsv (slide, fov, source, bytes, md5)
# Key deps    : Pillow
# =============================================================================
import argparse, hashlib, re
from pathlib import Path
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

DEFAULTS = {"tile_px": 2048, "quality": 90, "fov_regex": r"F(\d{5})"}

def run_pipeline(out_dir, source=None, slide=None, list_tsv=None, fovs=None, **kw):
    P = {**DEFAULTS, **kw}; out = Path(out_dir); out.mkdir(parents=True, exist_ok=True); jobs = []
    if list_tsv:
        for line in Path(list_tsv).read_text().splitlines()[1:]:
            if line.strip(): s, f, p = line.split("\t")[:3]; jobs.append((s, int(f), Path(p)))
    else:
        for p in sorted(Path(source).iterdir()):
            m = re.search(P["fov_regex"], p.name)
            if m and (not fovs or int(m.group(1)) in fovs): jobs.append((slide, int(m.group(1)), p))
    rows = []
    for s, f, p in jobs:
        im = Image.open(p).convert("RGB"); assert im.size[0] == im.size[1], f"{p}: FOV image must be square"
        o = out / f"{s}_F{f:05d}.jpg"; im.resize((P["tile_px"],) * 2, Image.LANCZOS).save(o, quality=P["quality"])
        rows.append(f"{s}\t{f}\t{p.name}\t{o.stat().st_size}\t{hashlib.md5(o.read_bytes()).hexdigest()}")
    (out / "manifest.tsv").write_text("slide\tfov\tsource\tbytes\tmd5\n" + "\n".join(rows) + "\n")
    print(f"[ok] {len(rows)} tiles at {P['tile_px']} px -> {out}"); return {"output_path": str(out), "params_used": P}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="resize FOV images into viewer morphology tiles")
    ap.add_argument("out_dir"); ap.add_argument("--source"); ap.add_argument("--slide"); ap.add_argument("--list")
    ap.add_argument("--fovs", help="comma-separated FOV ids to include"); ap.add_argument("--tile-px", type=int, default=DEFAULTS["tile_px"])
    a = ap.parse_args()
    run_pipeline(a.out_dir, a.source, a.slide, a.list, {int(x) for x in a.fovs.split(",")} if a.fovs else None, tile_px=a.tile_px)
