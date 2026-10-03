#!/usr/bin/env python
# =============================================================================
# Viewer alignment QC — overlay cells on the embedded morphology (viewer-exact)
# =============================================================================
# Self-contained: reads ONLY the built <prefix><unit>_sample_viewer.html (DATA has
# x/y, rawExtent, rawImage, per-cell Mean.PanCK), reproduces the VIEWER's exact
# DATA.global Y-flip cell mapping, overlays cells on the (already flipud-fixed)
# embedded mosaic, and reports a blob-proof quantitative check:
#   corr(per-cell PanCK, local green) must be POSITIVE and epithelial cells must
#   sit on brighter tissue than non-epithelial -> morphology is viewer-correct.
# (A flipped/mirrored mosaic gives a NEGATIVE corr even though blob cores still
#  visually "look" covered — that was the original bug.)
#
# Usage: verify_overlay.py <viewer.html> [out.png]
# Env  : <env>/bin/python  # e.g. a conda env named scvi  (PIL, numpy)
# =============================================================================
import sys, json, base64, io
import numpy as np
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None

PAL = {'Epithelial':(255,20,60),'Fibroblast':(154,99,36),'Tcell':(67,99,216),
       'Myeloid':(245,130,49),'Plasma':(145,30,180),'Bcell':(60,180,75),
       'Endothelial':(240,50,230),'Mast':(66,212,244)}


def extract_data(html_path):
    txt = open(html_path, encoding="utf-8", errors="replace").read()
    i = txt.index("const DATA=") + len("const DATA=")
    depth = 0; instr = False; esc = False
    for j in range(i, len(txt)):
        c = txt[j]
        if instr:
            if esc: esc = False
            elif c == "\\": esc = True
            elif c == '"': instr = False
        else:
            if c == '"': instr = True
            elif c == "{": depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0: return json.loads(txt[i:j+1])
    raise ValueError("DATA not found")


def main():
    html = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else html.replace(".html", "_align.png")
    D = extract_data(html)
    if not D.get("rawImage"):
        sys.exit("no rawImage embedded -> nothing to verify")
    img = Image.open(io.BytesIO(base64.b64decode(D["rawImage"].split(",", 1)[1]))).convert("RGB")
    W, H = img.size
    e = D["rawExtent"]
    x = np.array(D["x"], float); y = np.array(D["y"], float)
    typ = np.array(D["type"]); panck = np.array(D["channels"]["PanCK"], float)

    # ---- viewer-exact DATA.global Y-flip cell mapping ----
    a, b = y.min(), y.max()
    col = (x - e[0]) / (e[2] - e[0]) * W
    row = ((a + b - y) - e[1]) / (e[3] - e[1]) * H

    # ---- quantitative, blob-proof check ----
    arr = np.asarray(img, float)
    ci = np.clip(np.round(col).astype(int), 0, W - 1)
    ri = np.clip(np.round(row).astype(int), 0, H - 1)
    green = arr[ri, ci, 1]; bright = arr[ri, ci].sum(1)
    cg = np.corrcoef(panck, green)[0, 1]
    epi = typ == "Epithelial"
    ratio = bright[epi].mean() / max(1.0, bright[~epi].mean())
    on = float((bright > 30).mean())
    verdict = "PASS (viewer-correct)" if (cg > 0.05 and ratio > 1.0 and on > 0.7) else "FAIL (check flip)"
    print(f"corr(PanCK, local green) = {cg:+.4f}   [want > 0]")
    print(f"epi/non-epi local brightness = {ratio:.3f}   [want > 1]")
    print(f"fraction of cells on tissue (bright>30) = {on:.3f}   [want > 0.7]")
    print(f"VERDICT: {verdict}")

    # ---- overlay render ----
    canvas = img.copy(); dr = ImageDraw.Draw(canvas)
    for cx, cy, t in zip(col, row, typ):
        if 0 <= cx < W and 0 <= cy < H:
            dr.ellipse([cx-1.3, cy-1.3, cx+1.3, cy+1.3], fill=PAL.get(t, (200,200,200)))
    canvas.save(out, quality=85)
    print(f"overlay -> {out}  ({W}x{H})")


if __name__ == "__main__":
    main()
