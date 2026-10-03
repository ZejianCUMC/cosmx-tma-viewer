#!/usr/bin/env python
# =============================================================================
# Produce a de-identified example viewer from a built one
# =============================================================================
# Description : Masks every clinical / patient attribute in a built viewer's
#               DATA with a placeholder, so the file can be shared as an example
#               of the interface. Structural fields (coordinates, FOV numbers,
#               cell-type labels, core geometry) are preserved so the UI still
#               demonstrates what it does.
#
#               WHAT THIS DOES *NOT* REMOVE: the expression matrix, the cell
#               coordinates and the embedded tissue morphology image are all
#               real measured data from a real sample. They carry no clinical
#               labels after this, but they are still sample-derived. Decide
#               separately whether that is acceptable for your audience.
# Input       : a built viewer HTML
# Output      : a de-identified copy
# Conda env   : any python3
# Validated   : 2026-09-16
# Parameters  : See DEFAULTS dict below
# =============================================================================
import argparse, io, json, re

DEFAULTS = {
    "mask": "**",
    # DATA.patient values are all clinical; every one is masked.
    # Per-core fields that are clinical attributes rather than geometry:
    "core_fields": ["path_dx", "stage", "grade", "block"],
    # Keep the treatment timepoint: it drives the frame colour that the example
    # exists to demonstrate, and "pre-" vs "post-treatment" identifies no one.
    # Pass --mask-treatment to blank it too.
    "sample_label": "example unit",
}


def data_span(txt):
    m = re.search(r"const DATA\s*=\s*", txt)
    s = m.end(); d = 0; i = s
    while i < len(txt):
        c = txt[i]
        if c == '{':
            d += 1
        elif c == '}':
            d -= 1
            if d == 0:
                return s, i + 1
        i += 1
    raise SystemExit("[err] could not locate the DATA object")


def main():
    ap = argparse.ArgumentParser(description="De-identify a built viewer.")
    ap.add_argument("src")
    ap.add_argument("out")
    ap.add_argument("--mask", default=DEFAULTS["mask"])
    ap.add_argument("--mask-treatment", action="store_true",
                    help="also blank the treatment timepoint (frames go neutral grey)")
    ap.add_argument("--template", default=None,
                    help="splice into this shell instead of keeping the source's")
    a = ap.parse_args()
    M = a.mask

    txt = io.open(a.src, encoding="utf-8").read()
    s, e = data_span(txt)
    D = json.loads(txt[s:e])
    report = []

    # 1. patient-level chips: mask every value, keep the keys so labels still show
    if D.get("patient"):
        report.append(f"patient fields masked: {', '.join(D['patient'])}")
        D["patient"] = {k: M for k in D["patient"]}

    # 2. the sample caption carries the unit and patient id
    slide = str(D.get("sample", "")).split()[0] if D.get("sample") else ""
    D["sample"] = f"{slide} · {DEFAULTS['sample_label']}".strip(" ·")
    report.append(f"sample caption -> {D['sample']!r}")

    # 3. per-core clinical attributes
    for c in D.get("cores", []):
        for f in DEFAULTS["core_fields"]:
            if f in c:
                c[f] = M
        if a.mask_treatment:
            c["tp"], c["tpLabel"] = "", M
        c["clin"] = "; ".join([M, c.get("tpLabel", M), M, M])
    if D.get("cores"):
        report.append(f"{len(D['cores'])} core records masked "
                      f"({'including' if a.mask_treatment else 'keeping'} treatment timepoint)")

    # 4. per-cell pathology annotation
    if D.get("meta", {}).get("pathdx"):
        n = len(set(D["meta"]["pathdx"]))
        D["meta"]["pathdx"] = [M] * len(D["meta"]["pathdx"])
        report.append(f"per-cell pathology annotation masked ({n} distinct value(s))")

    # 5. this viewer's region key naming, if the target shell expects the newer one
    if "tlsRegion" in D and a.template:
        D["importantRegion"] = D.pop("tlsRegion")
        report.append("region key renamed to importantRegion for the target shell")

    blob = json.dumps(D, separators=(",", ":"))
    if a.template:
        tpl = io.open(a.template, encoding="utf-8").read()
        assert tpl.count("__DATA__") == 1, "template lost its __DATA__ placeholder"
        out = tpl.replace("__DATA__", blob)
    else:
        out = txt[:s] + blob + txt[e:]
    io.open(a.out, "w", encoding="utf-8").write(out)

    for r in report:
        print("  " + r)
    print(f"[deid] {a.out}  ({len(out)/1e6:.2f} MB)")
    print("  NOTE: expression values, coordinates and the morphology image are "
          "unchanged and remain sample-derived data.")


if __name__ == "__main__":
    main()
