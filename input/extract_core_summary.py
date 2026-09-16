#!/usr/bin/env python
# =============================================================================
# Extract a per-unit core/treatment/TTR summary from built viewers
# =============================================================================
# Description : Reads each built viewer's DATA and emits one row per unit:
#               core_id count, frame (sub-core) count, pre/post-treatment core
#               counts and TTR. Sourced from the viewers themselves so the
#               landing page cannot drift from what the viewers display.
#               Uses a targeted bracket scan rather than parsing the whole DATA
#               object, which is dominated by the ~35 MB gene block.
# Input       : dir of <prefix><unit>_sample_viewer.html
# Output      : viewer_core_summary.tsv
# Conda env   : any python3 (no third-party dependencies)
# Validated   : 2026-09-16
# =============================================================================
import argparse, csv, glob, io, json, os, re

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


_PREFIX = ((_cohort_cfg().get("naming") or {}).get("file_prefix") or "sample_")



def slice_value(txt, key, open_ch, close_ch):
    """Return the JSON value for `key` by scanning matching brackets."""
    m = re.search(r'"%s":\s*\%s' % (re.escape(key), open_ch), txt)
    if not m:
        return None
    s = m.end() - 1
    d = 0
    i = s
    while i < len(txt):
        c = txt[i]
        if c == open_ch:
            d += 1
        elif c == close_ch:
            d -= 1
            if d == 0:
                return json.loads(txt[s:i + 1])
        i += 1
    return None


def main():
    ap = argparse.ArgumentParser(description="Extract per-unit core/treatment/TTR summary from built viewers.")
    ap.add_argument("viewer_dir")
    ap.add_argument("--prefix", default=_PREFIX)
    ap.add_argument("--out", default="viewer_core_summary.tsv")
    a = ap.parse_args()

    rows = []
    for p in sorted(glob.glob(os.path.join(a.viewer_dir, a.prefix + "*_sample_viewer.html"))):
        unit = os.path.basename(p)[len(a.prefix):-len("_sample_viewer.html")]
        txt = io.open(p, encoding="utf-8").read()
        cores = slice_value(txt, "cores", "[", "]") or []
        pat = slice_value(txt, "patient", "{", "}") or {}
        by_id = {}
        for c in cores:
            by_id.setdefault(c["id"], c)               # one entry per core_id
        pre = sum(1 for c in by_id.values() if c.get("tp") == "pre_treatment")
        post = sum(1 for c in by_id.values() if c.get("tp") == "post_treatment")
        ttr = str(pat.get("TTR (months)", "")).strip()
        try:
            ttr = f"{float(ttr):.1f}"
        except ValueError:
            ttr = ""                                    # NA / blank stays blank
        rows.append({"unit_id": unit, "n_core": len(by_id), "n_frame": len(cores),
                     "n_pre": pre, "n_post": post, "ttr_months": ttr,
                     "unknown_tp": len(by_id) - pre - post})
        print(f"  {unit:14s} {len(by_id)} cores / {len(cores)} pieces   "
              f"pre {pre} post {post}   TTR {ttr or 'NA'}")

    with io.open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\n[summary] {a.out}  ({len(rows)} units, "
          f"{sum(r['n_core'] for r in rows)} cores -> {sum(r['n_frame'] for r in rows)} pieces, "
          f"pre {sum(r['n_pre'] for r in rows)} / post {sum(r['n_post'] for r in rows)}, "
          f"{sum(1 for r in rows if r['unknown_tp'])} units with an unknown timepoint)")


if __name__ == "__main__":
    main()
