#!/usr/bin/env python
# =============================================================================
# cosmx-tma-viewer : landing page
# =============================================================================
# Description : One searchable card per built viewer. Every figure on a card (cores, cells, FOVs, per-core text) is
#               read back out of the viewer's own DATA block, so the page cannot drift from what the viewers show.
#               Optional per-core flags (flags.tsv: unit, core, label, detail, tooltip) add a badge on the card and a
#               note on the core; flagged cards are found by searching the flag text.
# Input       : directory of viewers built by build_viewer.py; optional flags.tsv; optional config (index.title/subtitle)
# Output      : <viewer_dir>/index.html (or --out)
# Key deps    : pandas
# =============================================================================
import argparse, html, json, re
from pathlib import Path
import pandas as pd

DEFAULTS = {"title": "CosMx TMA viewers", "subtitle": "", "file_glob": "*.html", "out_name": "index.html"}
CSS = """:root{--bg:#f5f7fb;--panel:#fff;--ink:#172033;--muted:#667085;--line:#dfe4ec;--accent:#3157c8;--accent2:#e9efff}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:1180px;margin:auto;padding:40px 22px 60px}h1{font-size:28px;margin:0 0 6px}.subtitle{color:var(--muted);margin:0 0 22px}
.summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-bottom:22px}
.stat{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px}.stat b{display:block;font-size:22px}.stat span{color:var(--muted)}
input{width:100%;border:1px solid var(--line);border-radius:10px;padding:11px 13px;font:inherit;background:#fff;margin-bottom:18px}
#grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(min(100%,270px),1fr));gap:13px;align-items:start}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px;min-width:0}
.card-top{display:flex;justify-content:space-between;gap:10px;align-items:center}.card h2{font-size:19px;margin:0}.card a{color:var(--accent)}
.chip{font-size:12px;background:var(--accent2);border-radius:999px;padding:3px 8px;color:#294ba8}
.facts{display:flex;flex-wrap:wrap;gap:4px;margin:10px 0}.fact{font-size:10.5px;border-radius:5px;padding:2px 6px;font-weight:600;background:#eef1f4;color:#48515e}
.counts{display:flex;justify-content:space-between;flex-wrap:wrap;font-size:13px;border-top:1px solid var(--line);padding-top:10px;margin-top:6px}
.flagbadge{display:inline-block;font-size:11px;background:#fde2e1;color:#9b1c1c;border:1px solid #f3b4b0;border-radius:999px;padding:2px 8px;font-weight:700}
.flagrow{margin:6px 0 0}.coreflag{color:#9b1c1c;font-weight:600}.open{display:inline-block;margin-top:10px;font-weight:600;text-decoration:none}
details{margin-top:10px;font-size:13px}.core{border-top:1px solid var(--line);padding:6px 0}.core p{margin:2px 0;color:var(--muted)}.hidden{display:none}"""

def read_data(fn):
    t = fn.read_text(encoding="utf-8"); m = re.search(r"const DATA\s*=\s*", t)
    return json.JSONDecoder().raw_decode(t, m.end())[0] if m else None

def run_pipeline(viewer_dir, out=None, flags_tsv=None, config=None):
    P = {**DEFAULTS, **((config or {}).get("index", {}))}; vd = Path(viewer_dir)
    flags = pd.read_csv(flags_tsv, sep="\t", dtype=str).fillna("") if flags_tsv else pd.DataFrame(columns=["unit", "core", "label", "detail", "tooltip"])
    out = Path(out or vd / P["out_name"]); cards = []; tot = {"viewers": 0, "cores": 0, "cells": 0}
    for fn in sorted(vd.glob(P["file_glob"])):
        if fn.resolve() == out.resolve(): continue
        D = read_data(fn)
        if not D or "cores" not in D: continue
        unit, slide = (D.get("sample", fn.stem).split(" · ") + [""])[:2]
        fu = flags[flags.unit == unit]; e = html.escape
        facts = "".join(f'<span class="fact">{e(k)} {e(str(v))}</span>' for k, v in D.get("patient", {}).items())
        cores_html = ""
        for c in D["cores"]:
            fc = fu[fu.core == c["id"]]
            badge = "".join(f' <span class="flagbadge" title="{e(r.tooltip)}">{e(r.label)}</span>' for r in fc.itertuples())
            note = "".join(f'<p class="coreflag">{e(r.detail)}</p>' for r in fc.itertuples() if r.detail)
            cores_html += (f'<div class="core"><b>{e(c["id"])}</b>{badge}<p>{e(c.get("short", ""))}</p>'
                           f'<p>{c["n"]:,} cells · {c["nfov"]} FOVs</p>{note}</div>')
        flagrow = "".join(f'<div class="flagrow"><span class="flagbadge" title="{e(r.tooltip)}">{e(r.label)} · {e(r.core)}</span></div>' for r in fu.itertuples())
        search = " ".join([unit, slide] + [str(v) for v in D.get("patient", {}).values()] + [c["id"] + " " + c.get("short", "") for c in D["cores"]]
                          + list(fu.label) + list(fu.detail)).lower()
        cards.append(f'<article class="card" data-search="{e(search)}"><div class="card-top"><h2><a href="{e(fn.name)}">{e(unit)}</a></h2>'
                     f'<span class="chip">{e(slide)}</span></div>{flagrow}<div class="facts">{facts}</div>'
                     f'<div class="counts"><span>{len(D["cores"])} core{"s" if len(D["cores"]) != 1 else ""}</span><span>{D["n"]:,} cells</span><span>{D.get("fovs", "")} FOVs</span></div>'
                     f'<a class="open" href="{e(fn.name)}">Open viewer →</a><details><summary>Core details</summary>{cores_html}</details></article>')
        tot["viewers"] += 1; tot["cores"] += len(D["cores"]); tot["cells"] += D["n"]
    stats = "".join(f'<div class="stat"><b>{v:,}</b><span>{k}</span></div>' for k, v in tot.items())
    page = (f'<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{html.escape(P["title"])}</title><style>{CSS}</style></head><body><main><h1>{html.escape(P["title"])}</h1>'
            f'<p class="subtitle">{html.escape(P["subtitle"])}</p><div class="summary">{stats}</div>'
            f'<input id="q" placeholder="search unit, core, annotation or flag…"><div id="grid">{"".join(cards)}</div></main>'
            '<script>const q=document.getElementById("q");q.addEventListener("input",()=>{const v=q.value.toLowerCase().trim();'
            'document.querySelectorAll("article.card").forEach(a=>a.classList.toggle("hidden",!!v&&!a.dataset.search.includes(v)))})</script></body></html>')
    out.write_text(page, encoding="utf-8"); print(f"[ok] {out}: {tot['viewers']} viewers, {tot['cores']} cores, {tot['cells']:,} cells, {len(flags)} flags")
    return {"output_path": str(out), "params_used": P, "summary": tot}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="build the landing page for a directory of viewers")
    ap.add_argument("viewer_dir"); ap.add_argument("--out"); ap.add_argument("--flags"); ap.add_argument("--config")
    a = ap.parse_args(); run_pipeline(a.viewer_dir, a.out, a.flags, json.load(open(a.config)) if a.config else None)
