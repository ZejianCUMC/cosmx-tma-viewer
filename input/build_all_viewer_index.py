#!/usr/bin/env python3
"""Build the local landing page for all metadata-backed CGC viewers."""

import csv
import html
import os
from pathlib import Path

# -- cohort config (study data lives OUTSIDE the code) ------------------------
# Unit ids, per-unit pathology labels and cohort counts are study data, not
# code, so they are read from cohort_config.json instead of being hardcoded.
# See cohort_config.example.json. A missing file or key disables the feature
# that uses it; nothing else changes.
def _cohort_cfg():
    import json as _json
    for _d in (os.path.dirname(os.path.abspath(__file__)),
               os.path.dirname(os.path.dirname(os.path.abspath(__file__))), os.getcwd()):
        _p = os.path.join(_d, "cohort_config.json")
        if os.path.exists(_p):
            try:
                return _json.load(open(_p, encoding="utf-8"))
            except Exception as _e:
                print(f"[warn] bad cohort_config.json ({_e}); continuing without it")
                return {}
    return {}

_COHORT = _cohort_cfg()


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "viewer_units.tsv"
OUTDIR = HERE / "viewer"
OUTPUT = OUTDIR / "index_all_CGC_viewers.html"


def natural_key(unit_id):
    patient = unit_id.split("_", 1)[0]
    number = int(patient.removeprefix("CGC"))
    return number, unit_id


FOVMETA = HERE.parent / "fov_meta_qc.tsv"

# Pathology-label -> (background, foreground) chip colours for the landing page.
# The label vocabulary is study-specific, so it lives in cohort_config.json
# ("index": {"path_colors": {...}}). Unlisted labels fall back to neutral grey.
_PATHCOL = {k: tuple(v) for k, v in
            ((_COHORT.get("index") or {}).get("path_colors") or {}).items()
            if not k.startswith("_")}
def _pcol(dx):
    return _PATHCOL.get(dx, ("#e5e8ec", "#5a6472"))

# Units given a region badge on the landing page; study data, see cohort_config.json.
BADGE_UNITS = {k: v for k, v in ((_COHORT.get("index") or {}).get("region_badge_units") or {}).items()
             if not k.startswith("_")}

def load_meta():
    with FOVMETA.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def unit_pathchips(row, meta):
    from collections import Counter
    tma = row["slide"].replace("_", "")
    fovs = set(row["fovs"].split(","))
    core_dx = {}
    for m in meta:
        if m.get("tma") == tma and m.get("fov") in fovs:
            dx = (m.get("path_dx") or "").strip()
            if dx and "image missing" not in dx:
                core_dx.setdefault(m.get("core_id"), dx)
    order = sorted(Counter(core_dx.values()).items(), key=lambda kv: (-kv[1], kv[0]))
    if not order:
        return "", ""
    chips = "".join(
        f'<span class="ptag" style="background:{_pcol(dx)[0]};color:{_pcol(dx)[1]}">{html.escape(dx)} {n}</span>'
        for dx, n in order)
    return f'<div class="pathtags">{chips}</div>', " ".join(dx.lower() for dx, _ in order)


def main():
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = sorted(csv.DictReader(handle, delimiter="\t"), key=lambda r: natural_key(r["unit_id"]))

    expected = {f"CGC_{row['unit_id']}_sample_viewer.html" for row in rows}
    present = {path.name for path in OUTDIR.glob("CGC_*_sample_viewer.html")}
    if expected != present:
        raise RuntimeError(
            f"Viewer set mismatch; missing={sorted(expected - present)}, "
            f"unexpected={sorted(present - expected)}"
        )

    total_cells = sum(int(row["n_cells"]) for row in rows)
    total_fovs = sum(int(row["n_fovs"]) for row in rows)
    total_patients = len({row["patient"] for row in rows})
    meta = load_meta()
    cards = []
    for row in rows:
        unit = html.escape(row["unit_id"])
        patient = html.escape(row["patient"])
        slide = html.escape(row["slide"])
        filename = f"CGC_{row['unit_id']}_sample_viewer.html"
        ptags, psrch = unit_pathchips(row, meta)
        tlsb = '<span class="regionbadge" title="curated important region">\u2605 Region</span>' if row["unit_id"] in BADGE_UNITS else ''
        tsr = " region" if row["unit_id"] in BADGE_UNITS else ""
        cards.append(
            f'''<a class="card" href="{filename}" data-search="{unit.lower()} {patient.lower()} {slide.lower()} {psrch}{tsr}" data-slide="{slide}">
              <div class="card-top"><strong>{unit}</strong><span class="chip">{slide.replace('_', ' ')}</span>{tlsb}</div>
              <div class="patient">Patient {patient}</div>
              {ptags}
              <div class="counts"><span>{int(row['n_cells']):,} cells</span><span>{int(row['n_fovs'])} FOVs</span></div>
            </a>'''
        )

    document = f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>CGC sample viewers — all metadata-backed units</title>
  <style>
    :root{{--bg:#f5f7fb;--panel:#fff;--ink:#172033;--muted:#667085;--line:#dfe4ec;--accent:#3157c8;--accent2:#e9efff}}
    *{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}}
    main{{max-width:1180px;margin:auto;padding:42px 24px 64px}} h1{{font-size:30px;margin:0 0 8px}} .subtitle{{color:var(--muted);margin:0 0 24px}}
    .summary{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin-bottom:24px}}
    .stat{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:15px 17px}} .stat b{{display:block;font-size:22px}} .stat span{{color:var(--muted)}}
    .tools{{display:flex;gap:10px;margin-bottom:20px;flex-wrap:wrap}} input{{flex:1;min-width:230px;border:1px solid var(--line);border-radius:10px;padding:11px 13px;font:inherit;background:white}}
    button{{border:1px solid var(--line);border-radius:10px;padding:10px 14px;background:white;color:var(--ink);cursor:pointer}} button.active{{background:var(--accent);border-color:var(--accent);color:white}}
    #grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(235px,1fr));gap:13px}}
    .card{{display:block;background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px;text-decoration:none;color:inherit;transition:.14s transform,.14s border-color,.14s box-shadow}}
    .card:hover{{transform:translateY(-2px);border-color:#9badde;box-shadow:0 7px 20px #20306018}} .card-top,.counts{{display:flex;justify-content:space-between;gap:12px;align-items:center}}
    .card strong{{font-size:18px;color:var(--accent)}} .chip{{font-size:12px;background:var(--accent2);border-radius:999px;padding:3px 8px;color:#294ba8}}
    .patient{{color:var(--muted);margin:7px 0 8px}} .counts{{font-size:13px;border-top:1px solid var(--line);padding-top:11px}} .hidden{{display:none}}
    .pathtags{{display:flex;flex-wrap:wrap;gap:4px;margin:0 0 12px}} .ptag{{font-size:10.5px;border-radius:5px;padding:2px 6px;font-weight:600;white-space:nowrap}}
    .regionbadge{{font-size:11px;background:#ffe08a;color:#8a5a00;border-radius:999px;padding:2px 8px;font-weight:700;margin-left:4px}}
    .note{{color:var(--muted);font-size:13px;margin-top:22px}} @media(max-width:650px){{.summary{{grid-template-columns:repeat(2,1fr)}}main{{padding:25px 15px}}}}
  </style>
</head>
<body><main>
  <h1>CGC sample viewers</h1>
  <p class="subtitle">All metadata-backed viewer units. Select a card to open its self-contained spatial viewer.</p>
  <section class="summary">
    <div class="stat"><b>{len(rows)}</b><span>viewer units</span></div>
    <div class="stat"><b>{total_patients}</b><span>patients</span></div>
    <div class="stat"><b>{total_fovs:,}</b><span>FOVs</span></div>
    <div class="stat"><b>{total_cells:,}</b><span>cells</span></div>
  </section>
  <div class="tools">
    <input id="search" type="search" placeholder="Search unit, patient, or slide…" aria-label="Search viewers">
    <button class="active" data-filter="all">All</button><button data-filter="TMA_1">TMA 1</button><button data-filter="TMA_2">TMA 2</button>
  </div>
  <section id="grid">{''.join(cards)}</section>
  <p class="note">Generated from <code>viewer_units.tsv</code>. The individual viewers include embedded morphology and require no network connection.</p>
</main>
<script>
  const cards=[...document.querySelectorAll('.card')], search=document.getElementById('search'); let slide='all';
  function apply(){{const q=search.value.trim().toLowerCase(); cards.forEach(c=>c.classList.toggle('hidden',!(c.dataset.search.includes(q)&&(slide==='all'||c.dataset.slide===slide))))}}
  search.addEventListener('input',apply); document.querySelectorAll('button[data-filter]').forEach(b=>b.addEventListener('click',()=>{{slide=b.dataset.filter;document.querySelectorAll('button[data-filter]').forEach(x=>x.classList.toggle('active',x===b));apply()}}));
</script></body></html>'''
    OUTPUT.write_text(document, encoding="utf-8")
    print(f"[done] {len(rows)} viewers, {total_patients} patients, {total_fovs} FOVs, {total_cells} cells")
    print(f"[index] {OUTPUT}")


if __name__ == "__main__":
    main()
