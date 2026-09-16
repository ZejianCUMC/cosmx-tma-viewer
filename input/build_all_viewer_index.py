#!/usr/bin/env python3
"""Build the local landing page for all metadata-backed viewers."""

import argparse
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
_PREFIX = ((_COHORT.get("naming") or {}).get("file_prefix") or "sample_")
_INDEX = ((_COHORT.get("naming") or {}).get("index_file") or "index_all_viewers.html")


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "viewer_units.tsv"
OUTDIR = HERE / "viewer"
OUTPUT = OUTDIR / _INDEX


# Per-unit core / treatment / TTR facts, read back OUT of the built viewers by
# extract_core_summary.py, so this page cannot drift from what they display.
# Optional: without the file the cards simply omit those chips.
CORE_SUMMARY = HERE / "viewer_core_summary.tsv"


def load_core_summary(path):
    if not Path(path).exists():
        return {}
    with Path(path).open(newline="", encoding="utf-8") as h:
        return {r["unit_id"]: r for r in csv.DictReader(h, delimiter="\t")}


def aggregate_by_patient(summary, unit_to_patient):
    """Roll per-unit facts up to the PATIENT. A core belongs to exactly one slide,
    so summing across a patient's units double-counts nothing; TTR is a
    patient-level field and is identical across a patient's units."""
    agg = {}
    for unit, s in summary.items():
        pat = unit_to_patient.get(unit)
        if pat is None:
            continue
        a = agg.setdefault(pat, {"n_core": 0, "n_frame": 0, "n_pre": 0, "n_post": 0,
                                 "ttr_months": "", "units": 0})
        for k in ("n_core", "n_frame", "n_pre", "n_post"):
            a[k] += int(s[k])
        a["units"] += 1
        if s.get("ttr_months"):
            a["ttr_months"] = s["ttr_months"]
    return agg


def core_text(s):
    """'N cores (M pieces)'. The pieces note appears only when a core covers more
    than one physically separate piece of tissue."""
    if not s:
        return ""
    n_core, n_frame = int(s["n_core"]), int(s["n_frame"])
    return (f"{n_core} core" + ("s" if n_core != 1 else "")
            + (f" ({n_frame} pieces)" if n_frame != n_core else ""))


def fact_chips(s, this_slide=None):
    """Patient-level pre/post-treatment core counts and TTR. `this_slide` adds a
    note when the patient spans several viewers, so a card does not overstate
    what it opens."""
    if not s:
        return "", ""
    pre, post = int(s["n_pre"]), int(s["n_post"])
    ttr = s.get("ttr_months", "")
    out = []
    if pre:
        out.append(f'<span class="fact pre">Pre {pre}</span>')
    if post:
        out.append(f'<span class="fact post">Post {post}</span>')
    if ttr:
        out.append(f'<span class="fact ttr" title="time to recurrence">TTR {ttr} mo</span>')
    if this_slide and int(s.get("units", 1)) > 1:
        out.append(f'<span class="fact slide" title="this patient has cores on more than one '
                   f'slide; this viewer shows only this one">{int(this_slide["n_core"])} on this slide</span>')
    srch = " ".join(filter(None, [core_text(s), "pre" if pre else "", "post" if post else "",
                                  "multicore" if int(s["n_core"]) > 1 else "",
                                  "split" if int(s["n_frame"]) != int(s["n_core"]) else ""]))
    return (f'<div class="facts">{"".join(out)}</div>' if out else ""), srch


def natural_key(unit_id):
    patient = unit_id.split("_", 1)[0]
    number = int("".join(ch for ch in patient if ch.isdigit()) or 0)
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


def main(outdir=None, output=None, summary_path=None, manifest=None, fov_meta=None):
    global OUTDIR, OUTPUT, CORE_SUMMARY, MANIFEST, FOVMETA
    if outdir: OUTDIR = Path(outdir)
    OUTPUT = Path(output) if output else OUTDIR / OUTPUT.name
    if summary_path: CORE_SUMMARY = Path(summary_path)
    if manifest: MANIFEST = Path(manifest)
    if fov_meta: FOVMETA = Path(fov_meta)
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = sorted(csv.DictReader(handle, delimiter="\t"), key=lambda r: natural_key(r["unit_id"]))

    expected = {f"{_PREFIX}{row['unit_id']}_sample_viewer.html" for row in rows}
    present = {path.name for path in OUTDIR.glob(_PREFIX + "*_sample_viewer.html")}
    if expected != present:
        raise RuntimeError(
            f"Viewer set mismatch; missing={sorted(expected - present)}, "
            f"unexpected={sorted(present - expected)}"
        )

    total_cells = sum(int(row["n_cells"]) for row in rows)
    total_fovs = sum(int(row["n_fovs"]) for row in rows)
    total_patients = len({row["patient"] for row in rows})
    meta = load_meta()
    summary = load_core_summary(CORE_SUMMARY)
    by_patient = aggregate_by_patient(summary, {r["unit_id"]: r["patient"] for r in rows})
    tot_cores = sum(int(v["n_core"]) for v in summary.values())
    tot_frames = sum(int(v["n_frame"]) for v in summary.values())
    tot_pre = sum(int(v["n_pre"]) for v in summary.values())
    tot_post = sum(int(v["n_post"]) for v in summary.values())
    cards = []
    for row in rows:
        unit = html.escape(row["unit_id"])
        patient = html.escape(row["patient"])
        slide = html.escape(row["slide"])
        filename = f"{_PREFIX}{row['unit_id']}_sample_viewer.html"
        ptags, psrch = unit_pathchips(row, meta)
        facts, fsrch = fact_chips(by_patient.get(row["patient"]), summary.get(row["unit_id"]))
        tlsb = '<span class="regionbadge" title="curated important region">\u2605 Region</span>' if row["unit_id"] in BADGE_UNITS else ''
        tsr = " region" if row["unit_id"] in BADGE_UNITS else ""
        cards.append(
            f'''<a class="card" href="{filename}" data-search="{unit.lower()} {patient.lower()} {slide.lower()} {psrch}{tsr} {fsrch}" data-slide="{slide}">
              <div class="card-top"><strong>{unit}</strong><span class="chip">{slide.replace('_', ' ')}</span>{tlsb}</div>
              <div class="patient">Patient {patient}</div>
              {ptags}
              {facts}
              <div class="counts"><span>{int(row['n_cells']):,} cells</span><span>{core_text(by_patient.get(row["patient"]))}</span></div>
            </a>'''
        )

    document = f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>Sample viewers — all metadata-backed units</title>
  <style>
    :root{{--bg:#f5f7fb;--panel:#fff;--ink:#172033;--muted:#667085;--line:#dfe4ec;--accent:#3157c8;--accent2:#e9efff}}
    *{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}}
    main{{max-width:1180px;margin:auto;padding:42px 24px 64px}} h1{{font-size:30px;margin:0 0 8px}} .subtitle{{color:var(--muted);margin:0 0 24px}}
    .summary{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-bottom:24px}}
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
    .facts{{display:flex;flex-wrap:wrap;gap:4px;margin:0 0 12px}}
    .fact{{font-size:10.5px;border-radius:5px;padding:2px 6px;font-weight:600;white-space:nowrap;background:#eef1f4;color:#48515e}}
    .fact.pre{{background:#e4f1fb;color:#1c5c86}} .fact.post{{background:#f6e0de;color:#8a2016}}
    .fact.ttr{{background:#efeaf7;color:#57458a}} .fact.slide{{background:#fff;color:#7a8392;border:1px dashed #cfd6e0;font-weight:500}}
    .note{{color:var(--muted);font-size:13px;margin-top:22px}} @media(max-width:650px){{.summary{{grid-template-columns:repeat(2,1fr)}}main{{padding:25px 15px}}}}
  </style>
</head>
<body><main>
  <h1>Sample viewers</h1>
  <p class="subtitle">All metadata-backed viewer units. Core, treatment and TTR figures are per <strong>patient</strong>; the cell count is for the viewer you open. Select a card to open its self-contained spatial viewer.</p>
  <section class="summary">
    <div class="stat"><b>{len(rows)}</b><span>viewer units</span></div>
    <div class="stat"><b>{total_patients}</b><span>patients</span></div>
    <div class="stat"><b>{total_fovs:,}</b><span>FOVs</span></div>
    <div class="stat"><b>{total_cells:,}</b><span>cells</span></div>
    <div class="stat"><b>{tot_cores}</b><span>cores{f" ({tot_frames} pieces)" if tot_frames != tot_cores else ""}</span></div>
    <div class="stat"><b>{tot_pre} / {tot_post}</b><span>pre- / post-treatment cores</span></div>
    <div class="stat"><b>{sum(1 for v in by_patient.values() if v["n_pre"] and v["n_post"])}</b><span>patients with paired pre + post</span></div>
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
    ap = argparse.ArgumentParser(description="Build the viewer landing page.")
    ap.add_argument("--viewer-dir", default=None, help="dir holding the built viewers")
    ap.add_argument("--out", default=None, help="output html path")
    ap.add_argument("--summary", default=None, help="viewer_core_summary.tsv")
    ap.add_argument("--manifest", default=None, help="viewer_units.tsv")
    ap.add_argument("--fov-meta", default=None, help="fov_meta_qc.tsv")
    a = ap.parse_args()
    main(a.viewer_dir, a.out, a.summary, a.manifest, a.fov_meta)
