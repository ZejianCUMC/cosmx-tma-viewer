#!/usr/bin/env python
# =============================================================================
# viewer — attach per-CORE clinical metadata to DATA.cores
# =============================================================================
# Description : Adds core-level clinical fields (path_dx, treatment timepoint,
#               grade, stage, block, tma_location) to every entry of DATA.cores
#               written by regrid_by_core.py, and splits the header annotation
#               into genuinely patient-level fields (DATA.patient) vs core-level
#               ones. The deployed viewers put ONE core's timepoint/tma_location
#               in the patient chip row (obs "first value"), which is wrong for
#               multi-core patients: <unit> shows "post_treatment / 2F" although
#               5 of its 6 mapped cores are pre_treatment.
# Input       : viewer HTML (post regrid_by_core.py; DATA.cores present)
#               <core_clinical>.xlsx  (Core_Clinical + Patient_Summary)
# Output      : same HTML rewritten in place (DATA span only)
# Conda env   : <env>/bin/python  # e.g. a conda env named scvi
# Key deps    : pandas, openpyxl
# Validated   : 2026-09-15 (<unit>; workbook cross-checked against fov_meta_qc)
# Parameters  : See DEFAULTS dict below
# =============================================================================

import argparse, json, re
import pandas as pd

# -- DEFAULTS -----------------------------------------------------------------
DEFAULTS = {
    # obs-derived chips that are really CORE-level -> removed from the patient row
    "drop_patient_keys": ["Timepoint", "TMA loc", "Grade", "Stage"],
    "tp_labels": {"pre_treatment": "Pre-treatment", "post_treatment": "Post-treatment"},
}
# core_id -> Core_Clinical column -> viewer key
CORE_FIELDS = {"path_dx": "path_dx", "treatment_timepoint": "tp", "grade_enroll": "grade",
               "stage_enroll": "stage", "block": "block", "tma_location": "loc"}


def extract_data(txt):
    m = re.search(r"const DATA\s*=\s*", txt)
    s = m.end(); d = 0; i = s
    while i < len(txt):
        c = txt[i]
        if c == '{':
            d += 1
        elif c == '}':
            d -= 1
            if d == 0:
                i += 1
                break
        i += 1
    return s, i, json.loads(txt[s:i])


def _s(v):
    return "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v).strip()


def process(path, xlsx, out=None, **kw):
    P = {**DEFAULTS, **{k: v for k, v in kw.items() if v is not None}}
    txt = open(path, encoding="utf-8").read()
    s, e, D = extract_data(txt)
    if not D.get("cores"):
        raise SystemExit("[err] DATA.cores missing — run regrid_by_core.py first")

    patient = str(D["patient"]["Patient"])
    cc = pd.read_excel(xlsx, "Core_Clinical")
    cc = cc[cc["patient_id"].astype(str) == patient].set_index("TMA_core_id")

    n_ok = 0
    for co in D["cores"]:
        if co["id"] not in cc.index:
            raise SystemExit(f"[err] core {co['id']} not in Core_Clinical for {patient}")
        row = cc.loc[co["id"]]
        for src, dst in CORE_FIELDS.items():
            co[dst] = _s(row.get(src))
        co["tpLabel"] = P["tp_labels"].get(co["tp"], co["tp"] or "Unknown")
        # the tag/tooltip string the viewer renders verbatim
        co["clin"] = "; ".join([co["path_dx"], co["tpLabel"], co["stage"], co["grade"]])
        n_ok += 1

    # listed-but-unmapped cores (present in the workbook, no cells in this unit)
    mapped = {co["id"] for co in D["cores"]}
    unmapped = [
        {"id": cid, "path_dx": _s(r.get("path_dx")),
         "tp": _s(r.get("treatment_timepoint")),
         "tpLabel": P["tp_labels"].get(_s(r.get("treatment_timepoint")), "Unknown"),
         "stage": _s(r.get("stage_enroll")), "grade": _s(r.get("grade_enroll")),
         "block": _s(r.get("block"))}
        for cid, r in cc.iterrows() if cid not in mapped]
    D.pop("coresUnmapped", None)   # listed-but-no-cells cores are NOT shown in the viewer

    # patient chip row: drop the fields that are really per-core
    for k in P["drop_patient_keys"]:
        D["patient"].pop(k, None)

    open(out or path, "w", encoding="utf-8").write(
        txt[:s] + json.dumps(D, separators=(",", ":")) + txt[e:])
    print(f"[clinical] {path.split('/')[-1]}: {n_ok} mapped cores annotated, "
          f"{len(unmapped)} listed-unmapped (not shown); dropped patient chips "
          f"{P['drop_patient_keys']}")
    for co in D["cores"]:
        print(f"    {co['id']:10s} [{co['clin']}]  block={co['block']}")
    for co in unmapped:
        print(f"    {co['id']:10s} (no cells in this unit) {co['path_dx']}; {co['tpLabel']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Attach per-core clinical metadata to DATA.cores.")
    ap.add_argument("path")
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    process(a.path, a.xlsx, a.out)
