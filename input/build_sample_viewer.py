#!/usr/bin/env python
# =============================================================================
# CosMx per-unit spatial-viewer generator  (LOCAL, token-free)
# =============================================================================
# Description : Parameterized generalization of extract_sample.py + template
#               injection. For a given patient (and slide) it builds the DATA
#               json from the clean AnnData object, injects it into the golden
#               viewer template, and (if available) embeds the per-unit RAW
#               morphology mosaic (rawImage/rawExtent from rawimg_<unit>.json).
# Unit model  : a "unit" = (patient_id, slide). Single-slide patients -> one
#               viewer  <prefix><UNIT>_sample_viewer.html ; multi-slide units
#               -> one per slide  <prefix><UNIT>_<SLIDE>_sample_viewer.html.
# Input       : the cell object .h5ad  (obs index = barcode;
#                                                          X = CSR log1p)
#               subtype/cgc_hierarchical_labels_v6.tsv.gz    (L2/L3/tumor_putative)
#               sample_viewer_template.html               (DATA -> __DATA__)
#               [rawimg_<unit>.json]  {rawImage, rawExtent}  (from HPC stitch)
# Output      : <prefix><unit>_sample_viewer.html   (one per unit)
#               viewer_units.tsv                (with --dump-manifest)
# Conda env   : <env>/bin/python  # e.g. a conda env named scvi
# Key deps    : h5py, scipy, pandas, numpy   (NO anndata: backed read trips on
#                                             a null-encoded /uns/log1p entry)
# Validated   : 2026-08-26  (test unit <unit> / TMA_1, cells-on-tissue verified)
# Parameters  : See DEFAULTS dict below
# =============================================================================

import os, sys, json, argparse, warnings, base64
import numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
warnings.filterwarnings("ignore")

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
# Output filenames are "<prefix><unit>_sample_viewer.html". The prefix is
# study-specific, so it comes from cohort_config.json ("naming": {"file_prefix"}).
_PREFIX = ((_COHORT.get("naming") or {}).get("file_prefix") or "sample_")


# -- DEFAULTS (all tunable paths / constants) ---------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = ("<project_root>/"
         "<data_dir>")
DEFAULTS = {
    "h5ad":       f"{_ROOT}/object.h5ad",                    # obs index==barcode; X=CSR log1p
    "labels":     f"{_ROOT}/subtype/cgc_hierarchical_labels_v6.tsv.gz",  # cell,L1,L2,L3,tumor_putative
    "fov_meta_qc": f"{os.path.dirname(_ROOT)}/fov_meta_qc.tsv",   # AUTHORITATIVE (slide,fov)->patient
    "template":   os.path.join(HERE, "sample_viewer_template.html"),  # golden shell (DATA->__DATA__)
    "rawimg_dir": HERE,                        # where rawimg_<unit>.json live (after scp from HPC)
    "outdir":     os.path.join(HERE, "viewer"),
    "manifest":   os.path.join(HERE, "viewer_units.tsv"),
}
# per-cell IF channels (obs Mean.* columns) exposed in the viewer
CHANNELS = {"PanCK": "Mean.PanCK", "CD45": "Mean.CD45", "DAPI": "Mean.DAPI",
            "Membrane": "Mean.Membrane", "G": "Mean.G"}
# Empirical calibration of x_slide_mm against the per-slide FOV-position files.
# TMA_2 carries a constant +320.341 um X offset; TMA_1 is already aligned.
RAW_X_OFFSET_UM = {"TMA_1": 0.0, "TMA_2": 320.34095}

# -- important regions: curated per-unit overlay regions, loaded from cohort_config.json
#    (keys "<unit>|<slide>" -> cx_world/cy_world/r_um/pdx). Study data, not code.
_IMPORTANT_REGIONS = {tuple(k.split("|", 1)): v
                for k, v in (_COHORT.get("important_regions") or {}).items()
                if not k.startswith("_")}

# genes exposed (intersected with var_names; order preserved; absent ones dropped)
GENE_LIST = ["EPCAM","KRT20","KRT7","KRT8","KRT18","KRT19","KRT5","KRT6A","KRT14",
             "TP63","UPK2","GATA3","FOXA1","PPARG","MKI67","TOP2A","VIM","PTPRC",
             "CD3D","CD8A","CD4","IL7R","FOXP3","NKG7","MS4A1","CD79A","MZB1",
             "JCHAIN","LYZ","CD68","C1QA","CD163","FCN1","TPSAB1","CPA3","DCN",
             "LUM","ACTA2","POSTN","MYH11","PECAM1","VWF","CLDN5","LYVE1"]
# patient clinical fields: obs column -> viewer label
PATIENT_FIELDS = {
    "patient_id":"Patient","patient_sex":"Sex","patient_age_at_init":"Age",
    "patient_grade_enroll":"Grade","patient_stage_enroll":"Stage",
    "treatment_timepoint":"Timepoint","patient_prior_therapies":"Prior therapies",
    "patient_concomitant_cis":"Concomitant CIS","patient_hg_recurrence":"HG recurrence",
    "patient_ttr_months":"TTR (months)","patient_variant_histology":"Variant histology",
    "patient_progression":"Progression","patient_outcome":"Outcome","tma_location":"TMA loc"}





# -- lightweight h5ad reader (obs + selective CSR) ----------------------------
class H5adLite:
    def __init__(self, path):
        self.f = h5py.File(path, "r")
        self.n = self.f["X"].attrs["shape"][0]
        self._idx = np.array(self.f["obs"][self.f["obs"].attrs.get("_index", "_index")]).astype(str)
        vi = self.f["var"].attrs.get("_index", "_index")
        self.var_names = np.array(self.f["var"][vi]).astype(str)
        self._var_pos = {g: i for i, g in enumerate(self.var_names)}
        self._csr = None  # lazily loaded

    @property
    def index(self):
        return self._idx

    def col(self, k):
        """Return obs column k as a numpy array (categorical-decoded)."""
        d = self.f["obs"][k]
        if isinstance(d, h5py.Group):  # categorical
            return np.array(d["categories"]).astype(str)[np.array(d["codes"])]
        return np.array(d)

    def has(self, k):
        return k in self.f["obs"]

    def patient_auth(self, fov_meta_qc_path):
        """AUTHORITATIVE per-cell patient_id via (slide,fov)->patient from fov_meta_qc.
        The object's obs.patient_id is OVER-ASSIGNED by the tma_location mapping
        (e.g. <unit> gets 57% of all cells / 327 FOVs but is really ~21 valid FOVs).
        Cached after first call."""
        if getattr(self, "_apid", None) is None:
            m = pd.read_csv(fov_meta_qc_path, sep="\t")
            m = m[m["patient_id"].notna() & (m["patient_id"].astype(str).str.lower() != "nan")]
            s2 = np.where(m["slide"].str.contains("TMA_1"), "TMA_1",
                          np.where(m["slide"].str.contains("TMA_2"), "TMA_2", "?"))
            fmap = {(s, int(fv)): str(p) for s, fv, p in zip(s2, m["fov"], m["patient_id"])}
            slide = pd.Series(self.col("slide")).astype(str).values
            fovn = pd.to_numeric(pd.Series(self.col("fov_num")), errors="coerce").values
            self._apid = np.array([fmap.get((s, int(fn)), "nan") if fn == fn else "nan"
                                   for s, fn in zip(slide, fovn)])
        return self._apid

    def _load_csr(self):
        if self._csr is None:
            X = self.f["X"]
            self._csr = csr_matrix(
                (np.array(X["data"]), np.array(X["indices"]), np.array(X["indptr"])),
                shape=tuple(X.attrs["shape"]))
        return self._csr

    def genes(self, rows, gene_names):
        """Dense (len(rows) x len(gene_names)) log1p expression for the given rows."""
        cols = [self._var_pos[g] for g in gene_names]
        sub = self._load_csr()[rows][:, cols]
        return np.asarray(sub.todense())

    def gene_sparse(self, rows):
        """ALL-gene sparse quantized log1p expression for the given rows.
        Per-gene max scaling -> UInt8 codes (1..255; 0 = not stored). CSR (row=cell).
        Returns JSON-ready dict: names / geneMax / n / base64 typed arrays."""
        sub = self._load_csr()[rows].tocsr()
        gmax = np.asarray(sub.max(axis=0).todense()).ravel().astype(np.float32)
        gmax[gmax == 0] = 1.0
        scale = 255.0 / gmax[sub.indices]
        codes = np.clip(np.rint(sub.data * scale), 1, 255).astype(np.uint8)
        b64 = lambda a: base64.b64encode(a.tobytes()).decode("ascii")
        return {"names": self.var_names.tolist(),
                "geneMax": [round(float(x), 4) for x in gmax],
                "n": int(sub.shape[0]),
                "indptr": b64(sub.indptr.astype(np.int32)),
                "indices": b64(sub.indices.astype(np.uint16)),
                "data": b64(codes)}


# -- helpers ------------------------------------------------------------------
def unit_id(patient, slide, multi):
    return f"{patient}_{slide}" if multi else str(patient)


def build_manifest(A, apid):
    pid = pd.Series(apid).astype(str)                 # AUTHORITATIVE (fov_meta_qc), not obs.patient_id
    keep = ~pid.isin(["nan", "NA", "None", ""])
    df = pd.DataFrame({
        "patient": pid[keep].values,
        "slide":   pd.Series(A.col("slide")).astype(str)[keep].values,
        "fov":     pd.to_numeric(pd.Series(A.col("fov_num"))[keep], errors="coerce").astype("Int64").values,
    })
    nslide = df.groupby("patient")["slide"].nunique()
    rows = []
    for (patient, slide), g in df.groupby(["patient", "slide"]):
        fovs = sorted(int(v) for v in pd.unique(g["fov"].dropna()))
        multi = nslide[patient] > 1
        rows.append({"unit_id": unit_id(patient, slide, multi), "patient": patient,
                     "slide": slide, "n_fovs": len(fovs), "n_cells": len(g),
                     "fovs": ",".join(str(v) for v in fovs)})
    return pd.DataFrame(rows).sort_values(["patient", "slide"]).reset_index(drop=True)


def build_data_for_unit(A, labels, patient, slide, multi, rawimg_dir, apid):
    sel = ((pd.Series(apid).astype(str) == patient) &                 # AUTHORITATIVE, not obs.patient_id
           (pd.Series(A.col("slide")).astype(str) == slide)).values
    idx = np.where(sel)[0]
    if idx.size == 0:
        return None
    bc = pd.Series(A.index[idx])                       # cell barcodes (join key)
    def oc(k):                                         # obs column, subset to this unit
        return pd.Series(A.col(k)[idx]) if A.has(k) else pd.Series([np.nan] * idx.size)

    L1  = bc.map(labels["L1"]).fillna("LowQuality").astype(str)   # finalized major tier (L1)
    L3  = bc.map(labels["L3"]).fillna("LowQuality").astype(str)
    L2  = bc.map(labels["L2"]).fillna("LowQuality").astype(str)
    TUM = pd.to_numeric(bc.map(labels["tumor_putative"]), errors="coerce").fillna(0).astype(int)
    cov = float(bc.isin(labels.index).mean())

    gx = pd.to_numeric(oc("x_slide_mm"), errors="coerce") * 1000   # global microns (verified)
    gy = pd.to_numeric(oc("y_slide_mm"), errors="coerce") * 1000
    gene_sparse = A.gene_sparse(idx)   # ALL panel genes, sparse+quantized

    num = lambda k: pd.to_numeric(oc(k), errors="coerce").fillna(0).round(1).tolist()
    o_first = lambda k: (str(oc(k).dropna().iloc[0]) if oc(k).notna().any() else "NA")
    pat = {lab: o_first(c) for c, lab in PATIENT_FIELDS.items()}
    fov = pd.to_numeric(oc("fov_num"), errors="coerce").fillna(0).astype(int)

    data = {
        "sample": f"{slide} patient {patient}", "n": int(idx.size),
        "fovs": int(fov.nunique()), "global": True, "patient": pat,
        "x": [round(float(v), 1) for v in gx], "y": [round(float(v), 1) for v in gy],
        "fov": fov.tolist(),
        "type": L1.tolist(),   # finalized L1 (was obs.major_celltype)
        "subtype": L3.tolist(), "subL2": L2.tolist(), "tumor": TUM.tolist(),
        "channels": {name: num(c) for name, c in CHANNELS.items()},
        "meta": {"pathdx": oc("path_dx").astype(str).tolist(), "nCount": num("nCount_RNA"),
                 "typeScore": [round(float(v), 3) for v in
                               pd.to_numeric(oc("major_celltype_score"), errors="coerce").fillna(0)]},
        "geneSparse": gene_sparse,   # all 6175 genes selectable (sparse UInt8)
    }

    # -- important region important region cell-index membership (finalize computes polygon after post-processing) --
    _region_cfg = _IMPORTANT_REGIONS.get((patient, slide))
    if _region_cfg is not None:
        gx_arr = np.asarray(gx, dtype=float)
        gy_arr = np.asarray(gy, dtype=float)
        d2 = (gx_arr - _region_cfg["cx_world"])**2 + (gy_arr - _region_cfg["cy_world"])**2
        inside = np.where(d2 <= _region_cfg["r_um"]**2)[0].tolist()
        data["_importantRegion"] = {"pdx": _region_cfg["pdx"], "cellIdx": inside,
                        "cx_world": _region_cfg["cx_world"], "cy_world": _region_cfg["cy_world"],
                        "r_um": _region_cfg["r_um"]}
        data["importantRegion"] = None  # placeholder — finalize_region_polygon.py fills this in
    else:
        data["importantRegion"] = None

    uid = unit_id(patient, slide, multi)
    rawimg_path = os.path.join(rawimg_dir, f"rawimg_{uid}.json")
    have_raw = os.path.exists(rawimg_path)
    if have_raw:
        rj = json.load(open(rawimg_path))
        source_extent = list(rj["rawExtent"])
        x_offset = RAW_X_OFFSET_UM.get(slide, 0.0)
        y_sum = float(np.nanmin(gy) + np.nanmax(gy))
        # stitch_one flips the mosaic top-to-bottom, while the viewer flips
        # cells around the sample-specific min(y)+max(y). Reflect the displayed
        # extent around the same center and apply the measured slide X offset.
        display_extent = [source_extent[0] + x_offset,
                          y_sum - source_extent[3],
                          source_extent[2] + x_offset,
                          y_sum - source_extent[1]]
        data = {"rawExtent": display_extent, "rawImage": rj["rawImage"], **data}  # first -> gate honored

    info = {"unit_id": uid, "n": data["n"], "fovs": data["fovs"], "genes": len(gene_sparse["names"]),
            "label_cov": round(cov, 3), "morphology": have_raw, "outcome": pat.get("Outcome", "NA")}
    return data, info


def write_viewer(template, data, out_path):
    html = template.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return len(html)


# -- main ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Build per-unit spatial viewer(s).")
    ap.add_argument("patient", nargs="?", help="patient id (<unit>) or unit_id (<unit>)")
    ap.add_argument("--slide", default=None, help="restrict to one slide (TMA_1 / TMA_2)")
    ap.add_argument("--all", action="store_true", help="build every unit (loads matrix once)")
    ap.add_argument("--dump-manifest", action="store_true", help="write viewer_units.tsv and exit")
    for k, v in DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), default=v)
    args = ap.parse_args()

    A = H5adLite(args.h5ad)
    apid = A.patient_auth(args.fov_meta_qc)      # authoritative per-cell patient (obs is over-assigned)
    man = build_manifest(A, apid)

    if args.dump_manifest:
        man.to_csv(args.manifest, sep="\t", index=False)
        print(f"[manifest] {args.manifest}  ({man.shape[0]} units, {man.patient.nunique()} patients)")
        print(man[["unit_id", "patient", "slide", "n_fovs", "n_cells"]].to_string(index=False))
        return

    if args.all:
        rows = man
    elif args.patient:
        if args.patient in set(man["unit_id"]):
            rows = man[man["unit_id"] == args.patient]
        elif args.patient in set(man["patient"]):
            rows = man[man["patient"] == args.patient]
        else:
            ap.error(f"'{args.patient}' not found as patient or unit_id. See --dump-manifest.")
    else:
        ap.error("give a patient/unit_id, use --all, or --dump-manifest")
    if args.slide:
        rows = rows[rows["slide"] == args.slide]
    if rows.empty:
        print("[skip] no matching units"); return

    labels = pd.read_csv(args.labels, sep="\t")
    labels["cell"] = labels["cell"].astype(str)
    labels = labels.set_index("cell")
    template = open(args.template, encoding="utf-8").read()
    os.makedirs(args.outdir, exist_ok=True)
    multi_map = man.groupby("patient")["slide"].nunique()

    n_ok = n_skip = 0
    for _, r in rows.iterrows():
        patient, slide = r["patient"], r["slide"]
        res = build_data_for_unit(A, labels, patient, slide, multi_map[patient] > 1, args.rawimg_dir, apid)
        if res is None:
            print(f"[skip] {patient}/{slide}: 0 cells"); n_skip += 1; continue
        data, info = res
        out_path = os.path.join(args.outdir, f"{_PREFIX}{info['unit_id']}_sample_viewer.html")
        mb = write_viewer(template, data, out_path) / 1e6
        print(f"[ok] {out_path}  ({mb:.2f} MB) | cells {info['n']} fovs {info['fovs']} "
              f"genes {info['genes']} label_cov {info['label_cov']} "
              f"morphology {'YES' if info['morphology'] else 'no(rgbpath fallback)'} | {info['outcome']}")
        n_ok += 1
    print(f"[done] {n_ok} viewer(s) written, {n_skip} skipped -> {args.outdir}")


if __name__ == "__main__":
    main()
