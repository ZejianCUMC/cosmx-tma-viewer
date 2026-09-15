# Per-TMA-core annotation

Adds a per-core annotation layer to a built viewer: one frame per physical
tissue piece, coloured by treatment timepoint, with that core's own clinical
string in the header and a show/hide checkbox per core.

## Why it exists

The base pipeline annotates a viewer at the **unit** (patient, slide) level, and
takes core-level fields from the cell object's `obs` by "first value". For a unit
that contains several TMA cores this is wrong twice over:

- a field that varies per core (treatment timepoint, TMA location) gets one
  core's value applied to the whole unit;
- `regrid_viewer.py` groups tiles with a union-find heuristic over FOV
  centroids, which has no knowledge of `core_id`, so a tile can span two
  different cores.

This directory fixes both by making the authoritative FOV-to-core mapping the
basis for grouping, and by rendering every core-level field per core.

## What it adds to the viewer

1. **A frame per physical tissue piece**, with the `core_id` inline in the top
   border. The frame sits `frame_pad_um` outside the cell bounding box in world
   units and a further `dot radius + 2 px` at draw time, so it cannot touch a
   dot at any zoom.
2. **Frame colour by treatment timepoint** - configurable via
   `make_template_cores.py` DEFAULTS (`col_pre` / `col_post`). Pick values that
   do not collide with the categorical dot palette of the layer you view most.
3. **Hover a core label** -> a tooltip with that core's pathology, timepoint,
   stage, grade, block and cell / FOV counts.
4. **A per-core tag row** in the header: one tag per `core_id` carrying its own
   clinical string and a show/hide checkbox (default shown), plus All / None.
   Hiding a core removes its cells, its frames and its morphology tiles.
5. **Patient- vs core-level split.** Fields that are really per-core are dropped
   from the patient chip row and re-rendered as aggregates over the cores
   actually present, each with its core count, e.g. `Treatment  Pre(N cores),
   Post(M cores)`.

## One core_id can be two tissue pieces

A single TMA punch can appear as two physically separate tissue fragments
millimetres apart on the slide. Framing them as one box would draw a frame
across empty slide; grouping them by the old spatial heuristic alone could
instead merge FOVs belonging to two different cores.

`assign_subcores.py` does both correctly: it partitions FOVs by authoritative
`core_id` first, then splits each core into spatially connected components. A
component therefore never spans two cores and never merges two distant pieces.
Sub-cores of one punch **share the `core_id`** - same label, same clinical
record, same header tag, same checkbox. The `#k` suffix on the internal key
exists only so each piece gets its own frame.

**Component detection must run before `degap_viewer.py`.** Degap compacts every
empty band to roughly one FOV pitch, after which a multi-millimetre inter-piece
gap and a single tissue-free FOV hole are the same size and indistinguishable.

## Pipeline

```bash
PIPE=.
PY=python
V=work/<unit>_sample_viewer.html

$PY $PIPE/build_sample_viewer.py <unit> \
    --h5ad <object.h5ad> --labels <labels.tsv.gz> --fov-meta-qc <fov_meta_qc.tsv> \
    --template $PIPE/sample_viewer_template.html --rawimg-dir <rawimg_dir> --outdir work
$PY assign_subcores.py       $V --fov-meta-qc <fov_meta_qc.tsv>   # BEFORE degap
$PY $PIPE/degap_viewer.py    $V
$PY regrid_by_core.py        $V --ncol 3
$PY inject_core_clinical.py  $V --xlsx <core_clinical.xlsx>
$PY $PIPE/finalize_region_polygon.py $V
$PY make_template_cores.py --src $PIPE/sample_viewer_template.html
# splice: the new shell + the DATA span from $V  ->  the annotated viewer
```

`inject_core_clinical.py` expects a workbook with a `Core_Clinical` sheet keyed
by unit id + core id, carrying `path_dx`, `treatment_timepoint`, `grade_enroll`,
`stage_enroll`, `block` and `tma_location`. Validate that table against your FOV
metadata before trusting it: join its FOV sheet to `fov_meta_qc` on (slide, fov)
and require zero mismatches on core id, unit id, pathology and timepoint.

## Template edits

`make_template_cores.py` derives the annotated shell from the golden template by
applying edits anchored on exact substrings of the source, each asserted to match
exactly once, plus one whole-line replacement located in the source rather than
retyped - the minified `patinfo` builder, where a single space difference would
silently break matching. It then re-asserts every literal
`verify_viewer_features.py` greps for, so the UI gate cannot be broken silently.

Because only the `const DATA={...}` span differs between a built viewer and the
template, the annotated shell can be spliced onto viewers that were already
built, without rebuilding from the expression matrix. The exception is
`assign_subcores.py`, which needs pre-degap coordinates - adding core frames to
an existing viewer for the first time is a rebuild.

## Counts

Every number rendered comes from the viewer's own `DATA`, i.e. the cells that
survived QC - never from the clinical workbook's listed FOV counts, which
include FOVs dropped later in QC and will read high. Check where your object's
QC boundary actually sits before quoting a count: per-FOV QC flags and per-cell
quality labels are different axes, and a cell-level "low quality" class may
still be present in totals while hidden by default in the scatter.
