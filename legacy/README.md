# legacy/

The 1.x pipeline: `build_sample_viewer.py` (AnnData + labels → viewer) followed by the post-processing chain
(`assign_subcores.py` → `degap_viewer.py` → `regrid_by_core.py` → `inject_core_clinical.py` →
`finalize_region_polygon.py` → `verify_viewer_features.py`), the HPC morphology stitchers and the 1.x landing page.

It is kept for reference and is no longer maintained. New builds use `input/build_viewer.py` (2.x), which works from
plain tables, keeps every core at true microns and produces the current viewer behaviour (see `docs/features.md`).
Viewers built by 1.x use a different `DATA` layout and template; they are not compatible with the 2.x verifier.
