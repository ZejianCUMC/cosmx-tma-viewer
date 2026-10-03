# Changelog

## 2.0.0 (2026-10)

New table-driven pipeline; the previous AnnData-based chain moved to `legacy/`.

- `build_viewer.py`: one viewer per unit directly from cells / cores / units / FOV-position tables, with optional
  count export, overview and full-resolution FOV tiles, H&E images and regions. One translation per core places
  cells, tiles and regions together, at true microns.
- Viewer template rebuilt from the current viewer and made config-driven (wording, colours, core groups, header
  fields, hidden categories, H&E default orientation, region styling); no study-specific text in code.
- Viewer behaviour:
  - hover picking fixed after late layout changes (canvas size sync, mouse mapped to drawing coordinates, hit radius
    follows dot size);
  - hidden-by-default categories hidden in every view and never hovered;
  - no cell tooltip on morphology / H&E layers;
  - level-3 legend grouped by level 2 with group hide / show;
  - every FOV imaged, clipped to the core frame; full-resolution per-core morphology loaded lazily on zoom;
  - every cell inside its core frame; optional shared frame heights per row;
  - multiple regions with tiers, label placement clear of core labels, evidence tooltips;
  - H&E panel with per-core orientation controls.
- `build_index.py`: landing page with per-core flags.
- `make_fov_tiles.py`: raw FOV images → viewer tiles.
- `verify_viewer.py`: headless-Chrome checks of the behaviours above.
- `make_demo_data.py` + `run_demo.sh`: synthetic end-to-end demo.

## 1.x (2026-09)

AnnData-based builder with post-processing chain (degap, regrid, core annotation, region polygon) — see `legacy/`.
