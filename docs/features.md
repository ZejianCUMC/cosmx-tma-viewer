# Viewer features and how they are verified

`input/verify_viewer.py` opens each built viewer in headless Chrome and checks the items marked ✔. All checks use
real mouse events at the on-screen position of what is being tested.

| Feature | Behaviour | Check |
|---|---|---|
| Two synced panels | Left / right panel each pick any layer (morphology, H&E, cell types, genes, channels); pan/zoom is shared. | |
| Canvas sizing | The drawing buffer follows the panel's displayed size (ResizeObserver), so picking stays exact after late layout changes. | ✔ drawing size = displayed size |
| Hover picking | Mouse position is mapped into drawing coordinates; hit radius follows the dot size. | ✔ ≥95% of on-screen dot centres, zoomed in, largest dots |
| Hidden-by-default categories | Values in `hidden_default` are not drawn in any view and never hovered; cell-type legends can show them. | ✔ drawn dots = visible cells |
| Image layers | Morphology and H&E show no cell tooltip; core / region labels keep theirs. | ✔ no cell tooltip on image layers |
| Grouped legend | The level-3 legend is grouped under level-2 headers; clicking a header hides / shows the whole group. | ✔ every level-3 value listed, nothing overflows |
| Morphology | Every FOV of every core is imaged (overview), clipped to the core frame, same transform as the cells. | |
| Full resolution | Per-core full-resolution images, decoded only when zoomed in and on screen. | ✔ not loaded at fit view, loaded when zoomed |
| Core frames | Each frame contains all of its core's cells; optional shared row heights. | ✔ cells outside frames = 0 |
| Regions | Tier-styled outlines (solid / dashed / dotted), label placed clear of core labels, evidence tooltip, toggle. | ✔ labels inside the canvas, no overlap with core labels |
| H&E | One image per core fitted to its frame; rotate / flip per core, remembered in the browser. | |
| Editable palette | Click a legend swatch to recolour a category (remembered in the browser). | |
| Landing page | Cards read from the viewers themselves; per-core flags as badges + notes; search includes flag text. | |
