# contrib/

Optional helpers around the 2.0 pipeline in `input/`.

| file | what it does |
|---|---|
| `compare_viewers.py` | Regression check between two builds of the same viewer: per core, same cells (matched after a fitted per-axis scale + shift), labels, FOVs, decoded gene values and region contents; prints the fitted scale so geometry changes are visible. Exit code 1 on differences. |
| `slurm/make_fov_tiles.sbatch` | SLURM template for making full-resolution FOV tiles next to the raw images (fill the `<...>` placeholders). |
