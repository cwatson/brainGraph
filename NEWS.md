# brainGraph 0.51.0

## Minor changes
* Exported `plot_perm_diffs`
* Added argument checking for most functions
* `boot_global` added an OS check to get multicore functionality on Windows
* `brainGraph_init` now also returns the "tidied" dataset
* `count_interlobar` no longer takes `atlas.dt` as an argument
* `dti_create_mats` now accepts argument `P` for "number of samples"
* `edge_asymmetry` now works on Windows (changed from *mclapply* to *foreach*)
* `get.resid` got a complete overhaul; now works with *data.table* syntax
* `graph.efficiency` now works on Windows (changed from *mclapply* to *foreach*)
* `part.coeff` has a workaround to work on Windows
* `permute.group` no longer takes `atlas.dt` as an argument
* `vertex_attr_dt` is now essentially a wrapper for `igraph`'s function
    `as_data_frame`
