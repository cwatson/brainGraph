# brainGraph 0.52.1

## Minor changes
* `brainGraph_init` does a better job of dealing with subcortical gray matter data
* `boot_global` and `permute.group` updated to work with new version of `corr.matrix`
* `get.resid` now returns *data.table* of residuals with a *Study.ID* column
* `corr.matrix`: minor syntax change
* `get.resid`: fixed minor bug when `use.mean=FALSE` but *covars* has columns
    *mean.lh* and/or *mean.rh*
* `get.resid`: fixed bug when `use.mean=TRUE` (syntax error for RH vertices)
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


# brainGraph 0.48.0

2015-12-08

Initial CRAN acceptance
