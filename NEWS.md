# brainGraph 0.55.0

2015-12-24

Second CRAN version

## New functions
* `aop` and `loo` calculate measures of *individual contribution* (see Reference
    within the function help)
  * Now requires the package `ade4`
* `plot_boot`: new function based on the removed plotting code from `boot_global`
* `plot_rich_norm`: function to plot normalized rich club coefficient curves

## Minor changes
* `boot_global`:
  * added an OS check to get multicore functionality on Windows
  * removed the code that created some plots
  * updated to work with the newer version of `corr.matrix`
* `brainGraph_init`:
  * does a better job of dealing with subcortical gray matter data
  * now also returns the "tidied" dataset
* `corr.matrix`:
  * was basically reverted back for speed purposes
  * minor syntax change
* `count_interlobar` no longer takes `atlas.dt` as an argument
* `dti_create_mats` now accepts argument `P` for "number of samples"
* `edge_asymmetry` now works on Windows (changed from *mclapply* to *foreach*)
* `get.resid`:
  * got a complete overhaul; now works with *data.table* syntax
  * now returns *data.table* of residuals with a *Study.ID* column
  * fixed minor bug when `use.mean=FALSE` but *covars* has columns
    *mean.lh* and/or *mean.rh*; fixed minor bug w/ RH residual calculation
  * fixed bug when `use.mean=TRUE` (syntax error for RH vertices)
* `graph.efficiency`: now works on Windows (changed from *mclapply* to *foreach*)
* `part.coeff`: has a workaround to work on Windows
* `permute.group`:
  * updated to work with new version of `corr.matrix`
  * no longer takes `atlas.dt` as an argument
* `vertex_attr_dt` is now essentially a wrapper for `igraph`'s function
    `as_data_frame`

* Exported `plot_perm_diffs`
* Added argument checking for most functions

---
# brainGraph 0.48.0

2015-12-08

Initial CRAN acceptance
