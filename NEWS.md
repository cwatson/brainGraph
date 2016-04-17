# brainGraph 0.61.0

2016-04-17

## Bug fix
* `rich.club.norm` had a bug in calculating the p-values. If you have already
    gone through the process of creating random graphs and the object `phi.norm`,
    you can fix with the following code: (add another loop if you have
    single-subject graphs, e.g. DTI data)

```
for (i in seq_along(groups)) {
  for (j in seq_along(densities)) {
    max.deg <- max(V(g[[i]][[j]])$degree)
    phi.norm[[i]][[j]]$p <- sapply(seq_len(max.deg), function(x)
        sum(phi.norm[[i]][[j]]$phi.rand[, x] >= phi.norm[[i]][[j]]$phi.orig[x]) / N)
  }
}
```

where `N` is the number of random graphs generated.
* `dti_create_mats`: there was a bug when *sub.thresh* equals 0; it would take
    matrix entries, even if they were below the *mat.thresh* values. This has
    been fixed. Argument checking has also been added.

## Major changes
* Now requires the package `RcppEigen` for fast linear model calculations;
    resulted in major speed improvements
* Now requires the package `permute` for the `NBS` function
* `group.graph.diffs`:
  * Uses the function `fastLmPure` from `RcppEigen` for speed/efficiency
  * Can specify multiple alternative hypotheses
  * Linear model specification is more limited now, though (on *TODO* list)
* Added data table for the `destrieux.scgm` atlas

## New functions
* `NBS`: implements the network-based statistic
* `analysis_random_graphs`: perform all the steps for getting *small-world*
    parameters and normalized *rich-club* coefficients and p-values
* `plot_global`: create a line plot across all densities of global graph
    measures in the same figure
* `vertex_spatial_dist`: calculates the mean edge distance for all edges of a
    given vertex

## Minor changes
* `dti_create_mats`: changed a few arguments
* `edge_spatial_dist`: re-named from `spatial.dist`
* `group.graph.diffs`: returns a graph w/ spatial coord's for plotting
* `plot_brainGraph_list`:
  * You can now specify a condition for removing vertices (e.g. `hemi == "R"`
    will keep only right hemisphere vertices; includes complex logical
    expressions (i.e., with multiple '&' and '|' conditions)
  * Vertex sizing and coloring is a bit more flexible
* New vertex attribute `Lp` (average path length for each vertex)
* `plot_brainGraph_gui`:
  * Added a checkbox for displaying a color legend or not
  * Can color vertices by weighted community membership
  * Added an *Other* option for adjusting edge widths by a custom attribute
  * More options for adjusting vertex sizes when the graph is weighted
  * Made the GUI window more compact to fit lower screen resolutions
* `plot_rich_norm`:
  * New argument `facet.by` to group the plots by either "density" (default) or
    "threshold" (for multi-subject, e.g. DTI data)
* `set.brainGraph.attributes`: New calculations for weighted graphs:
  * *Modularity* and community membership
  * *Participation coefficient* and *within-module degree z-score*
  * Vertex-level *transitivity*
  * Vertex-level *shortest path lengths*


---
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
