# brainGraph 1.4.0

2017-06-10

## Bug fix
* `plot_brainGraph_gui` had multiple issues and a few features have been changed:
    * Overall execution should be faster than in previous versions
    * *Lobe*, *neighborhood*, and *community* selection are now in "scrolled windows" instead of drop-down lists. Multiple selections can be made either by pressing `Ctrl` and clicking, or by holding `Shift` and moving the arrow keys
    * Fixed problem with vertex colors
    * When choosing to plot *neighborhoods*, you can color the vertices based on which neighborhood they belong to (useful if multiple vertices are selected)
* `plot_brainGraph`: now returns `NA` (instead of throwing an error) if the specified *subgraph* expression results in a network with 0 vertices.
* `gateway_coeff` returned an error if the number of communities equals 1; this has been fixed
* `edge_asymmetry` contained a bug if the input graph had only one contralateral connection (usually only encountered in the GUI with neighborhood plots); this has been fixed for those situations.
* `create_mats`:
    * Fixed bug for deterministic tractography when the user would like to normalize the matrices by *ROI size*.
    * Fixed bug for when `threshold.by='density'`. Previously, it would keep the top *X*% for *each* subject
* `mtpc`: fixed a bug that would incorrectly calculate `A.crit`

## Major changes
* `brainGraph_GLM`:
    * new function argument *level*, which now allows you to perform inference for graph-level measures
    * new argument *perms*, to specify the permutation set if you would like to use the same one for multiple graph/vertex measures (to improve the speed of the new `mtpc` function).
* `create_mats`
    * From v1.1.1 onward, all `A.norm.sub` matrices will be symmetrized, regardless of the value of `threshold.by` (previously they were only symmetrized if using `threshold.by='density'`). This should not pose a problem as the default (to take the *maximum* of the off-diagonal elements) is also the default when creating graphs in `igraph`.
    * From v1.2.0 onward, you can specify `mean` as the value for `threshold.by`. This will threshold the matrices such that a connection will be kept if `mean(A_ij) + 2*sd(A_ij) > threshold`, for each of the `threshold` values specified. See the references in the User Guide for examples in the literature.
    * From v1.3.0 onward, the `threshold.by` default name has been changed from `raw` to `consensus`, as this functionality is what is called "consensus-based" thresholding in the literature.
    * From v1.3.0, you may also specify `threshold.by='consistency'` to perform *consistency-based* thresholding. See Roberts et al., 2017.

## New functions
* `apply_thresholds`: threshold an additional set of matrices (e.g., FA-weighted matrices in DTI tractography) based on a set of matrices that have already been thresholded (e.g., streamline-weighted matrices in DTI tractography)
* `centr_betw_comm`: calculate vertex *communicability betweenness centrality* (Estrada et al., 2009). This requires that the package `expm` is installed.
* `communicability`: calculate network *communicability* (Estrada & Hatano, 2008)
* `make_empty_brainGraph`: this is not a new function, but rather was not exported in previous versions. Now it is accessible without the "triple-colon" operator (i.e., it is no longer necessary to call with `brainGraph:::make_empty_brainGraph`).
* `mtpc`: the *multi-threshold permutation correction (MTPC)* method for statistical inference of either vertex- or graph-level measures (Drakesmith et al., 2015)
* `symmetrize_mats`: symmetrize a connectivity matrix by either the *maximum*, *minimum*, or *average* of the off-diagonal elements. You may select one of these three as an argument to the function `create_mats`.
* `s_core`: calculate the *s-core* membership of a graph's vertices (Eidsaa & Almaas, 2013); the graphs will have vertex attributes called `s.core`. This is analogous to the *k-core* but for weighted networks. The vertex attribute for *k-core* has been changed from `coreness` to `k.core`.

## Minor changes
* `analysis_random_graphs`: no longer requires a *covars* argument
* `get.resid`: no longer requires a *covars* argument, as it was redundant
* `sim.rand.graph.par`: the argument *clustering* is no longer TRUE by default
* Some function arguments have been slightly modified to reflect the object type (for example, changing `g` to `g.list` if the function requires a *list* object as input).
* `set_brainGraph_attr`: no longer calculates the graph's *clique number*. This operation is usually fast but takes exceedingly long in dense graphs and graphs with more vertices (e.g., `craddock200`)


----
# brainGraph 1.0.0

2017-04-10

First *major* release; Fifth CRAN version

## Bug fix
* `plot_perm_diffs` previously didn't work with a low number of permutations, but now will work with any number
* `sim.rand.graph.par` previously didn't work with graphs lacking a `degree` vertex attribute
* Fixed problem with `plot_brainGraph_GUI` when plotting in the sagittal view for neighborhood graphs

## Major changes
* Multiple functions now run significantly faster after I updated the code to be more efficient
* `permute.group.auc` has been removed, and now `permute.group` accepts multiple densities and returns the same results. It can still take a single density for the old behavior
* The `lobe` and `network` vertex attributes are now *character* vectors
* `NBS` now handles more complex designs and contrasts through `brainGraph_GLM_design` and `brainGraph_GLM_fit`. The function arguments are different from previous versions
* `SPM` has been removed and is replaced by `brainGraph_GLM`
* Added atlas `craddock200` (with coordinates from `DPABI/DPARSF`)

## New functions
* `brainGraph_GLM`: replaces `SPM` and allows for more complex designs and contrasts
* `brainGraph_GLM_design`: function that creates a design matrix from a `data.table`
* `brainGraph_GLM_fit`: function that calculates the statistics from a design matrix and response vector
* `create_mats`: replaces `dti_create_mats` and adds functionality for resting-state fMRI data; also can create matrices that will have a specific graph density
* `gateway_coeff`: calculate the *gateway coefficient* (Vargas & Wahl, 2014); graphs will have vertex attributes `GC` or `GC.wt` (if weighted graph)
* `plot_brainGraph_multi`: function to write a PNG file of 3-panel brain graphs (see User Guide for example)

## Minor changes
* `efficiency` replaces `graph.efficiency`; the old function name is still accessible (but may be removed eventually)
* `set_brainGraph_attr` replaces `set.brainGraph.attributes`; the old function name is still accessible (but may be removed eventually)
* `part_coeff` replaces `part.coeff`
* All of the `rich.` functions have been renamed. The period/point/dot in each of those functions is replaced by the *underscore*. So, `rich.club.norm` is now `rich_club_norm`, etc.
* `set_vertex_color` and `set_edge_color` replace `color.vertices` and `color.edges` (these functions are not exported, in any case)
* `contract_brainGraph` replaces `graph.contract.brain`
* `make_ego_brainGraph` replaces `graph_neighborhood_multiple` (so it is a similar name to *igraph*'s function `make_ego_graph`)
* `write_brainnet` replaces `write.brainnet`
* In the GUI, ordering of vertices in circle plots now more closely reflect their anatomical position, being ordered by y- and x-coordinates (in addition to still being grouped by *lobe*)

----
# brainGraph 0.72.0

2016-10-10

Fourth CRAN version

## Bug fix
* `sim.rand.graph.clust` previously returned a list; now it correctly returns an
    `igraph` graph object
* `aop` and `loo`: regional contributions were calculated incorrectly (without
    an absolute value)
* `rich.club.norm`: changed the p-value calculation again; this shouldn't affect
    many results, particularly if N=1,000 (random graphs)
* `NBS`:
  * the `t.stat` edge attribute was, under certain situations, incorrectly
    assigning the values; this has been fixed in the latest version
  * fixed bug when permutations didn't result in any connected components
  * fixed bug w/ data randomization; the bug didn't seem to affect the results
* `SPM`:
  * the permutation p-values were previously incorrect; has been fixed
  * added an argument to remove `NA` values
* `vec.transform`: fixed bug which occurred when the input vector is the same
    number repeated (i.e., when `range(x) = 0`)

## Major changes
* `dti_create_mats`: new function argument `algo` can be used to specify either 'probabilistic' or 'deterministic'. In the case of the latter, when dividing streamline count by ROI size, you can supply absolute streamline counts with the `mat.thresh` argument.
* Changed instances of `.parallel` to `use.parallel`; also, added it as an
    argument to `set.brainGraph.attributes` to control all of the functions that
    it calls; also added the argument to `part.coeff` and
    `within_module_deg_z_score`
* Added atlases `aal2.94`, `aal2.120`, and `dosenbach160`
* `plot_brainGraph`: can now specify the orientation plane, hemisphere to plot,
    showing a legend, and a character string of logical expressions for plotting
    subgraphs (previously was in `plot_brainGraph_list`)

## New functions
* `auc_diff`: calculates the area-under-the-curve across densities for two
    groups
* `cor.diff.test`: calculates the significance of the difference between
    correlation coefficients
* `permute.group.auc`: does permutation testing across all densities, and
    returns the permutation distributions for the difference in AUC between two
    groups
* `rich.club.attrs`: give a graph attributes based on rich-club analysis

## Minor changes
* Removed the `x`, `y`, and `z` columns from the atlas data files; now only the
    MNI coordinates are used. This should simplify adding a personal atlas to use
    with the package
* Added a column, `name.full` to some of the atlas data files
* `NBS`:
  * New edge attribute `p`, the p-value for that specific connection
  * Returns the `p.init` value for record-keeping
* `brainGraph_init`: can now provide a `covars` data table if you want to subset
    certain variables yourself, or if the file is named differently from
    `covars.csv`
* `plot_brainGraph`: can now manually specify a subtitle;
* `plot_brainGraph_gui`:
  * Option for specifying maximum values for edge widths
* `plot_corr_mat`: color cells based on weighted community or network
* `plot_global`:
  * legend position is now "bottom" by default
  * can specify `xvar` to be either "density" or "threshold"; if the latter, the
    x-axis is reversed
  * If data has a `Study.ID` column, the `ggplot2` function `stat_smooth` is used
    and the statistic is based on a generalized additive model
* `plot_perm_diffs`: added argument `auc` for using the area-under-the-curve
    across densities
* `plot_rich_norm`:
  * Added argument `fdr` to choose whether or not to use FDR-adjusted p-values
  * Should work for more than 2 groups
  * Now works with multi-subject data; collapses by *Group* and plots the group mean
* `plot_vertex_measures`: can facet by different variables (e.g., lobe, community, network, etc.)
* `set.brainGraph.attributes`:
  * calculate graph `strength`, which is the mean of vertex strength (weighted networks)
  * Invert edge weights for distance-based measures
* `write.brainnet`:
  * Now allows for writing weighted adjacency matrices, using the `edge.wt`
    function argument
  * Can color vertices by multiple variables

---
# brainGraph 0.62.0

2016-04-22

Third CRAN version

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
* `SPM`: new function that replaces and improves upon both `group.graph.diffs`
    and `permute.vertex`
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
