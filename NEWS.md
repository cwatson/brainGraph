# brainGraph 3.0.0

2020-09-28

## New functions/features
1. There are a few dozen new S3 methods for `bg_GLM` objects. See `methods(class='bg_GLM')` for the full list.
`coeff_determ` calculates the coefficient of determination.
`coeff_table` returns the coefficients table (same as `summary.lm(x)$coefficients`) for each region.
2. There are new GLM fitting functions (prefixed with `fastLmBG`) that are significantly faster and work with matrix/array inputs.
You can see these with the code `apropos('fastLm')`.
3. New functions `randomise` and `randomise_3d` can be called directly (although not recommended), and the `partition` function is now exported.
Each of these is for permutation-based analysis of linear models.
Furthermore, there are new permutation methods (`manly`, `draperStoneman`, and `stillWhite`).
4. New object `brainGraphList` for a collection of all graphs at a single density/threshold.
There are multiple S3 methods for this object, including the creation method `make_brainGraphList`.
5. `make_brainGraph` is now a S3 method.
6. There are several new matrix utility functions.
`inv` calculates the the "unscaled covariance" matrix used in linear models.
`pinv` calculates the *pseudoinverse*.
`qr` calculates the QR decomposition for each matrix in a 3D array.
`colMax`, `colMin`, and `colMaxAbs` calculate the max, min, and max of the absolute value across matrix columns.
`is_binary` determines if a matrix consists only of 0's and 1's.
`symmetrize` is now a S3 method. `symm_mean` symmetrizes a matrix using the mean of off-diagonal elements.
7. There are 4 new atlases: `hcp_mmp1.0` (HCP multimodal atlas), `power264`, `gordon333`, and `brainnetome`
8. New function `sim.rand.graph.hqs` generates random graphs from random covariance matrices for *structural covariance networks* using the HQS algorithm.
9. New plotting function `slicer` plots multiple graphs in a single figure.
10. Fewer package dependencies.
11. `mean_distance_wt` calculates weighted shortest path lengths.

## Removed/renamed functions
* `count_interlobar` is replaced by `count_inter`
* `make_mediate_brainGraph` is replaced by the `make_brainGraph` S3 method.
* `plot_brainGraph` is replaced by the `plot` S3 method for `brainGraph` objects.

# brainGraph 2.7.2

2019-10-20

## Bug fix
* The `mediation` package no longer exports `summary.mediate`, so it has to be removed from `brainGraph`
    * Move `mediation` to *Suggests*, as well


# brainGraph 2.7.1

2019-06-29

## Bug fix
* Fix bug in `import_scn` so the *Study.ID* column is always read as `character`
* Remove deprecated function `brainGraph_init`


# brainGraph 2.7.0

2018-12-15

## New functions/features
* `import_scn` replaces `brainGraph_init`, with a few changes in behavior:
    * It is no longer necessary to abbreviate region names yourself; the function does it automatically
    * Expects files with the name `${parcellation}_${hemi}_${modality}.csv` in the `datadir`
        - Here, `${parcellation}` could be `aparc`, for example
        - Also, `${modality}` could be `thickness`
    * If the *atlas* you are using includes `scgm`, there should be a `asegstats.csv` file
* `aop` and `loo` now return *S3* objects, with class name `IC`
    * These also have `summary` and `plot` methods
    * Furthermore, these objects return some more information

## Minor changes
* The `exclude` argument to `get.resid` is now `exclude.cov` to highlight that it is for specifying *covariates* to exclude from the GLM

# brainGraph 2.6.1

2018-12-07

## Bug fix
* Fix bug in `partition`, changing `method` to `part.method`
* Add some checks in `robustness` so it doesn't throw an error


# brainGraph 2.6.0

2018-09-04

## Bug fix
* Fixed bug in `count_homologous` that affected a subset of atlases
    * Performance is modestly improved (approx. 2-3x faster)

## New functions/features
* `count_inter` replaces `count_interlobar`; the new function calculates inter- and intra-group edge counts, where a group of vertices can be a *lobe*, *hemisphere*, *network* (for the `dosenbach160` atlas), or *class* (for the `destrieux` atlases)
    * The return object is now more informative; the function returns a matrix of all inter- and intra-group edge counts, in addition to a data.table containing a summary (that matches the output of previous versions)
* `rich_core` now calculates the rich core for weighted networks
    * In addition, the function runs *considerably* faster:
        - for smaller sparse graphs, it is ~40-80x faster
        - for larger dense graphs, it is more than 2,000x faster
* `robustness` now returns a data.table containing much more information (in addition to the max. connected component ratio)
    * This makes plotting outputs simpler; see Chapter 14 of the *User Guide*
    * When `type='edge'`, the function is about 2-3x faster than previous versions

## Minor changes
* `plot.mtpc`: the stats displayed in the caption have been "transposed", such that `S.crit` and `A.crit` are in the top row


# brainGraph 2.5.0

2018-09-01

## Bug fix
* Fixed regression bug in `NBS` (introduced by `v2.0.0`) which only occurred if `alternative='less'` when calculating the minimum statistic of permuted networks

## Minor changes
* Updated code that symmetrizes matrices:
    * Introduced new function, `symm_mean`, that more quickly symmetrizes a matrix about the diagonal by assigning `mean(c(A[i, j], A[j, i]))` to the off-diagonal elements
    * Uses `pmin` and `pmax` for symmetrizing matrices based on the off-diagonal minimum and maximum, respectively
* Optimized code in a few functions for faster execution:
    * `sim.rand.graph.clust` is about 2x faster due to improvement in the `choose.edges` helper function
    * `centr_lev` and `edge_asymmetry` are also faster


# brainGraph 2.4.0

2018-07-21

## New features
* `hubness`: new function for determining which vertices are hubs
* `set_brainGraph_attr`
    * New argument `clust.method` lets the user choose which clustering (community detection) method to use.
        1. The default is still the `louvain` algorithm.
        2. If you select `spinglass`, but the graph is unconnected, then `louvain` is used instead.
        3. If there are any negative edge weights, and you choose anything other than `walktrap` or `spinglass`, the `walktrap`  method is used.
    * Now calculates `num.hubs` using the new `hubness` function, and calculates separate values for weighted and unweighted networks


# brainGraph 2.3.4

2018-07-06

## Bug fix
* Fixed bugs in `rich_club_norm` that would throw an error if certain graph attributes weren't present

## New features
* `rich_club_all` - new function that is a wrapper for `rich_club_coeff`, applying over all possible degree values


# brainGraph 2.3.3

2018-06-25

## Bug fix
* Fixed regression bug in `plot.brainGraph`, which occurred when choosing `plane='sagittal'`


# brainGraph 2.3.2

2018-06-22

## Bug fix
* Fixed regression bug when fitting GLM models with a `F contrast`
* Fixed minor bug in `make_nbs_brainGraph` which did not properly assign the `p.nbs` attribute to all vertices

## Minor change
* The elements of the `NBS` output object, `p.mat` and `T.mat`, are now 3-dimensional arrays (with extent equal to the # of contrasts) instead of lists of matrices

# brainGraph 2.3.1

2018-06-20

## Bug fix
* Fixed a bug in `brainGraph_permute` that I didn't catch before

# brainGraph 2.3.0

2018-06-20

## Bug fix
* `brainGraph_boot` and `corr.matrix`:
    * Incorrectly calculated `E.global.wt` before; now it transforms edge weights
    * To do so, includes argument `xfm.type`
    * Fixed bug when calling `corr.matrix` (added `rand=TRUE`)
    * Also had to update the return object of `corr.matrix` for this purpose
* `mtpc`
    * Previously gave some incorrect results when `alt='less'`; fixed
    * The `plot` method also now gives correct values when `alt='less'`
* `brainGraph_GLM` now returns the correct *null.thresh* when `alt != 'greater'`
* `plot.brainGraph`: fixed bug that occurred when `plane='sagittal'` and a `hemi` value was not supplied
* `plot_rich_norm`: didn't plot values for all degrees present in the networks under certain scenarios

## New features
* `set_brainGraph_attr` now calculates a graph-level `Lp.wt`, which equals:
```
Lpv.wt <- distances(g)
Lpv.wt[is.infinite(Lpv.wt)] <- NA
g$Lp.wt <- mean(Lpv.wt[upper.tri(Lpv.wt)], na.rm=T)
```
* `plot_rich_norm`: new argument `smooth` lets you plot with a smoother in the case of single-subject data, as opposed to the previous default of a line plot for all subjects

## Minor changes
* GLM-related and other functions will now:
    1. Allow for the `Study.ID` column to be *numeric*; they will convert it to class *character*
    2. Creates a *character* vector of integers if `Study.ID` is not present in the data
* The `summary.mtpc` method now also prints the value of `clust.size`


# brainGraph 2.2.0

2018-05-28

## Minor changes
* Moved `RGtk2` and `cairoDevice` to *Suggests* (from *Depends*) to allow installation on headless servers
    * Thanks to `@michaelhallquist` for the pull request!
* Functions `boot_global`, `permute.group`, and `plot_group_means` are no longer accessible (deprecated since *v2.0.0*)

# brainGraph 2.1.0

2018-05-03 (mostly changes to *structural covariance network* functionality)

## Bug fix
* Fixed a bug in `mtpc` that was introduced in `v2.0.1`

## New functions/features
* `brainGraph_GLM_design` has a new argument `factorize` which specifies whether or not to convert all character columns (excluding *Study.ID*) to factor variables. The default is `TRUE`. Previously, character columns were ignored.
* `get.resid`
    * In the function call, you can choose whether or not to test a linear model for all groups together or separately, via the `method` argument
    * The `plot` method now returns a *list* of *ggplot* objects, and is similar to the `plot` methods for `bg_GLM` and `mtpc`
* `corr.matrix`
    * The `resids` argument must be the output of `get.resid` (not a *data.table* as before)
    * Correlations will be calculated separately for all subject groups (as this information is stored in the output of `get.resid`); you no longer need to loop (or `lapply`) across groups
    * In the function call, you can choose whether to correlate the residuals or raw structural values, via the `what` argument
    * The `exclusions` argument was renamed to `exclude.reg` to highlight that you should specify *region names* to be excluded (if any)
    * You can explicitly choose whether to calculate Pearson or Spearman correlations, via the `type` argument (previously, this behavior was "hidden")

## Minor changes
* `brainGraph_init`: the `modality` argument now will accept *any* character string; the default is still *thickness*. The files with the string you supply still must be present on your system.
* Due to `corr.matrix` expecting different input, the following functions also require, for their `resids` argument, the output of `get.resid` (instead of a *data.table*):
    * `aop`
    * `brainGraph_boot`
    * `brainGraph_permute`
    * `loo`


# brainGraph 2.0.4

2018-04-28

## Bug fix
* `gateway_coeff`: no longer throws an error for very sparse graphs; instead, it returns a vector with `NaN` values for unconnected vertices
* `make_mediate_brainGraph`: did not return correct values (for the treatment condition) when `INT=TRUE` (it recycled the values for the control condition)
* `make_intersection_brainGraph`
    * Previously exited with error if any of the input graphs did not contain vertices meeting the desired `subgraph` condition
    * Now returns an empty graph if none of the input graphs meet the `subgraph` condition
* `NBS`:
    * When getting the indices for which matrix elements to transpose (so that result is symmetric), the result was slightly wrong for `alt='greater'`
    * Calculation of edge counts in `summary` method contained an error

## Minor changes
* All `summary` methods now provide a `DT.sum` element in the returned list; previously it was inconsistent


# brainGraph 2.0.3

2018-04-26

## Bug fix
* In `mtpc`, the stats table that is returned previously was not always unique
* `mtpc` did not return a list with a named element `clust.size` (it was unnamed)
* In `plot.mtpc`, if the user selected a contrast other than the first, it would not plot the correct null statistics (green dots)


# brainGraph 2.0.2

2018-02-23

Release on CRAN; bugfix release.

## Bug fix
* Fixed a bug in `create_mats` in which the ordering (along the 3rd dimension) of the arrays in `A.norm.sub` did not match the ordering of the input matrix files (and therefore the ordering along the 3rd dimension of the arrays `A` and `A.norm`).
    * In the case that the input matrix files were already ordered by *Group* and *Study.ID*, then this is not a "bug", in that the ordering was already correct. So, if your subject groups are `groups <- c('Control', 'Patient')`, and the matrix files are separated on the filesystem by group, there is no change in behavior.
    * This bug only appeared when `threshold.by='consistency'` or `threshold.by='consensus'` (the default option).


# brainGraph 2.0.1

2018-02-07

## Bug fix
* Fixed error in `mtpc` when creating the MTPC statistics `data.table`


# brainGraph 2.0.0

2018-02-05

*2nd major release; 6th CRAN release*. (The previous CRAN release was at v1.0.0)

For other updates and bug fixes, see the minor release notes below.

## New functions/features
1. Mediation analysis is now possible through `brainGraph_mediate`.
2. I have introduced some simple *S3 classes* and *methods*. All of the classes have `plot` (except `NBS`) and `summary` methods.
The classes and corresponding "creation functions" are:

| Class             | Creation func.            | Description                       |
| -----             | -----                     | -----                             |
| brainGraph        | make_brainGraph           | Any graph with certain attributes |
| bg_GLM            | brainGraph_GLM            | Results of GLM analysis           |
| NBS               | NBS                       | Results of NBS analysis           |
| mtpc              | mtpc                      | Results of MTPC analysis          |
| brainGraph_GLM    | make_glm_brainGraph       | Graphs from GLM analysis          |
| brainGraph_NBS    | make_nbs_brainGraph       | Graphs from NBS analysis          |
| brainGraph_mtpc   | make_glm_brainGraph       | Graphs from MTPC analysis         |
| brainGraph_mediate| make_mediate_brainGraph   | Graphs from mediation analysis    |
| brainGraph_boot   | brainGraph_boot           | Results of bootstrap analysis     |
| brainGraph_permute| brainGraph_permute        | Results of permutation tests      |
| brainGraph_resids | get.resid                 | Residuals for covariance networks |

3. Multiple contrasts (in the same function call), as well as F-contrasts, are now allowed in the GLM-based functions: `brainGraph_GLM`, `mtpc`, `NBS`, and `get.resid`.
    * There is a new function argument, `con.type`, for this purpose.
    * Since both contrast types are now specified in the form of a *contrast matrix*, the argument `con.vec` has been replaced by `con.mat`.
4. Designs with 3-way interactions (e.g., `2 x 2 x 2`) are now allowed for GLM-based analyses.
5. Permutations for GLM-based analyses are now done using the *Freedman-Lane* method (the same as in FSL's *randomise* and in *PALM*).
6. Plot the "diagnostics" from GLM analyses through the `plot.bg_GLM` method to the output of `brainGraph_GLM`.
7. Plot the statistics from MTPC analyses through the `plot.mtpc` method for `mtpc` results.
8. `aop` has a new argument `control.value` allowing you to specify the control group; all comparisons will be to that group.
    * Removes the need to loop through patient groups in the console (if you have more than 1).
9. Most of the GLM-based functions have a new argument, `long`, which will not return all of the permutation results if `long=FALSE`.

## Removed/renamed functions
* `boot_global` was renamed to `brainGraph_boot`.
* `check.resid` was removed; you now just call the `plot` method to outputs of `get.resid`.
* `permute.group`:
    1. Function was renamed to `brainGraph_permute`.
    2. The arguments are slightly re-ordered
    3. Argument `permSet` was renamed to `perms`.
    4. New argument `auc` lets you explicitly define whether or not you want statistics for the *area under the curve (AUC)*.
* `plot_boot` was removed; you now just call the `plot` method to outputs of `brainGraph_boot`.
* `plot_brainGraph_mni` has been removed; this functionality can be changed by the `mni` argument to `plot.brainGraph` (i.e., the *plot method* for objects of class `brainGraph`)
* `plot_group_means` was renamed to `plot_volumetric`, as it works specifically for structural covariance networks.
* `plot_perm_diffs` was removed; you now just call the `plot` method to outputs of `brainGraph_permute`.

## Major changes
* `NBS` now automatically symmetrizes the input matrices. This is partly for speed and partly because `igraph` symmetrizes the matrices anyway.
    * There is a new function argument, `symm.by` (which is the same as that for `create_mats`) for this purpose.
* `corr.matrix`:
    * Now expects as its first input the residuals from `get.resid`.
    * You may specify multiple `densities` (or `thresholds`),
    * Returns a *list* including the binarized, thresholded matrices as an *array* (still named `r.thresh`).
* `get.resid` now allows for any design matrix for getting LM residuals (similar to `brainGraph_GLM`).
    * Must supply a `data.table` of covariates.
    * You may pass on arguments to `brainGraph_GLM_design` for creating the correct design matrix.
* `mtpc` accepts 2 new arguments (in addition to explicitly naming required arguments that pass on to `brainGraph_GLM`):
    1. `clust.size` lets you change the "cluster size", the number of consecutive thresholds needed to deem a result significant (default: `3`)
    2. `res.glm` lets you input the `res.glm` list element from a previous `mtpc` run. This is only useful if you would like to compare results with different values for `clust.size`.
* `permute.group` (see above section for changes)
* `rich_club_norm` now returns a `data.table`, which simplifies working with the data (and plotting).
* `set_brainGraph_attr`: multiple (explicit) arguments were removed; these are now passed on to `make_brainGraph` and can still be specified in the function call.
* I now use the `ggrepel` package for any `ggplot` objects with text labels.

# brainGraph 1.6.0

2017-09-14

## Bug fix
* `brainGraph_init`: fixed bug regarding the use of a custom atlas

## Minor changes
* Some function arguments have been modified to reflect the object type (e.g., changing `g` to `g.list` if the function requires a *list* object).
* `brainGraph_init`:
    * New argument `custom.atlas` allows you to use an atlas that is not in the package (you must also specify `atlas="custom"`).
    * This requires that the atlas you specify already be loaded into the R environment and meet the specifications of the package's atlases
        * It should be a `data.table`, and have columns *name*, *x.mni*, *y.mni*, *z.mni*, *lobe*, *hemi* (at a minimum).
* `permute.group`: can now calculate `ev.cent`


# brainGraph 1.5.0

2017-08-31

## Bug fix
* `boot_global`: fixed bug in *modularity* calculation

## Major changes
* `boot_global`:
    * can omit display of the progress bar (by setting `.progress=FALSE`)
    * can now create weighted networks; to do so, you must choose a weighted metric in the function argument `measure`
    * added some weighted metrics as options for `measure` (*strength*, *mod.wt*, *E.global.wt*)
    * can specify the confidence level (for calculating confidence intervals) via the `conf` argument (default: 0.95)
* `set_brainGraph_attr`:
    * New argument `xfm.type`, which allows you to choose how edge weights should be transformed for calculating distance-based metrics.
    * The default is the *reciprocal* (which is what was hard-coded in previous versions).
    * Other options are: `1-w` (subtract weights from 1); and `-log(w)` (take the negative natural logarithm of weights).

### New functions
* `symmetrize_array`: a convenience function that applies `symmetrize_mats` along the third dimension of an array
* `xfm.weights`: utility function to transform edge weights (necessary when calculating distance-based metrics).

## Minor changes
* `graph_attr_dt` and `vertex_attr_dt` will now include `weighting`, if present
* `set_brainGraph_attr` has 2 new arguments:
    1. `weighting` will create a graph-level attribute indicating how the edges are weighted (e.g., 'fa' for FA-weighted tractography networks)
    2. `threshold` will create a graph-level attribute indicating the (numeric) threshold used to create the network (if applicable)

# brainGraph 1.4.0

2017-06-10

## Bug fix
* `mtpc`: fixed a bug that would incorrectly calculate `A.crit`

## New functions
* `apply_thresholds`: threshold an additional set of matrices (e.g., FA-weighted matrices in DTI tractography) based on a set of matrices that have already been thresholded (e.g., streamline-weighted matrices in DTI tractography)

## Minor changes
* `analysis_random_graphs`: no longer requires a *covars* argument


# brainGraph 1.3.0

2017-04-30

## Bug fix
* `create_mats`
    * fixed bug for deterministic tractography when the user would like to normalize the matrices by *ROI size*.
    * Fixed bug for when `threshold.by='density'`. Previously, it would keep the top *X*% for *each* subject

## Major changes
* `create_mats`
    * `threshold.by='consensus'` is the name of the new default, as this is what is called "consensus-based" thresholding in the literature.
    * `threshold.by='consistency'` is a new option, for performing *consistency-based* thresholding. See Roberts et al., 2017.

## Minor changes
* `set_brainGraph_attr` no longer calculates the graph's *clique number*, which takes exceedingly long in denser and/or larger graphs (e.g., `craddock200`)

# brainGraph 1.2.0

2017-04-29

## Bug fix
* `plot_brainGraph`: now returns `NA` (instead of throwing an error) if the specified *subgraph* expression results in a network with 0 vertices.
* `edge_asymmetry` fixed bug when the input graph had only one contralateral connection (usually only encountered in the GUI with neighborhood plots)

## Major changes
* `create_mats`: you can specify `threshold.by='mean'`, which will threshold the matrices such that a connection will be kept if `mean(A_ij) + 2*sd(A_ij) > mat.thresh`, for each of `mat.thresh`.

## New functions
* `make_empty_brainGraph`: this is not a new function, but rather was not exported in previous versions
* `s_core`: calculate the *s-core* membership of a graph's vertices (Eidsaa & Almaas, 2013)
    * Adds a vertex attributes called `s.core` to the graph through `set_brainGraph_attr`.
    * Analogous to the *k-core* but for weighted networks.
    * The vertex attribute for *k-core* has been changed from `coreness` to `k.core` to distinguish these metrics.


# brainGraph 1.1.0

2017-04-22

## Bug fix
* `plot_brainGraph_gui` had multiple issues and a few features have been changed:
    * Overall execution should be faster than in previous versions
    * *Lobe*, *neighborhood*, and *community* selection are now in "scrolled windows" instead of drop-down lists. Multiple selections can be made either by pressing `Ctrl` and clicking, or by holding `Shift` and moving the arrow keys
    * Fixed problem with vertex colors
    * When choosing to plot *neighborhoods*, you can color the vertices based on which neighborhood they belong to (useful if multiple vertices are selected)
* `gateway_coeff` returned an error if the number of communities equals 1; this has been fixed

## New functions
* `centr_betw_comm`: calculate vertex *communicability betweenness centrality* (Estrada et al., 2009)
* `communicability`: calculate network *communicability* (Estrada & Hatano, 2008)
* `mtpc`: the *multi-threshold permutation correction (MTPC)* method for statistical inference of either vertex- or graph-level measures (Drakesmith et al., 2015)
* `symmetrize_mats`: symmetrize a connectivity matrix by either the *maximum*, *minimum*, or *average* of the off-diagonal elements. You may select one of these as an argument to `create_mats`.

## Major changes
* `brainGraph_GLM` has 2 new function arguments:
    * `level` allows you to perform inference for graph- or vertex-level measures
    * `perms` lets you specify the permutation set explicitly
* `create_mats`: All `A.norm.sub` matrices will be symmetrized, regardless of the value of `threshold.by` (previously they were only symmetrized if using `threshold.by='density'`).
    * This should not pose a problem, as the default (to take the *maximum* of the off-diagonal elements) is also the default when creating graphs in `igraph`.

## Minor changes
* `get.resid`: no longer requires a *covars* argument, as it was redundant
* `sim.rand.graph.par`: the argument *clustering* is no longer TRUE by default



# brainGraph 1.0.0

2017-04-10

*First major release; Fifth CRAN release*

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
* In the GUI, vertex order in circle plots now more closely reflect their anatomical position, being ordered by y- and x-coordinates (and within *lobe*)


# brainGraph 0.72.0

2016-10-10

*Fourth CRAN release*

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


# brainGraph 0.62.0

2016-04-22

*Third CRAN release*

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
  * Linear model specification is more limited now, though
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



# brainGraph 0.55.0

2015-12-24

*Second CRAN release*

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


# brainGraph 0.48.0

2015-12-08

*Initial CRAN release*
