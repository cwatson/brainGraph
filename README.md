# brainGraph
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Linux Build Status](https://travis-ci.org/cwatson/brainGraph.svg)](https://travis-ci.org/cwatson/brainGraph)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-ago/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)
[![Dependencies](https://tinyverse.netlify.com/badge/brainGraph)](https://cran.r-project.org/package=brainGraph)
[![GPL License](https://img.shields.io/cran/l/brainGraph.svg)](https://opensource.org/licenses/GPL-3.0/)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)
[![HitCount](http://hits.dwyl.io/cwatson/brainGraph.svg)](http://hits.dwyl.io/cwatson/brainGraph)

`brainGraph` ([RRID: SCR_017260](https://scicrunch.org/resolver/SCR_017260)) is an
[R package](https://cran.r-project.org/web/packages/brainGraph/index.html) for performing graph theory
analyses of brain MRI data. It is most useful in atlas-based analyses (e.g., using an atlas such as
[AAL](http://www.gin.cnrs.fr/en/tools/aal-aal2/),
or one from [Freesurfer](https://surfer.nmr.mgh.harvard.edu/)); however, many of
the computations (e.g., the *GLM*-based
functions and the network-based statistic) will work with any graph that
is compatible with [igraph](https://github.com/igraph/rigraph). The package will
perform analyses for *structural covariance networks (SCN)*, DTI tractography
(I use *probtrackx2* from [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)), and
resting-state fMRI covariance (I have used the Matlab-based [DPABI](http://rfmri.org/dpabi)
toolbox).

Table of Contents
====
<!-- vim-markdown-toc GFM -->

* [Requirements](#requirements)
    * [Operating Systems](#operating-systems)
    * [Multi-core processing](#multi-core-processing)
* [Compatibility](#compatibility)
    * [Neuroimaging software](#neuroimaging-software)
    * [Brain atlases](#brain-atlases)
        * [Other atlases](#other-atlases)
* [Installation](#installation)
    * [Multi-core processing](#multi-core-processing-1)
    * [GUI](#gui)
* [Usage - the User Guide](#usage---the-user-guide)
* [Graph measures](#graph-measures)
    * [Group analyses](#group-analyses)
        * [GLM-based](#glm-based)
        * [Non-GLM based](#non-glm-based)
    * [Null graph-related measures](#null-graph-related-measures)
    * [Other measures](#other-measures)
* [Visualization](#visualization)
* [Getting Help](#getting-help)
* [Future versions](#future-versions)

<!-- vim-markdown-toc -->
# Requirements

## Operating Systems
The package ***should*** work "out-of-the-box" on *Linux* systems (at least on Red
Hat-based systems; i.e., CentOS, RHEL, Scientific Linux, etc.) since almost all
development (and use, by me) has been on computers running CentOS 6 and (currently)
CentOS 7. I have also had success running it (and did some development) on
Windows 7, and have heard from users that it works on some versions of Mac OS
and on Ubuntu. Please see the User Guide (mentioned below) for more details.

## Multi-core processing
Many `brainGraph` functions utilize multiple CPU cores. This is primarily done
via the [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
package. Depending on your OS, you will need to install
[doMC](https://cran.r-project.org/web/packages/doMC/index.html) (*macOS* and *Linux*)
or [doSNOW](https://cran.r-project.org/web/packages/doSNOW/index.html)
(*Windows*).

# Compatibility
## Neuroimaging software
I mostly use Freesurfer and FSL, but the following software packages should be suitable.
Note that this is an incomplete list; any software that can output a connectivity matrix will work.
* [Freesurfer](https://surfer.nmr.mgh.harvard.edu)
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
* [DPARSF](http://rfmri.org/DPARSF)
* [PANDA](https://www.nitrc.org/projects/panda)
* [TrackVis](http://trackvis.org/)

## Brain atlases
There are several brain atlases for which the data are present in `brainGraph`.
Atlases containing `.scgm` in the name contain both cortical and *SubCortical Gray Matter (SCGM)* regions.
1. `dk` and `dk.scgm`: [Desikan-Killiany](https://dx.doi.org/10.1016/j.neuroimage.2006.01.021)
2. `dkt` and `dkt.scgm`: [Desikan-Killiany-Tourville](https://dx.doi.org/10.3389/fnins.2012.00171)
3. `destrieux` and `destrieux.scgm`: [Destrieux](https://dx.doi.org/10.1016/j.neuroimage.2010.06.010)
4. `aal90` and `aal116`: [Automated Anatomical Labeling atlas](https://dx.doi.org/10.1006/nimg.2001.0978)
5. `aal2.94` and `aal2.120`: [AAL-2](https://dx.doi.org/10.1016/j.neuroimage.2015.07.075)
6. `brainsuite`: [Brainsuite](https://dx.doi.org/10.1016/S1361-8415(02)00054-3)
7. `craddock200`: [Craddock-200](https://dx.doi.org/10.1002/hbm.21333)
8. `dosenbach160`: [Dosenbach-160](https://dx.doi.org/10.1126/science.1194144)
9. `hoa112`: [Harvard-Oxford atlas](https://dx.doi.org/10.1016/j.schres.2005.11.020)
10. `lpba40`: [LONI Probabilistic Brain Atlas](https://dx.doi.org/10.1016/j.neuroimage.2007.09.031)

### Other atlases
Some functions accept a `custom.atlas` argument, so that you can analyze data that is from an atlas not present in `brainGraph`.
Other atlases to be added in the future include the following (I would need specific coordinate, region name, and lobe and hemisphere information):
* [HCP-1mm](https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446)
    (see [Glasser et al., 2016](https://doi.org/10.1038/nature18933))
* [Power-264](https://dx.doi.org/10.1016%2Fj.neuron.2011.09.006)
* [Gordon-333](https://doi.org/10.1093/cercor/bhu239)
* [Shen-268](https://doi.org/10.1016/j.neuroimage.2013.05.081)
* [Von Economo-Koskinas](https://doi.org/10.1016/j.neuroimage.2016.12.069)
* [Brainnetome](https://doi.org/10.1093/cercor/bhw157)
* [Willard-499](https://findlab.stanford.edu/functional_ROIs.html) (see [Richiardi et al., 2015](https://doi.org/10.1126/science.1255905))
* [Schaefer-400](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal)
    (see [Schaefer et al., 2018](https://doi.org/10.1093/cercor/bhx179))

# Installation
There are (primarily) two ways to install this package:

1. Directly from CRAN: (use one of the following commands)
``` r
install.packages('brainGraph')
install.packages('brainGraph', dependencies=TRUE)
```

2. From the GitHub repo (for development versions). This requires that the
[devtools](https://cran.r-project.org/web/packages/devtools/index.html)
package be installed:
``` r
devtools::install_github('cwatson/brainGraph')
```
This should install all of the dependencies needed along with the package
itself. For more details, see the [User Guide (PDF link)](https://cwatson.github.io/files/brainGraph_UserGuide.pdf).

## Multi-core processing
To set up your R session for parallel processing, you can use the following code. Note that it is different for *Windows*.
This code should be run before any data processing. If you will always use a single OS, you can remove the unnecessary lines.
``` r
OS <- .Platform$OS.type
if (OS == 'windows') {
  library(snow)
  library(doSNOW)
  num.cores <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
  cl <- makeCluster(num.cores, type='SOCK')
  clusterExport(cl, 'sim.rand.graph.par')   # Or whatever functions you will use
  registerDoSNOW(cl)
} else {
  library(doMC)
  registerDoMC(detectCores() - 1L)  # Keep 1 core free
}
```

For example, I source the following simple script before I do any parallel processing with `brainGraph`:
``` r
pacman::p_load(brainGraph, doMC)
registerDoMC(detectCores())
```

## GUI
On some systems (e.g., *macOS* and *Windows*) it might be very difficult to
install the necessary packages/dependencies for the GUI functions. Since `v2.2.0` (released 2018-05-28),
the R packages [RGtk2](https://cran.r-project.org/web/packages/RGtk2/index.html)
and [cairoDevice](https://cran.r-project.org/web/packages/cairoDevice/index.html)
have been changed to *Suggests* (i.e., they are no longer required),
so it can be installed on a "headless" server.

If you are on *macOS* or *Windows* and would like GUI functionality, please see
[this GitHub Gist](https://gist.github.com/sebkopf/9405675#macos). The comments
contain more recent information.
You may also need to install a few additional packages, shown here:
``` r
install.packages('gWidgets', dependencies=TRUE)
install.packages('gWidgetsRGtk2', dependencies=TRUE)
install.packages('RGtk2Extras', dependencies=TRUE)
```

# Usage - the User Guide
I have a User Guide that contains *extensive* code examples for analyses common to
brain MRI studies. I also include some code for getting your data *into* R *from*
Freesurfer, FSL, and DPABI, and some suggestions for workflow organization.

The User Guide is the most complete documentation of this package.
If you are a beginner using R, I encourage you to read it thoroughly.
You may start with the [Preface](https://cwatson.github.io/files/brainGraph_UserGuide.pdf#chapter*.1)
or at whichever chapter is suitable for your analyses.

To access the User Guide, a PDF is available at
[this link.](https://cwatson.github.io/files/brainGraph_UserGuide.pdf)

# Graph measures
In addition to the extensive list of measures available in *igraph*, I have
functions for calculating/performing:

## Group analyses
There are several analyses based on the *General Linear Model (GLM)*, and others
that have different purposes.

### GLM-based
* Between-group differences in vertex- or graph-level measures (e.g., *degree*,
    *betweenness centrality*, *global efficiency*, etc.) using the GLM's.
    See [Chapter 8](https://cwatson.github.io/files/brainGraph_UserGuide.pdf#chapter.8) of the User Guide, which was partly modeled after the
    GLM page on the [FSL wiki](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM)
* The *multi-threshold permutation correction (MTPC)* method for statistical
    inference (see [Drakesmith et al., 2015](https://dx.doi.org/10.1016%2Fj.neuroimage.2015.05.011)
    and [Chapter 9](https://cwatson.github.io/files/brainGraph_UserGuide.pdf#chapter.9) of the User Guide)
* The *network-based statistic (NBS)* (see [Zalesky et al., 2010](https://doi.org/10.1016/j.neuroimage.2010.06.041)
    and [Chapter 10](https://cwatson.github.io/files/brainGraph_UserGuide.pdf#chapter.10) of the User Guide)
* Graph- and vertex-level *mediation* analysis (see [Chapter 11](https://cwatson.github.io/files/brainGraph_UserGuide.pdf#chapter.11)
    of the User Guide, and the [mediation](https://cran.r-project.org/web/packages/mediation/index.html) package in R)

### Non-GLM based
* Bootstrapping of graph-level metrics (e.g., *modularity*)
* Permutation analysis of between-group differences in vertex- or graph-level measures
* "Individual contributions" (*leave-one-out [LOO]* and *add-one-patient [AOP]*;
    see [Saggar et al., 2015](https://dx.doi.org/10.1016%2Fj.neuroimage.2015.07.006))

## Null graph-related measures
* Null/random graph generation (both the "standard" method, and also a method controlling for clustering;
    see [Bansal et al., 2009](https://doi.org/10.1186/1471-2105-10-405))
* Small-worldness: the "original" of [Watts & Strogatz, 1998](https://dx.doi.org/10.1515/9781400841356.301)
    and [Humphries et al., 2008](https://doi.org/10.1371/journal.pone.0002051),
    and "omega" introduced in [Telesford et al., 2011](https://doi.org/10.1089/brain.2011.0038)
* Rich-club coefficients and normalization (see [Zhou & Mondragon, 2004](https://doi.org/10.1109/LCOMM.2004.823426); and
    [Colizza et al., 2006](https://doi.org/10.1038/nphys209))

## Other measures
* Efficiency (global, nodal, and local; see [Latora & Marchiori, 2001](https://dx.doi.org/10.1103/PhysRevLett.87.198701))
* The *rich-core* (see [Ma & Mondragon, 2015](https://dx.doi.org/10.1371/journal.pone.0119678))
* Leverage centrality (see [Joyce et al., 2010](https://dx.doi.org/10.1371/journal.pone.0012200))
* Asymmetry index
* Robustness ("targeted attack" and "random failure") and vulnerability
* Euclidean distances of edges
* Participation coefficient and within-module degree z-score (see [Guimera & Amaral, 2005a](https://dx.doi.org/10.1038/nature03288) and [2005b](https://dx.doi.org/10.1088/1742-5468/2005/02/P02001))
* Gateway coefficient (see [Vargas & Wahl, 2014](https://dx.doi.org/10.1140/epjb/e2014-40800-7))
* *Communicability* and *communicability betweenness* (see [Estrada & Hatano, 2008](https://dx.doi.org/10.1103/PhysRevE.77.036111);
    [Estrada et al., 2009](https://dx.doi.org/10.1016/j.physa.2008.11.011);
    [Crofts & Higham, 2009](https://dx.doi.org/10.1098/rsif.2008.0484))
* Vertex *s-core* membership (see [Eidsaa & Almaas, 2013](https://dx.doi.org/10.1103/PhysRevE.88.062819))

# Visualization
There is a plotting GUI for fast and easy data exploration that will *not* work
without data from a standard atlas (ideally to be extended some time in the future).
You may use a custom atlas if you follow the same format as the other atlases in
the package (see [Chapter 4](https://cwatson.github.io/files/brainGraph_UserGuide.pdf#chapter.4)
of the *User Guide* for instructions).

![brainGraph GUI](http://i.imgur.com/KKZZqI2.png)

# Getting Help
For bug reports, feature requests, help with usage/code/etc., please join the
*Google Group*
[brainGraph-help](https://groups.google.com/forum/?hl=en#!forum/brainGraph-help).
You may also consult the User Guide, and you can
[open an issue](https://github.com/cwatson/brainGraph/issues) here on GitHub.

# Future versions
An *incomplete* list of features/functionality I plan on adding to future versions:
- Longitudinal modeling (with *linear mixed effects (LME)* models)
- Thresholding and graph creation using the *minimum spanning tree* as a base
- Thresholding and graph creation for resting-state fMRI using a technique such as the *graphical lasso*
- Write functions to print group analysis results in [xtable](https://cran.r-project.org/web/packages/xtable/index.html) format for `LaTeX` documents
