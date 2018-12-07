# brainGraph
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Linux Build Status](https://travis-ci.org/cwatson/brainGraph.svg)](https://travis-ci.org/cwatson/brainGraph)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-ago/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)
[![GPL License](https://img.shields.io/cran/l/brainGraph.svg)](https://opensource.org/licenses/GPL-3.0/)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)

`brainGraph` is an R package for performing graph theory analyses of brain MRI
data. It is most useful in atlas-based analyses (e.g., using an atlas such as
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

* [Installation](#installation)
    * [Installation commands](#installation-commands)
    * [Multi-core processing](#multi-core-processing)
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
# Installation
The package ***should*** work "out-of-the-box" on *Linux* systems (at least on Red
Hat-based systems; i.e., CentOS, RHEL, Scientific Linux, etc.) since almost all
development (and use, by me) has been on computers running CentOS 6 and (currently)
CentOS 7. I have also had success running it (and did some development) on
Windows 7, and have heard from users that it works on some versions of Mac OS
and on Ubuntu. Please see the User Guide (mentioned below) for more details.

## Installation commands
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
itself. For more details, see the User Guide (link to PDF in next section).

## Multi-core processing
Many `brainGraph` functions utilize multiple CPU cores. This is primarily done
via the [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
package. Depending on your OS, you may need to install
[doMC](https://cran.r-project.org/web/packages/doMC/index.html) (*macOS* and *Linux*)
or [doSNOW](https://cran.r-project.org/web/packages/doSNOW/index.html)
(*Windows*).

Then, to set up your R session for parallel processing:
``` r
OS <- .Platform$OS.type
if (OS == 'windows') {
  library(snow)
  library(doSNOW)
  num.cores <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
  cl <- makeCluster(num.cores, type='SOCK')
  clusterExport(cl, 'sim.rand.graph.par')   # Or whatever function you will use
  registerDoSNOW(cl)
} else {
  library(doMC)
  num.cores <- detectCores()
  registerDoMC(num.cores)
```

## GUI
On some systems (e.g., *macOS* and *Windows*) it might be very difficult to
install the necessary packages/dependencies for the GUI functions. Since `v2.2.0`,
the R packages [RGtk2](https://cran.r-project.org/web/packages/RGtk2/index.html)
and [cairoDevice](https://cran.r-project.org/web/packages/cairoDevice/index.html)
have been changed to *Suggests* (i.e., they are no longer required), so
installation should still be possible. This also means it should be possible to
install on a headless server.

# Usage - the User Guide
I have a User Guide that contains *extensive* code examples for analyses common to
brain MRI studies. I also include some code for getting your data *into* R *from*
Freesurfer, FSL, and DPABI, and some suggestions for workflow organization. The
User Guide is the most complete documentation of this package, so I encourage you
to read it thoroughly. Please start with the *Preface*.

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
    See Chapter 8 of the User Guide, which was partly modeled after the
    GLM page on the [FSL wiki](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM)
* The *multi-threshold permutation correction (MTPC)* method for statistical
    inference (see Drakesmith et al., 2015 and Chapter 9 of the User Guide)
* The *network-based statistic (NBS)* (see Zalesky et al., 2010 and Chapter 10 of the User Guide)
* Graph- and vertex-level *mediation* analysis (see Chapter 11 of the User Guide, and the [mediation](https://cran.r-project.org/web/packages/mediation/index.html) package in R)

### Non-GLM based
* Bootstrapping of graph-level metrics (e.g., *modularity*)
* Permutation analysis of between-group differences in vertex- or graph-level measures
* "Individual contributions" (*leave-one-out [LOO]* and *add-one-patient [AOP]*;
    see Saggar et al., 2015)

## Null graph-related measures
* Null/random graph generation (both the "standard" method, and also a method controlling for clustering; see Bansal et al., 2009)
* Small-worldness (the "original" of Watts & Strogatz, 1998 and Humphries et al., 2008; and "omega" introduced in Telesford et al., 2011)
* Rich-club coefficients and normalization (see Zhou & Mondragon, 2004; and
    Colizza et al., 2006)

## Other measures
* Efficiency (global, nodal, and local; see Latora & Marchiori, 2001)
* The *rich-core* (see Ma & Mondragon, 2015)
* Leverage centrality (see Joyce et al., 2010)
* Asymmetry index
* Robustness ("targeted attack" and "random failure") and vulnerability
* Euclidean distances of edges
* Participation coefficient and within-module degree z-score (see Guimera & Amaral, 2005)
* Gateway coefficient (see Vargas & Wahl, 2014)
* *Communicability* and *communicability betweenness* (see Estrada & Hatano, 2008; Estrada et al., 2009; Crofts & Higham, 2009)
* Vertex *s-core* membership (see Eidsaa & Almaas, 2013)

# Visualization
There is a plotting GUI for fast and easy data exploration that will *not* work
without data from a standard atlas (ideally to be fixed some time in the future).
You may use a custom atlas if you follow the same format as the other atlases in
the package (see Chapter 4 of the *User Guide* for instructions).

![brainGraph GUI](http://i.imgur.com/KKZZqI2.png)

# Getting Help
For bug reports, feature requests, help with usage/code/etc., please join the
*Google Group*
[brainGraph-help](https://groups.google.com/forum/?hl=en#!forum/brainGraph-help).
You may also consult the User Guide, and you can
[open an issue](https://github.com/cwatson/brainGraph/issues) here on GitHub.

# Future versions
An *incomplete* list of features/functionality I plan on adding to future versions:
- Easy import of *Freesurfer* data for structural covariance networks
- Longitudinal modeling (with *linear mixed effects (LME)* models)
- Thresholding and graph creation using the *minimum spanning tree* as a base
- Thresholding and graph creation for resting-state fMRI using a technique such as the *graphical lasso*
- Write functions to print group analysis results in [xtable](https://cran.r-project.org/web/packages/xtable/index.html) format for `LaTeX` documents
- Add data tables for multiple atlases/parcellations (e.g., *Gordon* [see Gordon et al., 2016, Cerebral Cortex], *HCP*, etc.)
