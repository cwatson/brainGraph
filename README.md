# brainGraph
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Linux Build Status](https://travis-ci.org/cwatson/brainGraph.svg)](https://travis-ci.org/cwatson/brainGraph)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-ago/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)
[![GPL License](https://img.shields.io/cran/l/brainGraph.svg)](https://opensource.org/licenses/GPL-3.0/)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)

*brainGraph* is an R package for performing graph theory analyses of brain MRI
data. It is most useful in atlas-based analyses (e.g. using an atlas such as *AAL*
or one from *Freesurfer*); however, many of the computations (e.g., the *GLM*
functions and the network-based statistic) will still work with any graph that
is compatible with [igraph](https://github.com/igraph/rigraph). I also have used
this for tractography data from *FSL*'s *probtrackx2* and resting-state fMRI data
from *DPABI*.

Table of Contents
====
<!-- vim-markdown-toc GFM -->

* [Installation](#installation)
* [Usage](#usage)
* [Graph measures](#graph-measures)
    * [Group analyses](#group-analyses)
    * [Null graph-related measures](#null-graph-related-measures)
    * [Other measures](#other-measures)
* [Visualization](#visualization)
* [Getting Help](#getting-help)
* [Future versions](#future-versions)

<!-- vim-markdown-toc -->
# Installation
The package ***should*** work "out-of-the-box" on Linux systems (at least on Red
Hat-based systems; i.e., CentOS, RHEL, Scientific Linux, etc.) since almost all
development (and use, by me) has been on computers running CentOS 6 and (currently)
CentOS 7. I have also had success running it (and did some development) on
Windows 7, and have heard from users that it works on some versions of Mac OS
and on Ubuntu.

There are two ways to install this package:

* Directly from CRAN:
``` r
install.packages('brainGraph')
```

* For development versions, using the *devtools* package (install first, if
necessary):
``` r
devtools::install_github('cwatson/brainGraph')
```
This should install all of the dependencies needed along with the package
itself. For more details, see the User Guide (link to PDF in next section).

# Usage
I have a User Guide that contains extensive code examples for analyses common to
brain MRI studies. I also include some code for getting your data *into* R *from*
Freesurfer, FSL, and DPABI, and some suggestions for workflow organization. To
access the *User Guide*, please use
[this link.](https://dl.dropboxusercontent.com/s/wmupawb39bcdho3/brainGraph_UserGuide.pdf)
(NOTE: you will be asked to download the PDF)

# Graph measures
In addition to the extensive list of measures available in *igraph*, I have
functions for calculating/performing:

## Group analyses
* Between-group differences in vertex- or graph-level measures (e.g., *degree*,
    *betweenness centrality*, *global efficiency*, etc.) using the General
    Linear Model. See Chapter 7 of the User Guide, which was modeled after the
    GLM help page on FSL's wiki
* The *multi-threshold permutation correction (MTPC)* method for statistical
    inference (see Drakesmith et al., 2015)
* The *network-based statistic (NBS)* (see Zalesky et al., 2010)
* Bootstrapping of graph-level metrics (e.g., *modularity*)
* Permutation analysis of between-group differences in vertex- or graph-level measures
* "Individual contributions (*leave-one-out [LOO]* and *add-one-patient [AOP]*;
    see Saggar et al., 2015)

## Null graph-related measures
* Null/random graph generation (both the "standard" method, and also a method controlling for clustering; see Bansal et al., 2009)
* Small-worldness (the "original" of Watts & Strogatz, 1998 and Humphries et al., 2008; and "omega" introduced in Telesford et al., 2011)
* Rich-club coefficients and normalization (see Zhou & Mondragon, 2004; and
    Colizza et al., 2006)

## Other measures
* Efficiency (global, nodal, and local; see Latora & Marchiori, 2001)
* The "rich-core" (see Ma & Mondragon, 2015)
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
the package (see Chapter 4 of the *User Guide*). A screenshot of the GUI is here:

![brainGraph GUI](http://i.imgur.com/KKZZqI2.png)

# Getting Help
For bug reports, feature requests, help with usage/code/etc., please join the
*Google Group*
[brainGraph-help](https://groups.google.com/forum/?hl=en#!forum/brainGraph-help).

# Future versions
An incomplete list of features/functionality I plan on adding to future versions:
* Mediation analysis at both the graph- and vertex-level (currently only simple linear models)
* Thresholding and graph creation using the *minimum spanning tree* as a base
* Create *methods* for objects of class `brainGraph`, making it simpler to, for example, show a text summary of a graph or to plot the graph over the MNI slice, etc.
    * `plot_brainGraph_mni` will be removed as it will be a redundant step
* Use *methods* for GLM-related functions/objects. For example, I am writing a function for `GLM diagnostics` which will mimic the functionality of the base R function `plot.lm`
