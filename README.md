# brainGraph
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Linux Build Status](https://travis-ci.org/cwatson/brainGraph.svg)](https://travis-ci.org/cwatson/brainGraph)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)

*brainGraph* is an R package for performing graph theory analyses of brain MRI
data. It is most useful in atlas-based analyses (e.g. using an atlas such as *AAL*
or one from *Freesurfer*); however, many of the computations (e.g., the *GLM*
functions and the network-based statistic) will still work with any graph that
is compatible with [igraph](https://github.com/igraph/rigraph). I also have used
this for tractography data from *FSL*'s *probtrackx2* and resting-state fMRI data
from *DPABI*.

# Installation
The package should work "out-of-the-box" on Linux systems, since almost all
development (and use) has been on computers running CentOS 6 and CentOS 7. I
have also had success running it (and did some development) on Windows 7, and
have heard from users that it works on some versions of Mac OS and on Ubuntu.

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
Freesurfer and FSL, and some suggestions for workflow organization. To access
the *User Guide*, please use
[this link.](https://dl.dropboxusercontent.com/s/wmupawb39bcdho3/brainGraph_UserGuide.pdf)
(NOTE: you will be asked to download the PDF)

# Graph measures
In addition to the extensive list of measures available in *igraph*, I have
functions for calculating/performing:
* Between-group differences in vertex measures (e.g., *degree*,
    *betweenness centrality*, etc.) using the General Linear Model. See Chapter
    7 of the User Guide, which was modeled after the GLM help page on FSL's wiki
* Leverage centrality (see Joyce et al., 2010)
* Asymmetry index
* Efficiency (global, nodal, and local; see Latora & Marchiori, 2001)
* "Individual contributions (*leave-one-out [LOO]* and *add-one-patient [AOP]*;
    see Saggar et al., 2015)
* The *network-based statistic (NBS)* (see Zalesky et al., 2010)
* Rich-club coefficients and normalization (see Zhou & Mondragon, 2004; and
    Colizza et al., 2006)
* The "rich-core" (see Ma & Mondragon, 2015)
* Robustness ("targeted attack" and "random failure")
* Small-worldness
* Euclidean distances of edges
* Participation coefficient and within-module degree z-score (see Guimera & Amaral, 2005)
* Gateway coefficient (see Vargas & Wahl, 2014)
* Vulnerability
* Random graph generation (also controlling for clustering; see Bansal et al., 2009)
* Bootstrapping of graph-level metrics (e.g., *modularity*)
* Permutation analysis

# Visualization
There is a plotting GUI for fast and easy data exploration that will *not* work
without data from a standard atlas (ideally to be fixed some time in the future).
You may use a custom atlas if you follow the same format as the other atlases in
the package (see Chapter 4 of the *User Guide*). A screenshot of the GUI is here:

<img src="https://www.dropbox.com/s/e0hdng7flrxsgd2/brainGraph_GUI.png?dl=1" width="843" height="372" />

# Getting Help
For bug reports, feature requests, help with usage/code/etc., please join the
*Google Group*
[brainGraph-help](https://groups.google.com/forum/?hl=en#!forum/brainGraph-help).
