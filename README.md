# brainGraph
[![Build Status](https://travis-ci.org/cwatson/brainGraph.svg)](https://travis-ci.org/cwatson/brainGraph)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)
[![](http://cranlogs.r-pkg.org/badges/brainGraph)](http://cran.rstudio.com/web/packages/brainGraph/index.html)

*brainGraph* is an R package that is used to perform graph theory analyses of
brain MRI data. It is currently most useful in atlas-based analyses (e.g. using
an atlas such as AAL or one from Freesurfer); however, the computations will
still work with any graph. There is a plotting GUI for fast and easy data
exploration that will *not* work without data from a standard atlas (ideally to
be fixed some time in the future).

# Installation
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
itself.

# Usage
To see the *User Guide*, please use
[this link.](https://www.dropbox.com/s/j831n3q9muyz1go/brainGraph_UserGuide.pdf?dl=0)
