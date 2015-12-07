# brainGraph
[![Build Status](https://travis-ci.org/cwatson/brainGraph.svg)](https://travis-ci.org/cwatson/brainGraph)

*brainGraph* is an R package that is used to perform graph theory analyses of
brain MRI data. It is currently most useful in ROI-based analyses (e.g. using an
atlas such as AAL or one from Freesurfer); however, the computations will still
work with any graph. There is a rudimentary plotting GUI that will *not* work
without data from a standard atlas (to be fixed some time in the future, but 
it wouldn't make too much sense for non-brain MRI data).

# Installation
If you do not have the *devtools* package, I recommend you install that first.
Then, all you need to type is
``` r
devtools::install_github('cwatson/brainGraph')
```
This should install all of the dependencies needed along with the package
itself.

# Usage
To see the *User Guide*, please use
[this link.](https://www.dropbox.com/s/j831n3q9muyz1go/brainGraph_UserGuide.pdf?dl=0)
