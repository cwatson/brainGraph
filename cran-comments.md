## Test environments
* local 64-bit CentOS 7.3 install, R version 3.3.3 (2017-03-06)
* win-builder (64-bit), R version 3.4.0 beta (2017-04-07 r72496)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE

    Maintainer: 'Christopher G. Watson <cgwatson@bu.edu>'

    Possibly mis-spelled words in DESCRIPTION:
        DPABI (11:62)
        FSL (11:39)
        Freesurfer (9:44)
        fMRI (11:47)
        gyrification (10:31)
        tractography (11:16)

The first 3 listed are MRI analysis software packages. "fMRI" is the acronym for
"functional MRI". The last two are neuroscience terms.

## Downstream dependencies
There are currently no downstream dependencies for this package.
