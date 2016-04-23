## Test environments
* local 64-bit CentOS 6.7 install, R 3.2.3
* Ubuntu 12.04 (64-bit, on travis-ci), R 3.2.4 Revised (2016-03-16 r70338)
* win-builder (64-bit), R 3.3.0 beta (2016-04-14 r70486)
* win-builder (64-bit), R 3.2.5

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking package dependencies ... NOTE
  Package in Depends/Imports which should probably only be in LinkingTo: 'RcppEigen'

  I attempted to include the package in LinkingTo, to no avail. Downloading the
  package as-is has not resulted in any errors for me.

## Downstream dependencies
There are currently no downstream dependencies for this package.
