## Test environments
* local 64-bit CentOS 6.7 install, R 3.2.2
* Ubuntu 12.04 (on travis-ci), R 3.2.2
* local 32-bit Windows 7 install, R 3.1.0
* local 32-bit Windows 7 install, R-devel 3.2.0 (r67580)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking R code for possible problems ... NOTE
  assign_lobes: no visible global function definition for '.'
  plot_group_means: no visible global function definition for '.'

  This notation is an alias to "list" in "data.table" code.

## Downstream dependencies
There are currently no downstream dependencies for this package.
