## Test environments
* local 64-bit CentOS 7.8.2003 install, R version 3.6.0 (2019-04-26 r76424)
* win-builder (64-bit), R version 4.0.2 (2020-06-22)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Christopher G. Watson <cgwatson@bu.edu>'

Found the following (possibly) invalid URLs:
  URL: https://dx.doi.org/10.1088/1742-5468/2005/02/P02001
    From: man/vertex_roles.Rd
    Status: Error
    Message: libcurl error code 35:
        schannel: next InitializeSecurityContext failed: SEC_E_ILLEGAL_MESSAGE (0x80090326) - This error usually occurs when a fatal SSL/TLS alert is received (e.g. handshake failed).

I have checked the URL and it works and leads to the correct reference.

## Downstream dependencies
There are currently no downstream dependencies for this package.
