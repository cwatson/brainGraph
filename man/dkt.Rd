\name{dkt}
\docType{data}
\alias{dkt}
\title{Coordinates for Freesurfer data}
\description{
  This is a list of spatial coordinates for the DKT atlas, along with
  indices for the major lobes of the brain.
}
\usage{dkt}
\format{List with 21 elements:
\describe{
    \item{coords}{Data frame with 62 rows and 3 columns of (x, y, z) coordinates}
    \item{brainnet.coords}{Data frame with 62 rows and 3 columns of (x, y, z)
        coordinates in MNI space}
    \item{lobe.names}{Character string with 6 elements naming the major lobes of
        the brain}
}
}
\references{
  Klein A. and Tourville J. (2012) \emph{101 labeled brain images and a
  consistent human cortical labeling protocol}. Front Neurosci,
  doi:10.3389/fnins.2012.00171
}
\keyword{datasets}
