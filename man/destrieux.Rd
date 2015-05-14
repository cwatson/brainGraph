\name{destrieux}
\docType{data}
\alias{destrieux}
\title{Coordinates for Freesurfer data}
\description{
  This is a list of spatial coordinates for the Destrieux atlas, along with
  indices for the major lobes of the brain.
}
\usage{dk}
\format{List with 21 elements:
\describe{
    \item{coords}{Data frame with 148 rows and 3 columns of (x, y, z)
        coordinates}
    \item{brainnet.coords}{Data frame with 148 rows and 3 columns of (x, y, z)
        coordinates in MNI space}
    \item{lobe.names}{Character string with 6 elements naming the major lobes of
        the brain}
}
}
\references{
  Destrieux C., Fischl B., Dale E. \& Halgren E. (2010) \emph{Automatic
  parcellation of human cortical gyri and sulci using standard anatomic
  nomenclature}. NeuroImage, 53(1):1-15.
}
\keyword{datasets}
