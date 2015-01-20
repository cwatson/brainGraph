\name{dk}
\docType{data}
\alias{dk}
\title{Coordinates for Freesurfer data}
\description{
  This is a list of spatial coordinates for the DK atlas, along with
  indices for the major lobes of the brain.
}
\usage{dk}
\format{List with 21 elements:
\describe{
    \item{coords}{Data frame with 68 rows and 3 columns of (x, y, z)
        coordinates}
    \item{brainnet.coords}{Data frame with 68 rows and 3 columns of (x, y, z)
        coordinates in MNI space}
    \item{lobe.names}{Character string with 6 elements naming the major lobes of
        the brain}
}
}
\references{
  Desikan R.S., Segonne F., Fischl B., et al. (2006) \emph{An automated
  labeling system for subdividing the human cerebral cortex on MRI scans into
  gyral based regions of interest}. NeuroImage, 31:968-980.
}
\keyword{datasets}
