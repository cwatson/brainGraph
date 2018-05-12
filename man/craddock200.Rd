\name{craddock200}
\alias{craddock200}
\docType{data}
\title{Coordinates for data from the Craddock200 atlas}
\description{
  This is a list of spatial coordinates for the Craddock200 atlas, along with
  indices for the major lobes of the brain.
}
\usage{data("craddock200")}
\format{
  A data frame with 200 observations on the following 8 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Cingulate} \code{SCGM} \code{Cerebellum} \code{Brainstem}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R} \code{Brainstem}}
    \item{\code{index}}{a numeric vector}
    \item{\code{name.full}}{a character vector of longer region names}
  }
}
\references{
  Craddock, R. C., James, G. A., Holtzheimer, P. E., Hu, X. P., Mayberg, H. S.
  (2012). \emph{A whole brain fMRI atlas generated via spatially constrained
  spectral clustering}. Human Brain Mapping, 2012, 33, 1914-1928.
  doi: 10.1002/hbm.21333.
}
\examples{
data(craddock200)
str(craddock200)
}
\keyword{datasets}
