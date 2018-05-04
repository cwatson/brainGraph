\name{dkt}
\alias{dkt}
\docType{data}
\title{Coordinates for data from the Desikan-Killiany-Tourville atlas}
\description{
  This is a list of spatial coordinates for the Desikan-Killiany-Tourville (DKT)
  atlas, along with indices for the major lobes of the brain.

}
\usage{data("dkt")}
\format{
  A data frame with 62 observations on the following 8 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Cingulate}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
    \item{\code{name.full}}{a character vector of full region names}
  }
}
\references{
  Klein A. and Tourville J. (2012) \emph{101 labeled brain images and a
  consistent human cortical labeling protocol}. Front Neurosci,
  doi:10.3389/fnins.2012.00171
}
\examples{
data(dkt)
str(dkt)
}
\keyword{datasets}
