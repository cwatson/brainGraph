\name{dkt.scgm}
\alias{dkt.scgm}
\docType{data}
\title{Coordinates for data from the Desikan-Killiany-Tourville atlas}
\description{
  This is a list of spatial coordinates for the Desikan-Killiany-Tourville (DKT)
  atlas, along with indices for the major lobes of the brain and subcortical
  gray matter structures.

}
\usage{data("dkt.scgm")}
\format{
  A data frame with 76 observations on the following 10 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x}}{a numeric vector of x-coordinates (internal to \code{brainGraph})}
    \item{\code{y}}{a numeric vector of y-coordinates (internal to \code{brainGraph})}
    \item{\code{z}}{a numeric vector of z-coordinates (internal to \code{brainGraph})}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Cingulate}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
  }
}
\source{
  Klein A. and Tourville J. (2012) \emph{101 labeled brain images and a
  consistent human cortical labeling protocol}. Front Neurosci,
  doi:10.3389/fnins.2012.00171
}
\references{
  Klein A. and Tourville J. (2012) \emph{101 labeled brain images and a
  consistent human cortical labeling protocol}. Front Neurosci,
  doi:10.3389/fnins.2012.00171
}
\examples{
data(dkt.scgm)
str(dkt.scgm)
}
\keyword{datasets}
