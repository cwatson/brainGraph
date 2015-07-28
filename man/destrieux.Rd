\name{destrieux}
\alias{destrieux}
\docType{data}
\title{Coordinates for data from the Destrieux atlas}
\description{
  This is a list of spatial coordinates for the Destrieux atlas, along with
  indices for the major lobes of the brain.
}
\usage{data("destrieux")}
\format{
  A data frame with 148 observations on the following 11 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x}}{a numeric vector of x-coordinates (internal to \code{brainGraph})}
    \item{\code{y}}{a numeric vector of y-coordinates (internal to \code{brainGraph})}
    \item{\code{z}}{a numeric vector of z-coordinates (internal to \code{brainGraph})}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Limbic}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
    \item{\code{class}}{a factor with levels \code{G} \code{G_and_S} \code{S}}
  }
}
\source{
  Destrieux C., Fischl B., Dale E. \& Halgren E. (2010) \emph{Automatic
  parcellation of human cortical gyri and sulci using standard anatomic
  nomenclature}. NeuroImage, 53(1):1-15.
}
\references{
  Destrieux C., Fischl B., Dale E. \& Halgren E. (2010) \emph{Automatic
  parcellation of human cortical gyri and sulci using standard anatomic
  nomenclature}. NeuroImage, 53(1):1-15.
}
\examples{
data(destrieux)
str(destrieux)
}
\keyword{datasets}
