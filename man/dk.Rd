\name{dk}
\alias{dk}
\docType{data}
\title{Coordinates for data from the Desikan-Killiany atlas}
\description{
  This is a list of spatial coordinates for the Desikan-Killiany (DK) atlas,
  along with indices for the major lobes of the brain.

}
\usage{data("dk")}
\format{
  A data frame with 68 observations on the following 11 variables.
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
    \item{\code{name.full}}{a character vector of full region names}
  }
}
\source{
  Desikan R.S., Segonne F., Fischl B., et al. (2006) \emph{An automated
  labeling system for subdividing the human cerebral cortex on MRI scans into
  gyral based regions of interest}. NeuroImage, 31:968-980.
}
\references{
  Desikan R.S., Segonne F., Fischl B., et al. (2006) \emph{An automated
  labeling system for subdividing the human cerebral cortex on MRI scans into
  gyral based regions of interest}. NeuroImage, 31:968-980.
}
\examples{
data(dk)
str(dk)
}
\keyword{datasets}
