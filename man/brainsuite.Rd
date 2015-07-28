\name{brainsuite}
\alias{brainsuite}
\docType{data}
\title{Coordinates for data from BrainSuite atlas}
\description{
  This is a list of spatial coordinates for the BrainSuite software, along with
  indices for the major lobes of the brain.
}
\usage{data("brainsuite")}
\format{
  A data frame with 74 observations on the following 10 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x}}{a numeric vector of x-coordinates (internal to \code{brainGraph})}
    \item{\code{y}}{a numeric vector of y-coordinates (internal to \code{brainGraph})}
    \item{\code{z}}{a numeric vector of z-coordinates (internal to \code{brainGraph})}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Cingulate} \code{SCGM}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
  }
}
\source{
  Shattuck DW and Leahy RM (2002) \emph{BrainSuite: an automated cortical
  surface identification tool}. Medical Image Analysis, 8(2):129-142.
}
\references{
  Pantazis D, Joshi AA, Jintao J, Shattuck DW, Bernstein LE, Damasio H, and
  Leahy RM. (2009) \emph{Comparison of landmark-based and automatic methods
      for cortical surface registration}. NeuroImage, 49(3):2479-2493.
}
\examples{
data(brainsuite)
str(brainsuite)
}
\keyword{datasets}
