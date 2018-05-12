\name{lpba40}
\alias{lpba40}
\docType{data}
\title{Coordinates for data from the LONI probabilistic brain atlas}
\description{
  This is a list of spatial coordinates for the LPBA40 atlas, along with
  indices for the major lobes of the brain. The coordinates were obtained
  from some colleagues.
}
\usage{data("lpba40")}
\format{
  A data frame with 56 observations on the following 7 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Cingulate} \code{SCGM}}
    \item{\code{hemi}}{a factor with levels \code{} \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
  }
}
\references{
  Shattuck DW, Mirza M, Adisetiyo V, Hojatkashani C, Salamon G, Narr KL,
  Poldrack RA, Bilder RM, Toga AW (2007) \emph{Construction of a 3D probabilistic
  atlas of human cortical structures}. NeuroImage,
  doi:10.1016/j.neuroimage.2007.09.031
}
\examples{
data(lpba40)
str(lpba40)
}
\keyword{datasets}
