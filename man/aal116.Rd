\name{aal116}
\alias{aal116}
\docType{data}
\title{Coordinates for data from the AAl116 atlas}
\description{
  This is a list of spatial coordinates for the AAL116 atlas, along with
  indices for the major lobes of the brain.
}
\usage{data("aal116")}
\format{
  A data frame with 116 observations on the following 7 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Limbic} \code{SCGM} \code{Cerebellum}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
  }
}
\references{
  Tzourio-Mazoyer N., Landeau B., Papathanassiou D., Crivello F., Etard O.,
  Delcroix N., Mazoyer B., Joliot M. (2002) \emph{Automated anatomical labeling
  of activations in SPM using a macroscopic anatomical parcellation of the
  MNI MRI single-subject brain}. NeuroImage, 15(1):273-289.
}
\examples{
data(aal116)
str(aal116)
}
\keyword{datasets}
