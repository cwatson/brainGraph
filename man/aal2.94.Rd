\name{aal2.94}
\alias{aal2.94}
\docType{data}
\title{Coordinates for data from the AAL2 atlas}
\description{
  This is a list of spatial coordinates for the AAL2 atlas, along with
  indices for the major lobes of the brain.
}
\usage{data("aal2.94")}
\format{
  A data frame with 94 observations on the following 10 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x}}{a numeric vector of x-coordinates (internal to \code{brainGraph})}
    \item{\code{y}}{a numeric vector of y-coordinates (internal to \code{brainGraph})}
    \item{\code{z}}{a numeric vector of z-coordinates (internal to \code{brainGraph})}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Central} \code{Limbic} \code{SCGM}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
  }
}
\source{
  Rolls E.T., Joliot M., Tzourio-Mazoyer N. (2015) \emph{Implementation of a new
  parcellation of the orbitofrontal cortex in the automated anatomical labelling
  atlas}. NeuroImage, 122:1-5.
  Tzourio-Mazoyer N., Landeau B., Papathanassiou D., Crivello F., Etard O.,
  Delcroix N., Mazoyer B., Joliot M. (2002) \emph{Automated anatomical labeling
  of activations in SPM using a macroscopic anatomical parcellation of the
  MNI MRI single-subject brain}. NeuroImage, 15(1):273-289.
}
\references{
  Rolls E.T., Joliot M., Tzourio-Mazoyer N. (2015) \emph{Implementation of a new
  parcellation of the orbitofrontal cortex in the automated anatomical labelling
  atlas}. NeuroImage, 122:1-5.
  Tzourio-Mazoyer N., Landeau B., Papathanassiou D., Crivello F., Etard O.,
  Delcroix N., Mazoyer B., Joliot M. (2002) \emph{Automated anatomical labeling
  of activations in SPM using a macroscopic anatomical parcellation of the
  MNI MRI single-subject brain}. NeuroImage, 15(1):273-289.
}
\examples{
data(aal2.94)
str(aal2.94)
}
\keyword{datasets}
