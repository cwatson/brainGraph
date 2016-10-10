\name{hoa112}
\alias{hoa112}
\docType{data}
\title{Coordinates for data from Harvard-Oxford atlas}
\description{
  This is a list of spatial coordinates for the Harvard-Oxford atlas, along with
  indices for the major lobes of the brain.
}
\usage{data("hoa112")}
\format{
  A data frame with 112 observations on the following 7 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Cingulate} \code{SCGM}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R}}
    \item{\code{index}}{a numeric vector}
  }
}
\source{
    Makris N., Goldstein J.M., Kennedy D. et al. (2006) \emph{Decreased volume
    of left and total anterior insular lobule in schizophrenia}. Schizophr Res,
    83(2-3):155-171.
}
\references{
    Makris N., Goldstein J.M., Kennedy D. et al. (2006) \emph{Decreased volume
    of left and total anterior insular lobule in schizophrenia}. Schizophr Res,
    83(2-3):155-171.
}
\examples{
data(hoa112)
str(hoa112)
}
\keyword{datasets}
