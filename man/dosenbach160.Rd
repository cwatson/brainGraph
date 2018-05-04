\name{dosenbach160}
\alias{dosenbach160}
\docType{data}
\title{Coordinates for data from the Dosenbach160 atlas}
\description{
  This is a list of spatial coordinates for the Dosenbach160 atlas, along with
  indices for the major lobes of the brain and functional networks they specify
  in their manuscript.
}
\usage{data("dosenbach160")}
\format{
  A data frame with 160 observations on the following 8 variables.
  \describe{
    \item{\code{name}}{a character vector of region names}
    \item{\code{x.mni}}{a numeric vector of x-coordinates (in MNI space)}
    \item{\code{y.mni}}{a numeric vector of y-coordinates (in MNI space)}
    \item{\code{z.mni}}{a numeric vector of z-coordinates (in MNI space)}
    \item{\code{lobe}}{a factor with levels \code{Frontal} \code{Parietal} \code{Temporal} \code{Occipital} \code{Insula} \code{Cingulate} \code{SCGM} \code{Cerebellum}}
    \item{\code{hemi}}{a factor with levels \code{L} \code{R} \code{B}}
    \item{\code{index}}{a numeric vector}
    \item{\code{network}}{a factor with levels \code{default} \code{fronto-parietal} \code{cingulo-opercular} \code{sensorimotor} \code{cerebellum} \code{occipital}}
  }
}
\references{
  Dosenbach, N. U., Nardos, B., Cohen, A. L., Fair, D. A., Power, J. D., Church,
  J. A., Nelson, S.M., Wig, G.S., Vogel, A.C., Lessov-Schlaggar, C.N., Barnes,
  K. A. (2010). \emph{Prediction of individual brain maturity using fMRI}.
  Science, 329(5997), 1358-1361.
}
\examples{
data(dosenbach160)
str(dosenbach160)
}
\keyword{datasets}
