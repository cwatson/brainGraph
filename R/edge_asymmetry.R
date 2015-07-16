#' Calculate an asymmetry index based on edge counts
#'
#' This function will calculate an asymmetry index that is a measure of whether
#' or not more edges are present in the left or right hemisphere of an igraph
#' graph object for brain MRI data.
#'
#' The equation is:
#'\deqn{A = \frac{E_{lh} - E_{rh}}{0.5 \times (E_{lh} + E_{rh})}}
#' where \emph{lh} and \emph{rh} are left and right hemispheres, respectively.
#' The range of this measure is \eqn{[-1, 1]}, with negative numbers
#' indicating more edges in the left hemisphere, and a value of 0 indicating
#' equal number of edges in each hemisphere.
#'
#' @param g The igraph graph object
#' @export
#'
#' @return A data table with edge counts for both hemispheres and the asymmetry
#' index
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

edge_asymmetry <- function(g) {
  if (!'hemi' %in% vertex_attr_names(g)) {
    stop(sprintf('Graph "%s" does not have vertex attribute "hemi"!',
                 deparse(substitute(g))))
  }
  lh <- length(E(g)[which(V(g)$hemi == 'L') %--% which(V(g)$hemi == 'L')])
  rh <- length(E(g)[which(V(g)$hemi == 'R') %--% which(V(g)$hemi == 'R')])
  asymm <- (lh - rh) / .Internal(mean(c(lh, rh)))

  return(data.table(lh=lh, rh=rh, asymm=asymm))
}
