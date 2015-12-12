#' Calculate an asymmetry index based on edge counts
#'
#' This function will calculate an asymmetry index that is a measure of whether
#' or not more edges are present in the left or right hemisphere of a graph for
#' brain MRI data. You can choose a value for each vertex, or for the whole
#' hemisphere.
#'
#' The equation is:
#'\deqn{A = \frac{E_{lh} - E_{rh}}{0.5 \times (E_{lh} + E_{rh})}}
#' where \emph{lh} and \emph{rh} are left and right hemispheres, respectively.
#' The range of this measure is \eqn{[-2, 2]} (although the limits will only be
#' reached if all edges are in one hemisphere), with negative numbers
#' indicating more edges in the right hemisphere, and a value of 0 indicating
#' equal number of edges in each hemisphere.
#'
#' @param g The igraph graph object
#' @param level A character string indicating whether to calculate asymmetry for
#' each region, or the hemisphere as a whole (default: 'hemi')
#' @param .parallel Logical indicating whether or not to use \emph{foreach}
#' (default: TRUE)
#' @export
#'
#' @return A data table with edge counts for both hemispheres and the asymmetry
#' index; if \emph{level} is 'vertex', the data table will have \emph{vcount(g)}
#' rows.
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

edge_asymmetry <- function(g, level=c('hemi', 'vertex'), .parallel=TRUE) {
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if (!'hemi' %in% vertex_attr_names(g)) {
    stop(sprintf('Graph "%s" does not have vertex attribute "hemi"',
                 deparse(substitute(g))))
  }
  i <- region <- NULL
  level <- match.arg(level)
  if (level == 'hemi') {
    lh <- length(E(g)[which(V(g)$hemi == 'L') %--% which(V(g)$hemi == 'L')])
    rh <- length(E(g)[which(V(g)$hemi == 'R') %--% which(V(g)$hemi == 'R')])
    asymm <- data.table(region='all', lh=lh, rh=rh)

  } else if (level == 'vertex') {
    inds <- which(V(g)$degree > 0)
    asymm <- data.frame(lh=rep(0, vcount(g)), rh=rep(0, vcount(g)))
    if (isTRUE(.parallel)) {
      asymm[inds, ] <- foreach (i=inds, .combine='rbind') %dopar% {
        lh <- length(E(g)[i %--% which(V(g)$hemi == 'L')])
        rh <- length(E(g)[i %--% which(V(g)$hemi == 'R')])
        c(lh, rh)
      }
    } else {
      asymm$lh[inds] <- vapply(V(g)[inds], function(x)
                               length(E(g)[x %--% which(V(g)$hemi == 'L')]),
                               numeric(1))
      asymm$rh[inds] <- vapply(V(g)[inds], function(x)
                               length(E(g)[x %--% which(V(g)$hemi == 'R')]),
                               numeric(1))
    }
    setDT(asymm)
    asymm[, region := V(g)$name]
  }
  asymm[, asymm := 2 * (lh - rh) / (lh + rh)]
  return(asymm)
}
