#' Calculate Euclidean distance between vertices (MNI space)
#'
#' This function calculates the Euclidean distance of edges between vertices of
#' a graph. The distances are in mm and based on MNI space. The distances are
#' \emph{NOT} along the cortical surface, so can only be considered
#' approximations, particularly concerning inter-hemispheric connections. The
#' input graph must have \emph{atlas} as a graph-level attribute.
#'
#' @param g An igraph graph object
#' @export
#'
#' @return A numeric vector with length equal to the edge count of the input
#' graph
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

spatial.dist <- function(g) {
  if (! 'atlas' %in% graph_attr_names(g)) {
    stop(sprintf('Error: Input graph "%s" does not have an "atlas" attribute!',
                 deparse(substitute(g))))
  }
  name <- index <- NULL
  atlas.dt <- eval(parse(text=data(list=g$atlas)))
  coords <- atlas.dt[, c('x.mni', 'y.mni', 'z.mni'), with=F]
  es <- get.edgelist(g)
  dists <- sqrt(rowSums((coords[atlas.dt[name == es[, 2], index], ] -
                         coords[atlas.dt[name == es[, 1], index], ])^2))
}
