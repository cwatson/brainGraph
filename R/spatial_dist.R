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
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if (! 'atlas' %in% graph_attr_names(g)) {
    stop(sprintf('Input graph "%s" does not have an "atlas" attribute',
                 deparse(substitute(g))))
  }
  name <- x.mni <- y.mni <- z.mni <- NULL
  coords <- eval(parse(text=g$atlas))[, list(name, x.mni, y.mni, z.mni)]
  setkey(coords, name)
  es <- get.edgelist(g)
  dists <- sqrt(rowSums((coords[es[, 2], list(x.mni, y.mni, z.mni)] -
                         coords[es[, 1], list(x.mni, y.mni, z.mni)])^2))
}
