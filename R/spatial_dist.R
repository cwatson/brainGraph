#' Calculate Euclidean distance of edges
#'
#' This function calculates the Euclidean distance of edges between vertices of
#' a graph. The distances are in \emph{mm} and based on MNI space. The distances
#' are \emph{NOT} along the cortical surface, so can only be considered
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

edge_spatial_dist <- function(g) {
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

#' Calculate average Euclidean distance for each vertex
#'
#' This function calculates, for each vertex of a graph, the average Euclidean
#' distance across all of that vertex's connections.
#'
#' @param g An \code{igraph} graph object
#' @export
#'
#' @return A named numeric vector of average distance (in \emph{mm})
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Alexander-Bloch A.F., Vertes P.E., Stidd R. et al. (2013)
#'   \emph{The anatomical distance of functional connections predicts brain
#'   network topology in health and schizophrenia}. Cerebral Cortex, 23:127-138.

vertex_spatial_dist <- function(g) {
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if (! 'dist' %in% edge_attr_names(g)) {
    stop(sprintf('Input graph "%s" does not have a "dist" edge attribute',
                 deparse(substitute(g))))
  }

  Nv <- vcount(g)
  dists <- sapply(V(g)$name, function(x) mean(E(g)[x %--% seq_len(Nv)]$dist))
  return(dists)
}
