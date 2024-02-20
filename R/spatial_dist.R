#' Calculate Euclidean distance of edges and vertices
#'
#' \code{edge_spatial_dist} calculates the Euclidean distance of an
#' \code{igraph} graph object's edges. The distances are in \emph{mm} and based
#' on MNI space. These distances are \emph{NOT} along the cortical surface, so
#' can only be considered approximations, particularly concerning
#' inter-hemispheric connections. The input graph must have \emph{atlas} as a
#' graph-level attribute.
#'
#' @param g An \code{igraph} graph object
#' @export
#'
#' @return \code{edge_spatial_dist} - a numeric vector with length equal to the
#'   edge count of the input graph, consisting of the Euclidean distance (in
#'   \emph{mm}) of each edge
#'
#' @name Graph Distances
#' @rdname spatial_dist
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

edge_spatial_dist <- function(g) {
  stopifnot(is_igraph(g), 'atlas' %in% graph_attr_names(g))

  name <- x.mni <- y.mni <- z.mni <- NULL
  coords <- get(g$atlas)[, list(name, x.mni, y.mni, z.mni)]
  setkey(coords, name)
  es <- as_edgelist(g)
  dists <- sqrt(rowSums((coords[es[, 2L], list(x.mni, y.mni, z.mni)] -
                         coords[es[, 1L], list(x.mni, y.mni, z.mni)])^2))
}

#' Calculate average Euclidean distance for each vertex
#'
#' \code{vertex_spatial_dist} calculates, for each vertex of a graph, the
#' average Euclidean distance across all of that vertex's connections.
#'
#' @inheritParams edge_spatial_dist
#' @export
#' @importFrom Matrix colSums
#'
#' @return \code{vertex_spatial_dist} - a named numeric vector with length equal
#'   to the number of vertices, consisting of the average distance (in
#'   \emph{mm}) for each vertex
#'
#' @rdname spatial_dist
#' @references Alexander-Bloch, A.F. and Vertes, P.E. and Stidd, R. et al.
#'   (2013) The anatomical distance of functional connections predicts brain
#'   network topology in health and schizophrenia. \emph{Cerebral Cortex},
#'   \bold{23}, 127--138. \doi{10.1093/cercor/bhr388}

vertex_spatial_dist <- function(g) {
  stopifnot(is_igraph(g), 'dist' %in% edge_attr_names(g))

  A <- as_adj(g, names=FALSE, attr='dist')
  dists <- colSums(A) / colSums(A != 0)
  return(dists)
}
