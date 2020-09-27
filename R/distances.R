#' Calculate weighted shortest path lengths
#'
#' Calculate graph or vertex average shortest path lengths. For vertices, this
#' is just the row means of the distance matrix. For the graph-level, it is the
#' overall mean of the distance matrix.
#'
#' By default, edge weights are not transformed (e.g., inverted). However, if
#' set to \code{TRUE}, then the input graph must have a graph-level attribute
#' called \code{'xfm.type'} or you must supply a value in the function call. If
#' you supply a distance matrix (the \code{D} argument), it is not necessary to
#' transform edge weights, as it is assumed the the distance matrix was
#' calculated from a graph with transformed edge weights already.
#'
#' @param level Character string indicating whether to calculate vertex- or
#'   graph-level shortest path length. Default: \code{'graph'}
#' @inheritParams efficiency
#' @inheritParams xfm.weights
#' @export
#' @return Numeric vector (if \code{level='vertex'}) of each vertex's shortest
#'   path length, or a single number for the graph-level average

mean_distance_wt <- function(g, level=c('graph', 'vertex'), weights=NULL,
                             xfm=FALSE, xfm.type=NULL, D=NULL) {
  if (isTRUE(xfm)) {
    if (is.null(xfm.type)) {
      xfm.type <- if ('xfm.type' %in% graph_attr_names(g)) g$xfm.type else '1/w'
    }
    g <- xfm.weights(g, xfm.type)
  }
  weights <- check_weights(g, weights)
  level <- match.arg(level)
  if (is.null(D)) D <- distances(g, weights=weights)
  D[is.infinite(D)] <- diag(D) <- NA
  Lp <- switch(level,
               vertex=rowMeans(D, na.rm=TRUE),
               graph=mean(D, na.rm=TRUE))
  Lp[is.nan(Lp)] <- 0
  return(Lp)
}
