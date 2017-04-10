#' Create a graph of the union of multiple vertex neighborhoods
#'
#' This function accepts multiple vertices, creates graphs of their
#' neighborhoods (of order 1), and returns the union of those graphs.
#'
#' @param g An \code{igraph} graph object
#' @param vs Either a character or integer vector (vertex names or indices,
#' respectively) for the vertices of interest
#' @export
#'
#' @return An \code{igraph} graph object containing the union of all edges and
#'   vertices in the neighborhoods of the input vertices; only the vertex
#'   attribute \emph{name} will be present
#' @seealso \code{\link[igraph]{make_ego_graph}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' subg <- graph_neighborhood_multiple(g1[[N]], c(24, 58))
#' subg <- graph_neighborhood_multiple(g1[[N]], c('lPCUN', 'rPCUN'))
#' }

make_ego_brainGraph <- function(g, vs) {

  subgs <- make_ego_graph(g, order=1, nodes=vs)
  if (is.character(vs)) vs <- which(V(g)$name %in% vs)

  for (i in seq_along(vs)) {
    subgs[[i]] <- delete_all_attr(subgs[[i]], keep.names=TRUE)
  }

  combine_graphs <- function(x, y) {
    n <- length(x)
    if (n < 2) {
      res <- x[[1]] %u% y
    } else {
      y <- x[[n]] %u% y
      x <- x[-n]
      res <- combine_graphs(x, y)
    }
    return(res)
  }

  inds <- unique(c(vs, unlist(lapply(vs, function(x) neighbors(g, x)))))
  subg.all <- combine_graphs(subgs, make_empty_graph(directed=F) + vertices(V(g)$name[inds]))
  return(subg.all)
}
