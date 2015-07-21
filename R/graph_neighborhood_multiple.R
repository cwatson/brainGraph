#' Take the union of multiple neighborhood graphs
#'
#' This function takes multiple vertices, creates graphs of their neighborhoods
#' (of order 1), and takes the union of those graphs.
#'
#' @param g The igraph graph object
#' @param vs Either a character or integer vector (vertex names or indices,
#' respectively) for the vertices of interest
#' @export
#'
#' @return An igraph graph object containing the union of all edges and vertices
#' in the neighborhoods of the input vertices; only the vertex attribute
#' \emph{name} will be present
#' @seealso \code{\link[igraph]{make_ego_graph}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' subg <- graph_neighborhood_multiple(g1[[N]], c(24, 58))
#' subg <- graph_neighborhood_multiple(g1[[N]], c('lPCUN', 'rPCUN'))
#' }

graph_neighborhood_multiple <- function(g, vs) {

  subgs <- make_ego_graph(g, order=1, nodes=vs)
  if (is.character(vs)) {
    vs <- which(V(g)$name %in% vs)
  }

  inds <- unique(c(vs, unlist(lapply(vs, function(x) neighbors(g, x)))))

  for (i in seq_along(vs)) {
    for (att in graph_attr_names(subgs[[i]])) {
      subgs[[i]] <- delete_graph_attr(subgs[[i]], att)
    }
    for (att in vertex_attr_names(subgs[[i]])[!grepl('name', vertex_attr_names(subgs[[i]]))]) {
      subgs[[i]] <- delete_vertex_attr(subgs[[i]], att)
    }
    for (att in edge_attr_names(subgs[[i]])) {
      subgs[[i]] <- delete_edge_attr(subgs[[i]], att)
    }
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

  subg.all <- combine_graphs(subgs, make_empty_graph(directed=F) + vertices(V(g)$name[inds]))
  return(subg.all)
}
