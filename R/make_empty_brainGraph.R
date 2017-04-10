#' Create an empty graph with attributes for brainGraph
#'
#' This function creates an empty graph and includes some graph-, vertex-, and
#' edge-level attributes that are important for \code{brainGraph} functions.
#' Basically a wrapper for \code{\link[igraph]{make_empty_graph}}.
#'
#' The input graph must have the graph attribute \emph{atlas} and vertex
#' attribute \emph{name} already.
#'
#' @param g An \code{igraph} graph object
#'
#' @return An empty \code{igraph} graph object with several additional
#'   attributes
#' @seealso \code{\link[igraph]{make_empty_graph}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

make_empty_brainGraph <- function(g) {
  g.new <- make_empty_graph(vcount(g), directed=FALSE)
  g.new$atlas <- g$atlas
  V(g.new)$name <- V(g)$name
  g.new <- assign_lobes(g.new)

  return(g.new)
}
