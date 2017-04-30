#' Create an empty graph with attributes for brainGraph
#'
#' This function creates an empty undirected graph with vertex count equal to
#' the atlas specified, and includes some graph-, vertex-, and
#' edge-level attributes that are important for \code{brainGraph} functions.
#' Basically a wrapper for \code{\link[igraph]{make_empty_graph}}.
#'
#' The input graph must have the graph attribute \emph{atlas} and vertex
#' attribute \emph{name} already.
#'
#' @param atlas Character string of the atlas to create a graph from
#' @export
#'
#' @return An empty \code{igraph} graph object with several additional
#'   attributes
#' @seealso \code{\link[igraph]{make_empty_graph}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

make_empty_brainGraph <- function(atlas) {
  atlas <- match.arg(atlas, choices=data(package='brainGraph')$results[, 3])
  atlas.dt <- eval(parse(text=atlas))
  g.new <- make_empty_graph(nrow(atlas.dt), directed=FALSE)
  g.new$atlas <- atlas
  V(g.new)$name <- atlas.dt$name
  g.new <- assign_lobes(g.new)

  return(g.new)
}
