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
  V(g.new)$color.lobe <- group.cols[V(g.new)$lobe]
  E(g.new)$color.lobe <- color.edges(g.new, V(g.new)$lobe)
  x <- y <- z <- x.mni <- y.mni <- z.mni <- NULL
  atlas.dt <- eval(parse(text=g.new$atlas))
  vorder <- match(V(g.new)$name, atlas.dt$name)
  V(g.new)$x <- atlas.dt[vorder, x]
  V(g.new)$y <- atlas.dt[vorder, y]
  V(g.new)$z <- atlas.dt[vorder, z]
  V(g.new)$x.mni <- atlas.dt[vorder, x.mni]
  V(g.new)$y.mni <- atlas.dt[vorder, y.mni]
  V(g.new)$z.mni <- atlas.dt[vorder, z.mni]

  return(g.new)
}
