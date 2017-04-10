#' Color graph vertices and edges
#'
#' \code{set_vertex_color} takes an integer vector (representing membership of
#' a community or component, etc.) and creates a character vector of colors for
#' each community/module, component, etc. This only assigns a color to groups
#' with at least 2 members; isolated vertices will be colored \emph{gray}.
#'
#' @param memb An integer vector representing membership of e.g. a community
#'
#' @return \code{set_vertex_color} - a character vector of colors with length
#'   equal to the number of vertex groups
#'
#' @name GraphColors
#' @aliases set_vertex_color
#' @rdname color_vertices_edges

set_vertex_color <- function(memb) {
  big.modules <- which(as.integer(table(memb)) > 1)

  mod.colors.memb <- rep('gray', length=max(memb))
  mod.colors.memb[big.modules] <- group.cols[big.modules]

  return(mod.colors.memb)
}

#' Color graph edges
#'
#' \code{set_edge_color} additionally takes as input an \code{igraph} graph
#' object. It assigns to the edges a specific color (the same as the vertex
#' membership colors). Edges that connect vertices of two different groups are
#' colored gray.
#'
#' @param g An \code{igraph} graph object
#' @inheritParams set_vertex_color
#'
#' @return \code{set_edge_color} - a character vector of colors with length
#'   equal to the edge count of \code{g}
#'
#' @aliases set_edge_color
#' @rdname color_vertices_edges

set_edge_color <- function(g, memb) {
  stopifnot(length(memb) == vcount(g))
  big.modules <- as.integer(names(which(table(memb) > 1)))

  newcols <- rep('gray50', length=ecount(g))
  tmp <- vector('list', length=max(big.modules))
  for (i in big.modules) {
    x <- which(memb == i)
    tmp[[i]] <- as.vector(E(g)[x %--% x])
    if (!is.null(tmp[[i]])) newcols[tmp[[i]]] <- group.cols[i]
  }

  return(newcols)
}
