#' Plots an adjacency matrix with a specific spatial layout.
#'
#' This function plots an adjacency matrix when the spatial layout of the nodes
#' is important. It is really just a wrapper for \code{\link{plot}}, but it uses
#' as a default spatial layout coordinates that work for the brain. This means
#' that \code{\link{set.brainGraph.attributes}} needs to be run on the graph,
#' and a valid set of coordinates provided for the vertices. Most of the
#' parameters here are for \code{\link{plot.igraph}}.
#'
#' @param g The graph to plot
#' @param vertex.size The size of the vertices
#' @param edge.width The size of the edges
#' @param rescale A logical constant, whether to rescale the coordinates
#' @param ylim A vector giving limits for the vertical axis
#' @param ... Other parameters (passed to \code{\link{plot}}).
#' @export

plot.adj <- function(g, vertex.size=10, edge.width=2,
                    rescale=F, ylim=c(-1.5, 1.5), ...) {

  plot(g, vertex.size=vertex.size,
    edge.width=edge.width,
    asp=0, rescale=rescale, ylim=ylim, ...)

  if (!is.null(g$density)) {
    subt <- sprintf('%s: %1.3f', 'Density', g$density)
    par(new=T, mar=c(5, 0, 3, 0)+0.1)
    title(sub=subt, col.sub='white')
  }
}
