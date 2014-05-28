#' Plots an adjacency matrix with a specific spatial layout.
#'
#' This function plots an adjacency matrix when the spatial layout of the nodes
#' is important. It is really just a wrapper for \code{\link{plot}}, but it uses
#' as a default spatial layout coordinates that work for the brain. This means
#' that \code{\link{set.brainGraph.attributes}} needs to be run on the graph,
#' and a valid set of coordinates provided for the vertices. Most of the
#' parameters here are for \code{\link{plot.igraph}}.
#'
#' @param adj.graph The adjacency graph to plot
#' @param vertex.size The size of the vertices
#' @param vertex.color The color of the vertices
#' @param vertex.label.cex The scaling of vertex label text
#' @param edge.color The color of the edges
#' @param edge.width The size of the edges
#' @param rescale A logical constant, whether to rescale the coordinates
#' @param ylim A vector giving limits for the vertical axis
#' @param ... Other parameters (passed to \code{\link{plot}}).
#' @export

plot.adj <- function(adj.graph, vertex.size=10, vertex.color='lightblue3',
                    vertex.label.cex=0.75, edge.color='red', edge.width=2,
                    rescale=F, ylim=c(-1.5, 1.5), ...) {

  plot(adj.graph, vertex.size=vertex.size,
    vertex.color=vertex.color, vertex.label.cex=vertex.label.cex,
    edge.width=edge.width, edge.color=edge.color,
    asp=0, rescale=rescale, ylim=ylim, ...)
}
