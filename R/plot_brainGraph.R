#' Plot a graph with a specific spatial layout
#'
#' This function plots a graph when the spatial layout of the nodes is important
#' (e.g. in the brain). It is really just a wrapper for
#' \code{\link{plot.igraph}}, with some options pre-specified that work for
#' plotting in the brain's layout. This means that
#' \code{\link{set.brainGraph.attributes}} needs to be run on the graph,
#' and a valid set of coordinates provided for the vertices. Most of the
#' parameters valid here can be seen in \code{\link{igraph.plotting}}.
#'
#' @param g The graph to plot
#' @param rescale A logical, whether to rescale the coordinates
#' @param ylim A vector giving limits for the vertical axis
#' @param asp A numeric constant for the aspect ratio
#' @param main Character string for the main title
#' @param ... Other parameters (passed to \code{\link{plot}}).
#' @export

plot_brainGraph <- function(g, rescale=F, ylim=c(-1.5, 1.5),
                            asp=0, main=NULL, ...) {

  plot(g, asp=asp, rescale=rescale, ylim=ylim, ...)

  Nv <- vcount(g)
  Ne <- ecount(g)
  g.density <- round(graph.density(g), digits=3)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
  subt <- paste('# vertices: ', Nv, '# edges: ', Ne, '\n',
                'Density: ', g.density)
  title(main=main, sub=subt, col.main='white', col.sub='white')
}
