#' Plot the neighborhood of a given vertex.
#'
#' This function plots only the neighborhood of a single given vertex, if you
#' want to see e.g. the neighborhood of the left precuneus and nothing else.
#'
#' @param g The graph to be plotted
#' @param v Character string of the vertex name whose neighborhood is to be
#'  plotted
#' @param ... Other parameters (passed to \code{\link{plot.adj}})
#' @export

plot.neighborhood <- function(g, v, ...) {
  g.sub <- graph.neighborhood(g, nodes=v, order=1)[[1]]

  args <- list(...)
  if (hasArg(vertex.label)) {
    vertex.label <- args[['vertex.label']]
    if (!is.na(vertex.label)) {
      vertex.label <- V(g.sub)$name
      vertex.label.cex <- 0.75
    } else {
      vertex.label.cex <- NA
    }
  } else {
    vertex.label <- NA
    vertex.label.cex <- NA
  }

  if (hasArg(vertex.color)) {
    vertex.color <- args[['vertex.color']]
    edge.color <- args[['edge.color']]
    if (length(vertex.color) > 1) {
      vertex.color <- V(g.sub)$color
      edge.color <- E(g.sub)$color
    }
  } else {
    vertex.color <- 'lightblue3'
    edge.color <- 'red'
  }

  plot.adj(g.sub,
           vertex.color=vertex.color, edge.color=edge.color,
           vertex.label=vertex.label, vertex.label.cex=0.75, ...)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
  title(paste('Neighborhood of', v))
}
