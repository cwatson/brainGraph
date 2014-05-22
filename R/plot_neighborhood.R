#' Plot the neighborhood of a given vertex.
#'
#' This function plots only the neighborhood of a single given vertex, if you
#' want to see e.g. the neighborhood of the left precuneus and nothing else.
#'
#' @param g The total adjacency matrix
#' @param v Character string of the vertex name whose neighborhood is to be
#'  plotted
#' @export

plot.neighborhood <- function(g, v) {
  sub.g <- graph.neighborhood(g, nodes=v, order=1)[[1]]

  plot.over.brain(0)
  plot.adj(sub.g,
           layout=as.matrix(coords.cur[V(sub.g)$name, 1:2]) %*% coord.scale,
           rescale=F)
  par(new=T, mar=c(5, 0, 3, 0)+0.1)
  title(paste('Neighborhood of', v))
}
