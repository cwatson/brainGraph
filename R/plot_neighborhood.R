#' Plot the neighborhood of a given vertex.
#'
#' This function plots only the neighborhood of a single given vertex.
#'
#' @param g the total adjacency matrix
#' @param v character string of the vertex name whose neighborhood is to be
#'  plotted
#' @export

plot.neighborhood <- function(g, v) {
  sub.g <- graph.neighborhood(g, nodes=v, order=1)[[1]]

  plot.over.brain(0)
  plot.adj(sub.g,
           layout=as.matrix(coords[V(sub.g)$name, 1:2]) %*% diag(c(1/8, 1/10)),
           rescale=F)
  title(paste('Neighborhood of', v))
}
