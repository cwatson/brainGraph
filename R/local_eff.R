#' Calculate the local efficiency of a graph.
#'
#' This function calculates the local efficiency of each node of a graph.
#'
#' @param g the adjacency graph on which to calculate efficiency
#' @export
#'
#' @return A vector of the local efficiencies for each node of the graph.

local.eff <- function(g) {
  eff <- numeric(vcount(g))
  nodes <- which(degree(g) != 0)
  eff[nodes] <- vapply(nodes,
      function(x) global.eff(graph.neighborhood(g, nodes=x, order=1)[[1]]),
      numeric(1))
  return(eff)
}
