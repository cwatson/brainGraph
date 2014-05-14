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
  nodes <- which(degree(g) > 1)
  eff[nodes] <- vapply(nodes,
    function(x) global.eff(induced.subgraph(g, neighbors(g, v=x))),
      numeric(1))
  return(eff)
}
