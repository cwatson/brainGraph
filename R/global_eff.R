#' Calculate the global efficiency of a graph.
#'
#' This function calculates the global efficiency of a graph.
#'
#' @param g The adjacency graph on which to calculate efficiency
#' @export
#'
#' @return A number; the global efficiency of the graph.
#'
#' @references Latora V., Marchiori M. (2001) \emph{Efficient behavior of
#' small-world networks}. Phys Rev Lett, 87.19:198701.

global.eff <- function(g) {
  Nv <- vcount(g)
  nodes <- V(g)[degree(g) != 0]
  paths <- shortest.paths(g, v=nodes, to=nodes)
  paths <- paths[upper.tri(paths)]

  eff <- 2 / Nv / (Nv - 1) * sum(1 / paths[paths != 0])
  return(eff)
}
