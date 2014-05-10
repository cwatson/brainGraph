#' Calculate the rich club of a graph.
#'
#' This function calculates the rich club of a graph, both the coefficient
#' phi and the nodes that make up this subgraph.
#'
#' @param graph the adjacency graph of interest
#' @param k the minimum degree for including a node
#' @export
#'
#' @return A list with the following components:
#' \item{coeff}{The rich club coefficient, phi.}
#' \item{graph}{A subgraph containing only the rich club nodes.}
#' \item{Nk}{The number of vertices in the rich club graph.}
#' \item{Ek}{The number of edges in the rich club graph.}

rich.club.coeff <- function(graph, k) {
  Nv <- vcount(graph)
  Nk <- sum(degree(graph) > k)
  if (Nk == 0) {
    list(coeff=NaN, graph=graph.empty(), Nk=0, Ek=0)
  } else {
    rich.club.nodes <- sort(degree(graph), index.return=T)$ix[(Nv - Nk + 1):Nv]
    rich.club.graph <- induced.subgraph(graph, rich.club.nodes)
    Ek <- ecount(rich.club.graph)
    phi <- (2 * Ek) / (Nk * (Nk - 1))
    list(coeff=phi, graph=rich.club.graph, Nk=Nk, Ek=Ek)
  }
}
