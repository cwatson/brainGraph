#' Calculate the rich club of a graph.
#'
#' This function calculates the rich club of a graph, both the coefficient
#' phi and the nodes that make up this subgraph.
#'
#' @param g The graph of interest
#' @param k The minimum degree for including a node
#' @export
#'
#' @return A list with the following components:
#' \item{coeff}{The rich club coefficient, phi.}
#' \item{graph}{A subgraph containing only the rich club nodes.}
#' \item{Nk}{The number of vertices in the rich club graph.}
#' \item{Ek}{The number of edges in the rich club graph.}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Zhou S., Mondragon R.J. (2004) \emph{The rich-club phenomenon
#' in the internet topology}. IEEE Comm Lett, 8:180-182.

rich.club.coeff <- function(g, k=1) {
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  Nv <- vcount(g)
  Nk <- sum(degs > k)
  if (Nk == 0) {
    list(coeff=NaN, graph=graph.empty(), Nk=0, Ek=0)
  } else {
    rich.club.nodes <- order(degs)[(Nv - Nk + 1):Nv]
    rich.club.graph <- induced.subgraph(g, rich.club.nodes)
    Ek <- ecount(rich.club.graph)
    phi <- graph.density(rich.club.graph)
    list(coeff=phi, graph=rich.club.graph, Nk=Nk, Ek=Ek)
  }
}
