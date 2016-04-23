#' Calculate the rich club of a graph
#'
#' This function calculates the rich club of a graph, both the coefficient
#' \eqn{\phi} and the nodes that make up this subgraph.
#'
#' @param g The graph of interest
#' @param k The minimum degree for including a vertex (default: 1)
#' @param weighted A logical indicating whether or not edge weights should be
#' used (default: FALSE)
#' @export
#'
#' @return A list with the following components:
#' \item{phi}{The rich club coefficient, \eqn{\phi}.}
#' \item{graph}{A subgraph containing only the rich club nodes.}
#' \item{Nk}{The number of vertices in the rich club graph.}
#' \item{Ek}{The number of edges in the rich club graph.}
#'
#' @seealso \code{\link{rich.club.norm}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Zhou S., Mondragon R.J. (2004) \emph{The rich-club phenomenon
#' in the internet topology}. IEEE Comm Lett, 8:180-182.
#' @references Opsahl T., Colizza V., Panzarasa P., Ramasco J.J. (2008)
#' \emph{Prominence and control: the weighted rich-club effect}. Physical Review
#' Letters, 101.16:168702.

rich.club.coeff <- function(g, k=1, weighted=FALSE) {
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  Nv <- vcount(g)
  Nk <- sum(degs > k)
  if (Nk == 0) {
    return(list(phi=NaN, graph=make_empty_graph(), Nk=0, Ek=0))
  } else {
    rich.club.nodes <- order(degs)[(Nv - Nk + 1):Nv]
    rich.club.graph <- induced.subgraph(g, rich.club.nodes)
    Ek <- ecount(rich.club.graph)

    if (isTRUE(weighted)) {
      Wr <- sum(E(rich.club.graph)$weight)
      weights <- sort(E(g)$weight, decreasing=T)[1:Ek]
      phi <- Wr / sum(weights)
    } else {
      phi <- graph.density(rich.club.graph)
    }

    return(list(phi=phi, graph=rich.club.graph, Nk=Nk, Ek=Ek))
  }
}
