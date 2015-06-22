#' Analysis of network robustness
#'
#' This function performs a "targeted attack" of a graph or a "random failure"
#' analysis, calculating the size of the largest component after edge or vertex
#' removal.
#'
#' In a targeted attack, it will sort the vertices by either degree or
#' betweenness centrality (or sort edges by betweenness), and successively
#' remove the top vertices/edges. Then it calculates the size of the largest
#' component.
#'
#' In a random failure analysis, vertices/edges are removed in a random order.
#'
#' @param g The igraph graph object of interest
#' @param type A character string; either 'vertex' or 'edge' removals
#' @param measure A character string; sort by either 'btwn.cent' or 'degree', or
#' choose 'random'
#' @export
#'
#' @return A vector representing the ratio of maximal component size after each
#' removal to the graph's original maximal component
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Albert R., Jeong H., Barabasi A. (2000) \emph{Error and attack
#' tolerance of complex networks}. Nature, 406:378-381.

robustness <- function(g, type=c('vertex', 'edge'),
                       measure=c('btwn.cent', 'degree', 'random')) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  if (type == 'vertex') {
    n <- vcount(g)

    if (measure == 'random') {
      ord <- V(g)$name[sample(n)]
    } else {
      ord <- V(g)$name[order(vertex_attr(g, measure), decreasing=T)]
    }
    max.comp <- vector('integer', length=n+1)
    max.comp[1] <- g$conn.comp[1, 1]
    for (i in seq_len(n - 1)) {
      g <- delete.vertices(g, ord[i])
      max.comp[i+1] <- max(components(g)$csize)
    }

  } else {
    if (measure == 'degree') {
      stop('For edge attacks, must choose "btwn.cent" or "random"!')
    } else if (measure == 'random') {
      m <- ecount(g)
      ord <- sample(m)
    } else {
      ord <- order(E(g)$btwn, decreasing=T)
    }

    verts <- as_edgelist(g)[ord, ]
    max.comp <- vector('integer', length=length(ord)+1)
    max.comp[1] <- g$conn.comp[1, 1]
    for (i in seq_along(ord)) {
      g <- delete.edges(g, E(g)[verts[ord[i], 1] %--% verts[ord[i], 2]])
      max.comp[i+1] <- max(components(g)$csize)
    }
  }
  return(max.comp / max.comp[1])
}
