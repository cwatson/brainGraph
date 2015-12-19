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
#' @param N Integer; the number of iterations if \emph{random} is chosen
#' @export
#'
#' @return A vector representing the ratio of maximal component size after each
#' removal to the graph's original maximal component
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Albert R., Jeong H., Barabasi A. (2000) \emph{Error and attack
#' tolerance of complex networks}. Nature, 406:378-381.

robustness <- function(g, type=c('vertex', 'edge'),
                       measure=c('btwn.cent', 'degree', 'random'), N=1e3) {
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  type <- match.arg(type)
  measure <- match.arg(measure)
  max.comp.orig <- g$max.comp
  if (type == 'vertex') {
    n <- vcount(g)

    if (measure == 'random') {
      max.comp <- matrix(nrow=n+1, ncol=N)
      max.comp <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        g.new <- g
        ord <- V(g.new)$name[sample(n)]
        tmp <- vector('integer', length=n)
        for (j in seq_len(n - 1)) {
          g.new <- delete.vertices(g.new, ord[j])
          tmp[j] <- max(components(g.new)$csize)
        }
        tmp
      }
      max.comp <- cbind(1, max.comp)
      max.comp.removed <- colMeans(max.comp)

    } else {
      ord <- V(g)$name[order(vertex_attr(g, measure), decreasing=T)]
      max.comp.removed <- vector('integer', length=n+1)
      max.comp.removed[1] <- 1
      for (i in seq_len(n - 1)) {
        g <- delete.vertices(g, ord[i])
        max.comp.removed[i+1] <- max(components(g)$csize)
      }
    }

  } else {
    m <- ecount(g)
    if (measure == 'degree') {
      stop('For edge attacks, must choose "btwn.cent" or "random"!')
    } else if (measure == 'random') {
      max.comp <- matrix(nrow=m+1, ncol=N)
      max.comp <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        g.new <- g
        ord <- sample(m)
        verts <- as_edgelist(g.new)[ord, ]
        tmp <- vector('integer', length=m)
        for (j in seq_along(ord)) {
          g.new <- delete.edges(g.new, E(g.new)[verts[ord[j], 1] %--% verts[ord[j], 2]])
          tmp[j] <- max(components(g.new)$csize)
        }
        tmp
      }
      max.comp <- cbind(1, max.comp)
      max.comp.removed <- colMeans(max.comp)

    } else {
      ord <- order(E(g)$btwn, decreasing=T)
      verts <- as_edgelist(g)[ord, ]
      max.comp.removed <- vector('integer', length=m+1)
      max.comp.removed[1] <- 1
      for (i in seq_along(ord)) {
        g <- delete.edges(g, E(g)[verts[ord[i], 1] %--% verts[ord[i], 2]])
        max.comp.removed[i+1] <- max(components(g)$csize)
      }
    }

  }
  comp.ratio <- c(max.comp.removed[1], max.comp.removed[-1] / max.comp.orig)
  return(comp.ratio)
}
