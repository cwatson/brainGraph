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
#' @param g An \code{igraph} graph object
#' @param type Character string; either 'vertex' or 'edge' removals (default:
#'   \code{vertex})
#' @param measure Character string; sort by either 'btwn.cent' or 'degree', or
#'   choose 'random' (default: \code{btwn.cent})
#' @param N Integer; the number of iterations if \emph{random} is chosen
#'   (default: \code{1e3})
#' @export
#'
#' @return Numeric vector representing the ratio of maximal component size after
#' each removal to the observed graph's maximal component size
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Albert R., Jeong H., Barabasi A. (2000) \emph{Error and attack
#' tolerance of complex networks}. Nature, 406:378-381.

robustness <- function(g, type=c('vertex', 'edge'),
                       measure=c('btwn.cent', 'degree', 'random'), N=1e3) {
  stopifnot(is_igraph(g))
  type <- match.arg(type)
  measure <- match.arg(measure)
  max.comp.orig <- g$max.comp
  if (type == 'vertex') {
    n <- vcount(g)

    if (measure == 'random') {
      max.comp <- matrix(nrow=n+1, ncol=N)
      rand <- matrix(rep(1:n, N), byrow=TRUE, nrow=N)
      index <- t(apply(rand, 1, sample))
      max.comp <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        g.new <- g
        ord <- V(g.new)$name[index[i, ]]
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
      ord <- V(g)$name[order(vertex_attr(g, measure), decreasing=TRUE)]
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
      rand <- matrix(rep(1:m, N), byrow=TRUE, nrow=N)
      index <- t(apply(rand, 1, sample))
      max.comp <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        g.new <- g
        ord <- index[i, ]
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
      ord <- order(E(g)$btwn, decreasing=TRUE)
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
