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
#' @return Data table with elements:
#'   \item{type}{Character string describing the type of analysis performed}
#'   \item{measure}{The input argument \code{measure}}
#'   \item{comp.pct}{Numeric vector of the ratio of maximal component size after
#'     each removal to the observed graph's maximal component size}
#'   \item{removed.pct}{Numeric vector of the ratio of vertices/edges removed}
#'   \item{Group}{Character string indicating the subject group, if applicable}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Albert R., Jeong H., Barabasi A. (2000) \emph{Error and attack
#' tolerance of complex networks}. Nature, 406:378-381.
#' \url{https://dx.doi.org/10.1515/9781400841356.503}

robustness <- function(g, type=c('vertex', 'edge'),
                       measure=c('btwn.cent', 'degree', 'random'), N=1e3) {
  stopifnot(is_igraph(g))
  group <- NULL
  if ('Group' %in% graph_attr_names(g)) group <- g$Group
  type <- match.arg(type)
  measure <- match.arg(measure)
  if ('max.comp' %in% graph_attr_names(g)) {
    max.comp.orig <- g$max.comp
  } else {
    max.comp.orig <- max(components(g)$csize)
  }
  if (type == 'vertex') {
    n <- vcount(g)
    removed.pct <- seq(0, 1, length=n+1)

    if (measure == 'random') {
      type <- 'Random vertex removal'
      max.comp <- matrix(nrow=n+1, ncol=N)
      rand <- matrix(rep(1:n, N), byrow=TRUE, nrow=N)
      index <- t(apply(rand, 1, sample))
      max.comp <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        g.rand <- g
        ord <- V(g.rand)$name[index[i, ]]
        tmp <- rep(max.comp.orig, n+1)
        for (j in seq_len(n - 1)) {
          g.rand <- delete_vertices(g.rand, ord[j])
          tmp[j+1] <- max(components(g.rand)$csize)
        }
        tmp
      }
      max.comp.removed <- colMeans(max.comp)

    } else {
      if (measure == 'btwn.cent') {
        if (measure %in% vertex_attr_names(g)) {
          val <- vertex_attr(g, measure)
        } else {
          val <- centr_betw(g)$res
        }
      } else if (measure == 'degree') {
        val <- check_degree(g)
      }
      type <- 'Targeted vertex attack'
      ord <- V(g)$name[order(val, decreasing=TRUE)]
      max.comp.removed <- rep(max.comp.orig, n+1)
      for (i in seq_len(n - 1)) {
        g <- delete_vertices(g, ord[i])
        max.comp.removed[i+1] <- max(components(g)$csize)
      }
    }

  } else {
    m <- ecount(g)
    removed.pct <- seq(0, 1, length=m+1)
    if (measure == 'degree') {
      stop('For edge attacks, must choose "btwn.cent" or "random"!')
    } else if (measure == 'random') {
      type <- 'Random edge removal'
      max.comp <- matrix(nrow=m+1, ncol=N)
      rand <- matrix(rep(1:m, N), byrow=TRUE, nrow=N)
      index <- t(apply(rand, 1, sample))
      max.comp <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        g.rand <- g
        ord <- index[i, ]
        el <- as_edgelist(g.rand)[ord, ]
        tmp <- rep(max.comp.orig, m+1)
        for (j in seq_len(m - 1)) {
          g.rand <- graph_from_edgelist(el[-seq_len(j), , drop=FALSE], directed=FALSE)
          tmp[j+1] <- max(components(g.rand)$csize)
        }
        tmp
      }
      max.comp.removed <- colMeans(max.comp)

    } else {
      type <- 'Targeted edge attack'
      ord <- order(E(g)$btwn, decreasing=TRUE)
      el <- as_edgelist(g)[ord, ]
      max.comp.removed <- rep(max.comp.orig, length=m+1)
      for (i in seq_len(m - 1)) {
        g <- graph_from_edgelist(el[-seq_len(i), , drop=FALSE], directed=FALSE)
        max.comp.removed[i+1] <- max(components(g)$csize)
      }
    }

  }
  comp.pct <- max.comp.removed / max.comp.orig
  comp.pct[length(comp.pct)] <- 0
  out <- data.table(type=type, measure=measure, comp.pct=comp.pct, removed.pct=removed.pct)
  if (!is.null(group)) out$Group <- group
  return(out)
}
