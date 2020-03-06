#' Calculate the s-core of a network
#'
#' Calculates the \emph{s-core} decomposition of a network. This is analogous to
#' the \emph{k-core} decomposition, but takes into account the \emph{strength}
#' of vertices (i.e., in weighted networks). If an unweighted network is
#' supplied, then the output of the function \code{\link[igraph]{coreness}} is
#' returned.
#'
#' The \emph{s-core} consists of all vertices \eqn{i} with \eqn{s_i > s}, where
#' \eqn{s} is some threshold value. The \eqn{s_0} core is the entire network,
#' and the threshold value of the \eqn{s_{n}} core is
#' \deqn{s_{n-1} = min_i s_i}
#' for all vertices \eqn{i} in the \eqn{s_{n-1}} core.
#'
#' Note that in networks with a wide distribution of vertex strengths, in which
#' there are almost as many unique values as there are vertices, then several
#' separate cores will have a single vertex. See the reference provided below.
#'
#' @param W Numeric matrix of edge weights (default: \code{NULL})
#' @inheritParams efficiency
#' @export
#'
#' @return Integer vector of the vertices' \emph{s-core} membership
#'
#' @seealso \code{\link[igraph]{coreness}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Eidsaa, M and Almaas, E. (2013) s-core network decomposition: a
#'   generalization of k-core analysis to weighted networks. \emph{Physical
#'   Review E}, \bold{88}, 062819.
#'   \url{https://dx.doi.org/10.1103/PhysRevE.88.062819}

s_core <- function(g, W=NULL) {
  stopifnot(is_igraph(g))
  if (!is_weighted(g)) return(coreness(g))

  if (is.null(W)) W <- as_adj(g, names=FALSE, sparse=FALSE, attr='weight')
  ct <- 1L
  s.core <- vector('integer', length=dim(W)[1L])
  repeat {
    str.tmp <- colSums(W)
    s.thr <- min(str.tmp[which(str.tmp > 0)])
    v.remove <- which(str.tmp <= s.thr & str.tmp > 0)
    if (length(v.remove) > 0L) {
      s.core[v.remove] <- ct
      W[v.remove, ] <- W[, v.remove] <- 0
      ct <- ct + 1L
    }
    if (sum(colSums(W) > 0) == 0L) break
  }
  return(s.core)
}
