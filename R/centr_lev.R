#' Calculate a vertex's leverage centrality
#'
#' Calculates the leverage centrality of each vertex in a graph.
#'
#' The leverage centrality relates a vertex's degree with the degree of its
#' neighbors. The equation is:
#' \deqn{l_i = \frac{1}{k_i} \sum_{j \in N_i} \frac{k_i - k_j}{k_i + k_j}}
#' where \eqn{k_i} is the degree of the \eqn{i^{th}} vertex and \eqn{N_i} is the
#' set of neighbors of \emph{i}. This function replaces \emph{NaN} with
#' \emph{NA} (for functions that have the argument \emph{na.rm}).
#'
#' @inheritParams efficiency
#' @export
#'
#' @return A vector of the leverage centrality for all vertices.
#'
#' @family Centrality functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Joyce, K.E. and Laurienti P.J. and Burdette J.H. and Hayasaka S.
#' (2010) A new measure of centrality for brain networks. \emph{PLoS One},
#' \bold{5(8)}, e12200. \url{https://dx.doi.org/10.1371/journal.pone.0012200}

centr_lev <- function(g, A=NULL) {
  stopifnot(is_igraph(g))

  if (is.null(A)) A <- as_adj(g, sparse=FALSE, names=FALSE)
  k <- colSums(A)
  lev.cent <- rep.int(NA, dim(A)[1L])
  # This is a tiny bit slower for larger graphs
#  for (i in which(k > 0)) {
#    lev.cent[i] <- sum((k[i] - k[A[i, ] == 1]) / (k[i] + k[A[i, ] == 1])) / k[i]
#  }
  dd <- lapply(seq_along(k), function(x) k[A[x, ] == 1])
  lev.cent[which(k > 0)] <- vapply(which(k > 0), function(x) sum((k[x] - dd[[x]]) / (k[x] + dd[[x]])) / k[x], numeric(1L))

  lev.cent[is.nan(lev.cent)] <- NA
  return(lev.cent)
}
