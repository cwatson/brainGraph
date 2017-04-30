#' Calculate communicability
#'
#' \code{communicability} calculates the communicability of a network, a measure
#' which takes into account all possible paths (including non-shortest paths)
#' between vertex pairs.
#'
#' The communicability \eqn{G_{pq}} is a weighted sum of the number of walks
#' from vertex \emph{p} to \emph{q} and is calculated by taking the exponential
#' of the adjacency matrix \emph{A}:
#' \deqn{G_{pq} = \sum_{k=0}^{\infty} \frac{(\mathbf{A}^k)_{pq}}{k!} =
#' (e^{\mathbf{A}})_{pq}}
#' where \eqn{k} is \emph{walk} length.
#'
#' For weighted graphs with \eqn{D = diag(d_i)} a diagonal matrix of vertex
#' strength,
#' \deqn{G_{pq} = (e^{\mathbf{D}^{-1/2} \mathbf{A} \mathbf{D}^{-1/2}})_{pq}}
#'
#' @inheritParams efficiency
#' @export
#' @importFrom expm expm
#'
#' @return A numeric matrix of the communicability
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Estrada E. & Hatano N. (2008) \emph{Communicability in complex
#'   networks}. Physical Review E, 77:036111.
#' @references Crofts J.J. & Higham D.J. (2009) \emph{A weighted communicability
#'   measure applied to complex brain networks}. J. R. Soc. Interface, 6:411-414.

communicability <- function(g, weights=NULL) {
  stopifnot(is_igraph(g), is_connected(g))
  if (is.null(weights) && 'weight' %in% edge_attr_names(g)) {
    if ('strength' %in% vertex_attr_names(g)) {
      d <- V(g)$strength
    } else {
      d <- strength(g)
    }
    S <- diag(1 / sqrt(d))
    A <- as.matrix(as_adj(g, names=FALSE, attr='weight'))
  } else {
    S <- diag(vcount(g))
    A <- as_adj(g, names=FALSE, sparse=FALSE)
  }

  C <- expm(S %*% A %*% S)
  return(C)
}

#' Calculate communicability betweenness centrality
#'
#' \code{centr_betw_comm} calculates the \emph{communicability betweenness} of
#' the vertices of a graph. The centrality for vertex \code{r} is
#' \deqn{\omega_r = \frac{1}{C} \sum_p \sum_q \frac{(e^{\mathbf{A}})_{pq} -
#' (e^{\mathbf{A} + \mathbf{E}(r)})_{pq}}{(e^{\mathbf{A}})_{pq}}}
#' where \eqn{C = (n - 1)^2 - (n - 1)} is a normalization factor.
#'
#' @param g An \code{igraph} graph object
#' @param A Numeric matrix, the graph's adjacency matrix (default: \code{NULL})
#' @export
#' @importFrom expm expm
#'
#' @return A numeric vector of the centrality for each vertex
#'
#' @family Centrality functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Estrada E., Higham D.J., Hatano N. (2009) \emph{Communicability
#'   betweenness in complex networks}. Physica A, 388:764-774.

centr_betw_comm <- function(g, A=NULL) {
  stopifnot(is_igraph(g), is_connected(g))
  if (is.null(A)) A <- as_adj(g, names=FALSE, sparse=FALSE)

  C <- expm(A)
  N <- nrow(A)
  n <- (N - 1)^2 - (N - 1)
  Wr <- rep(0, N)
  for (i in seq_len(N)) {
    Er <- A
    Er[i, ] <- Er[, i] <- 0
    G <- (C - expm(Er)) / C
    Wr[i] <- (sum(G[-i, -i]) - sum(diag(G[-i, -i]))) / n
  }
  return(Wr)
}
