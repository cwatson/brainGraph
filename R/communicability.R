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
#'
#' @return A numeric matrix of the communicability
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Estrada, E. and Hatano, N. (2008) Communicability in complex
#'   networks. \emph{Physical Review E}. \bold{77}, 036111.
#'   \doi{10.1103/PhysRevE.77.036111}
#' @references Crofts, J.J. and Higham, D.J. (2009) A weighted communicability
#'   measure applied to complex brain networks. \emph{J. R. Soc. Interface}.
#'   \bold{6}, 411--414. \doi{10.1098/rsif.2008.0484}

communicability <- function(g, weights=NULL) {
  if (!requireNamespace('expm', quietly=TRUE)) stop('Must install the "expm" package.')
  stopifnot(is_igraph(g), is_connected(g))
  if (is.null(weights) && 'weight' %in% edge_attr_names(g)) {
    d <- check_strength(g)
    S <- diag(1 / sqrt(d))
    A <- as.matrix(as_adj(g, names=FALSE, attr='weight'))
  } else {
    S <- diag(vcount(g))
    A <- as_adj(g, names=FALSE, sparse=FALSE)
  }

  C <- expm::expm(S %*% A %*% S)
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
#' @inheritParams efficiency
#' @export
#'
#' @return A numeric vector of the centrality for each vertex
#'
#' @family Centrality functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Estrada, E. and Higham, D.J. and Hatano N. (2009) Communicability
#'   betweenness in complex networks. \emph{Physica A}, \bold{388}, 764--774.
#'   \doi{10.1016/j.physa.2008.11.011}

centr_betw_comm <- function(g, A=NULL) {
  stopifnot(is_igraph(g), is_connected(g))
  if (!requireNamespace('expm', quietly=TRUE)) stop('Must install the "expm" package.')
  if (is.null(A)) A <- as_adj(g, names=FALSE, sparse=FALSE)

  C <- expm::expm(A)
  N <- dim(A)[1L]
  n <- (N - 1)^2 - (N - 1)
  Wr <- rep.int(0, N)
  for (i in seq_len(N)) {
    Er <- A
    Er[i, ] <- Er[, i] <- 0
    G <- (C - expm::expm(Er)) / C
    Wr[i] <- (sum(G[-i, -i]) - sum(diag(G[-i, -i]))) / n
  }
  return(Wr)
}
