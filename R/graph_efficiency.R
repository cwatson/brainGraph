#' Calculate graph global, local, or nodal efficiency
#'
#' This function calculates the global efficiency of a graph or the local or
#' nodal efficiency of each vertex of a graph.
#'
#' Local efficiency for vertex \emph{i} is:
#' \deqn{E_{local}(i) = \frac{1}{N} \sum_{i \in G} E_{global}(G_i)}
#' where \eqn{G_i} is the subgraph of neighbors of \emph{i}, and \emph{N} is the
#' number of vertices in that subgraph.
#'
#' Nodal efficiency for vertex \emph{i} is:
#' \deqn{E_{nodal}(i) = \frac{1}{N-1} \sum_{j \in G} \frac{1}{d_{ij}}}
#'
#' Global efficiency for graph \emph{G} with \emph{N} vertices is:
#' \deqn{E_{global}(G) = \frac{1}{N(N-1)} \sum_{i \ne j \in G} \frac{1}{d_{ij}}}
#' where \eqn{d_{ij}} is the shortest path length between vertices \emph{i} and
#' \emph{j}. Alternatively, global efficiency is equal to the mean of all nodal
#' efficiencies.
#'
#' @param g An \code{igraph} graph object
#' @param type Character string; either \code{local}, \code{nodal}, or
#'   \code{global} (default: \code{local})
#' @param weights Numeric vector of edge weights; if \code{NULL} (the default),
#'   and if the graph has edge attribute \code{weight}, then that will be used.
#'   To avoid using weights, this should be \code{NA}.
#' @param use.parallel Logical indicating whether or not to use \code{foreach}
#'   (default: \code{TRUE})
#' @param A Numeric matrix; the (weighted or unweighted) adjacency matrix of the
#'   input graph (default: \code{NULL})
#' @export
#' @importFrom Matrix rowSums
#'
#' @return A numeric vector of the efficiencies for each vertex of the graph
#'   (if \emph{type} is \code{local|nodal}) or a single number (if \emph{type}
#'   is \code{global}).
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Latora, V. and Marchiori, M. (2001) Efficient behavior of
#'   small-world networks. \emph{Phys Rev Lett}, \bold{87.19}, 198701.
#'   \url{https://dx.doi.org/10.1103/PhysRevLett.87.198701}
#' @references Latora, V. and Marchiori, M. (2003) Economic small-world
#'   behavior in weighted networks. \emph{Eur Phys J B}, \bold{32}, 249--263.
#'   \url{https://dx.doi.org/10.1140/epjb/e2003-00095-5}

efficiency <- function(g, type=c('local', 'nodal', 'global'), weights=NULL,
                       use.parallel=TRUE, A=NULL) {
  stopifnot(is_igraph(g))
  i <- NULL
  weights <- check_weights(g, weights)

  type <- match.arg(type)
  if (type == 'local') {
    if (is.null(weights)) {
      if (is.null(A)) A <- as_adj(g, names=FALSE, attr='weight')
      weighted <- TRUE
    } else {
      A <- as_adj(g, names=FALSE, sparse=FALSE)
      weighted <- NULL
    }
    eff <- rep(0, nrow(A))
    nodes <- which(rowSums((A > 0) + 0) > 1)
    X <- apply(A, 1, function(x) which(x > 0))
    if (is.matrix(X)) X <- as.list(data.frame(X))   # If the graph is complete

    if (length(nodes) > 0) {
      if (isTRUE(use.parallel)) {
        eff[nodes] <- foreach (i=nodes, .combine='c') %dopar% {
          g.sub <- graph_from_adjacency_matrix(A[X[[i]], X[[i]]], mode='undirected', weighted=weighted)
          efficiency(g.sub, 'global', weights=weights)
        }
      } else {
        for (i in nodes) {
          g.sub <- graph_from_adjacency_matrix(A[X[[i]], X[[i]]], mode='undirected', weighted=weighted)
          eff[i] <- efficiency(g.sub, 'global', weights=weights)
        }
      }
    }
  } else {
    D <- distances(g, weights=weights)
    Nv <- nrow(D)
    Dinv <- 1 / D
    eff <- colSums(Dinv * is.finite(Dinv), na.rm=T) / (Nv - 1)
    if (type == 'global') eff <- sum(eff) / length(eff)
  }
  return(eff)
}
