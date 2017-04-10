#' Calculate graph global, local, or nodal efficiency
#'
#' This function calculates the global efficiency of a graph or the local or
#' nodal efficiency of each vertex of a graph. The global efficiency is equal
#' to the mean of all nodal efficiencies.
#'
#' Global efficiency for graph \emph{G} with \emph{N} vertices is:
#' \deqn{E_{global}(G) = \frac{1}{N(N-1)} \sum_{i \ne j \in G} \frac{1}{d_{ij}}}
#' where \eqn{d_{ij}} is the shortest path length between vertices \emph{i} and
#' \emph{j}.
#'
#' Local efficiency for vertex \emph{i} is:
#' \deqn{E_{local}(i) = \frac{1}{N} \sum_{i \in G} E_{global}(G_i)}
#' where \eqn{G_i} is the subgraph of neighbors of \emph{i}, and \emph{N} is the
#' number of vertices in that subgraph.
#'
#' Nodal efficiency for vertex \emph{i} is:
#' \deqn{E_{nodal}(i) = \frac{1}{N-1} \sum_{j \in G} \frac{1}{d_{ij}}}
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
#'
#' @return A numeric vector of the local efficiencies for each vertex of the
#'   graph (if \emph{type} is \code{local|nodal}) or a single number (if
#'   \emph{type} is \code{global}).
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Latora V., Marchiori M. (2001) \emph{Efficient behavior of
#' small-world networks}. Phys Rev Lett, 87.19:198701.

efficiency <- function(g, type=c('local', 'nodal', 'global'), weights=NULL,
                       use.parallel=TRUE, A=NULL) {
  stopifnot(is_igraph(g))
  i <- NULL
  if (is.null(weights) && 'weight' %in% edge_attr_names(g)) {
    weights <- NULL
  } else {
    if (!is.null(weights) && any(!is.na(weights))) {
      weights <- as.numeric(weights)
    } else {
      weights <- NA
    }
  }

  type <- match.arg(type)
  if (type == 'local') {
    if ('degree' %in% vertex_attr_names(g)) {
      degs <- V(g)$degree
    } else {
      degs <- degree(g)
    }

    if (is.null(weights)) {
      if (is.null(A)) A <- as_adj(g, names=FALSE, attr='weight')
      weighted <- TRUE
    } else {
      A <- as_adj(g, names=FALSE, sparse=FALSE)
      weighted <- NULL
    }
    eff <- rep(0, length(degs))
    nodes <- which(degs > 1)
    X <- apply(A, 1, function(x) which(x > 0))

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
    Nv <- vcount(g)
    eff <- apply(distances(g, weights=weights), 2, function(x)
                 sum(1 / x[x != 0]) / (Nv - 1))
    if (type == 'global') eff <- sum(eff) / length(eff)
  }
  return(eff)
}

#' @inheritParams efficiency
#' @export
#' @rdname efficiency

graph.efficiency <- function(g, type=c('local', 'nodal', 'global'), weights=NULL,
                             use.parallel=TRUE, A=NULL) {
  .Deprecated('efficiency')
  efficiency(g, type=type, weights=weights, use.parallel=use.parallel, A=A)
}
