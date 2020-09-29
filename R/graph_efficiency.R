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
#'   \code{global}. Default: \code{local}
#' @param weights Numeric vector of edge weights; if \code{NULL} (the default),
#'   and if the graph has edge attribute \code{weight}, then that will be used.
#'   To avoid using weights, this should be \code{NA}.
#' @param xfm Logical indicating whether to transform the edge weights. Default:
#' \code{FALSE}
#' @param use.parallel Logical indicating whether or not to use \code{foreach}.
#'   Default: \code{TRUE}
#' @param A Numeric matrix; the adjacency matrix of the input graph. Default:
#'   \code{NULL}
#' @param D Numeric matrix; the graph's \dQuote{distance matrix}
#' @inheritParams xfm.weights
#' @export
#' @importFrom Matrix rowSums
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
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
                       xfm=FALSE, xfm.type=NULL, use.parallel=TRUE, A=NULL, D=NULL) {
  stopifnot(is_igraph(g))
  i <- NULL
  if (isTRUE(xfm)) {
    if (is.null(xfm.type)) {
      xfm.type <- if ('xfm.type' %in% graph_attr_names(g)) g$xfm.type else '1/w'
    }
    g <- xfm.weights(g, xfm.type)
  }
  weights <- check_weights(g, weights)

  type <- match.arg(type)
  if (type == 'local') {
    e.attr <- weighted <- NULL
    if (is.null(weights)) {
      e.attr <- 'weight'
      weighted <- TRUE
    }
    if (is.null(A)) A <- as_adj(g, names=FALSE, sparse=FALSE, attr=e.attr)
    eff <- rep.int(0, dim(A)[1L])
    verts <- which(rowSums((A > 0) + 0) > 1)
    X <- apply(A, 1L, function(x) which(x > 0))
    if (is.matrix(X)) X <- as.list(data.frame(X))   # If the graph is complete

    if (length(verts) > 0L) {
      `%d%` <- `%do%`
      if (isTRUE(use.parallel)) {
        if (!getDoParRegistered()) {
          cl <- makeCluster(getOption('bg.ncpus'))
          registerDoParallel(cl)
        }
        `%d%` <- `%dopar%`
      }
      eff[verts] <- foreach(i=verts, .combine='c') %d% {
        g.sub <- graph_from_adjacency_matrix(A[X[[i]], X[[i]]], mode='undirected', weighted=weighted)
        efficiency(g.sub, 'global', weights=weights)
      }
    }
  } else {
    Nv <- vcount(g)
    if (is.null(D)) {
      if (Nv > 650) {
        D <- foreach(i=seq_len(Nv), .combine=rbind) %dopar% {
          distances(g, v=i, weights=weights)
        }
      } else {
        D <- distances(g, weights=weights)
      }
    }
    Dinv <- 1 / D
    eff <- colSums(Dinv * is.finite(Dinv), na.rm=TRUE) / (Nv - 1L)
    if (type == 'global') eff <- sum(eff) / length(eff)
  }
  return(eff)
}
