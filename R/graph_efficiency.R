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
#' @param g The graph on which to calculate efficiency
#' @param type A character string; either 'local', 'nodal', or 'global'
#' @param weights A numeric vector of edge weights; if 'NULL', and if the graph
#' has edge attribute 'weight', then that will be used. To avoid using weights,
#' this should be 'NA'
#' @export
#'
#' @return A vector of the local efficiencies for each node of the graph (if 
#' type='local|nodal') or a number (if type='global').
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Latora V., Marchiori M. (2001) \emph{Efficient behavior of
#' small-world networks}. Phys Rev Lett, 87.19:198701.

graph.efficiency <- function(g, type=c('local', 'nodal', 'global'),
                             weights=NULL) {
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }

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
    eff <- numeric(length(degs))
    nodes <- which(degs > 1)
  
    eff[nodes] <- simplify2array(mclapply(nodes, function(x) {
      neighbs <- neighbors(g, v=x)
      g.sub <- induced.subgraph(g, neighbs)
      Nv <- vcount(g.sub)
  
      paths <- shortest.paths(g.sub, weights=weights)
      paths <- paths[upper.tri(paths)]
      2 / Nv / (Nv - 1) * sum(1 / paths[paths != 0])
      }, mc.cores=detectCores())
    )
  } else if (type == 'nodal') {
    Nv <- length(degs)
    paths <- shortest.paths(g, weights=weights)
    eff <- simplify2array(mclapply(seq_len(Nv), function(x) {
      path <- paths[x, ]
      eff <- sum(1 / path[path != 0]) / (Nv - 1)
      }, mc.cores=detectCores())
    )
  } else {
    Nv <- length(degs)
    nodes <- V(g)[degs != 0]
    paths <- shortest.paths(g, v=nodes, to=nodes, weights=weights)
    paths <- paths[upper.tri(paths)]

    eff <- 2 / Nv / (Nv - 1) * sum(1 / paths[paths != 0])
  }
  return(eff)
}
