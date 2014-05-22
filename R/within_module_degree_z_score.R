#' Calculate the within-module degree z-score of each node.
#'
#' This function calculates the within-module degree z-score of each node in a
#' graph. The graph must be split into modules beforehand. This is a measure of
#' the connectivity from a given node to other nodes in its module. See Guimera
#' et al., J Stat Mech, 2005.
#'
#' @param adj.mat the adjacency matrix
#' @param adj.graph the adjacency graph
#' @param community.mem the community membership indices of each node
#' @export
#'
#' @return A vector of the within-module degree z-scores for each node of the
#' graph.
#'
#' @references Guimera, R. and Amaral, L.A.N. (2005) Cartography of complex
#' networks: modules and universal roles, Journal of Statistical Mechanics:
#' Theory and Experiment, 02, P02001.

within.module.deg.z.score <- function(adj.mat, adj.graph, community.mem) {
  Nv <- vcount(adj.graph)
  vertices <- unname(which(degree(adj.graph) != 0))

  zs <- vector(length=Nv)
  zs2 <- vector('list', length=length(vertices))

  zs2 <- foreach (i=1:length(vertices)) %dopar% {
    si <- which(community.mem==community.mem[vertices[i]])
    K <- vapply(si, function(x) length(E(adj.graph)[x %--% si]), numeric(1))
    Ki <- length(E(adj.graph)[vertices[i] %--% si])
    (Ki - mean(K)) / sd(K)
  }

  zs[vertices] <- unlist(zs2)
  zs <- ifelse(is.na(zs), 0, zs)
  zs <- ifelse(is.infinite(zs), 0, zs)
  names(zs) <- V(adj.graph)$name
  return(zs)
}
