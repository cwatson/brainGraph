#' Calculate the participation coefficient of each node
#'
#' This function calculates the participation coefficient of each node in a
#' graph. The graph must be split into modules beforehand. The coefficient
#' should not exceed 1 - (1 / # modules).
#'
#' @param adj.mat the adjacency matrix
#' @param adj.graph the adjacency graph
#' @param community.mem the community membership indices of each node
#' @export
#'
#' @return A vector of the participation coeff's for each node of the graph.
#' @references Guimera, R. and Amaral, L.A.N. (2005) Cartography of complex
#' networks: modules and universal roles, Journal of Statistical Mechanics:
#' Theory and Experiment, 02, P02001.

part.coeff <- function(adj.mat, adj.graph, community.mem) {
  Nv <- vcount(adj.graph)
  Nc <- max(community.mem)
  vertices <- unname(which(degree(adj.graph) != 0))

  PC <- vector(length=Nv)
  PC2 <- vector('list', length=length(vertices))

  PC2 <- foreach (i=1:length(vertices)) %dopar% {
    Kis <- vapply(1:Nc, function(x)
                  length(E(adj.graph)[vertices[i] %--% which(community.mem==x)]),
                  integer(1))
    Ki <- degree(adj.graph)[vertices[i]]
    1 - sum((Kis/Ki)^2)
  }

  PC[vertices] <- unlist(PC2)
  PC <- ifelse(is.nan(PC), 0, PC)
  names(PC) <- V(adj.graph)$name
  return(PC)

}
