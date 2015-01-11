#' Calculate the participation coefficient of each node
#'
#' This function calculates the participation coefficient of each node in a
#' graph. The graph must be split into modules beforehand. The coefficient
#' should not exceed 1 - (1 / # modules).
#'
#' @param g The graph
#' @param community.mem The community membership indices of each node
#' @export
#'
#' @return A vector of the participation coeff's for each node of the graph.
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Guimera, R. and Amaral, L.A.N. (2005) Cartography of complex
#' networks: modules and universal roles, Journal of Statistical Mechanics:
#' Theory and Experiment, 02, P02001.

part.coeff <- function(g, community.mem) {
  if ('degree' %in% vertex_attr_names(g)) {
    degrees <- V(g)$degree
  } else {
    degrees <- degree(g)
  }
  Nv <- length(degrees)
  Nc <- max(community.mem)
  vertices <- unname(which(degrees != 0))

  PC <- vector(length=Nv)
  PC2 <- vector('list', length=length(vertices))

  edges <- E(g)
  PC2 <- foreach (i=seq_along(vertices)) %dopar% {
    Kis <- vapply(seq_len(Nc), function(x)
                  length(edges[vertices[i] %--% which(community.mem==x)]),
                  integer(1))
    Ki <- degrees[vertices[i]]
    1 - sum((Kis/Ki)^2)
  }

  PC[vertices] <- unlist(PC2)
  PC <- ifelse(is.nan(PC), 0, PC)
  names(PC) <- V(g)$name
  return(PC)
}
