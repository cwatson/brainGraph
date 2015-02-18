#' Calculate vertex participation coefficient
#'
#' This function calculates the participation coefficient of each vertex in a
#' graph, based on community membership.
#'
#' The participation coefficient \eqn{P_i} of vertex \emph{i} is:
#' \deqn{P_i = 1 - \sum_{s=1}^{N_M} \left ( \frac{\kappa_{is}}{\kappa_i} \right )^2}
#' where \eqn{\kappa_{is}} is the number of edges from vertex \emph{i} to
#' vertices in module \emph{s}, and \eqn{\kappa_s} is the degree of vertex
#' \emph{i}. \eqn{N_M} equals the number of modules.
#'
#' As discussed in Guimera et al., \eqn{P_i = 0} if vertex \emph{i} is connected
#' only to vertices in the same module, and \eqn{P_i = 1} if vertex \emph{i} is
#' equally connected to all other modules.
#'
#' @param g The graph
#' @param memb The community membership indices of each vertex
#' @export
#'
#' @return A vector of the participation coeff's for each vertex of the graph.
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Guimera, R. and Amaral, L.A.N. (2005) Cartography of complex
#' networks: modules and universal roles, Journal of Statistical Mechanics:
#' Theory and Experiment, 02, P02001.

part.coeff <- function(g, memb) {
  i <- NULL
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  edges <- E(g)
  Nv <- length(degs)
  Nc <- max(memb)
  vs <- unname(which(degs > 0))

  PC <- vector(length=Nv)
  PC2 <- vector('list', length=length(vs))


  PC2 <- foreach (i=seq_along(vs)) %dopar% {
    Kis <- vapply(seq_len(Nc), function(x)
                  length(edges[vs[i] %--% which(memb==x)]),
                  integer(1))
    Ki <- degs[vs[i]]
    1 - sum((Kis/Ki)^2)
  }

  PC[vs] <- unlist(PC2)
  PC <- ifelse(is.nan(PC), 0, PC)
  names(PC) <- V(g)$name
  return(PC)
}
