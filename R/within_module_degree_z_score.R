#' Calculate the within-module degree z-score of each node.
#'
#' This function calculates the within-module degree z-score of each node in a
#' graph. The graph must be split into modules beforehand. This is a measure of
#' the connectivity from a given node to other nodes in its module. See Guimera
#' et al., J Stat Mech, 2005.
#'
#' @param g The graph
#' @param community.mem The community membership indices of each node
#' @export
#'
#' @return A vector of the within-module degree z-scores for each node of the
#' graph.
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Guimera, R. and Amaral, L.A.N. (2005) Cartography of complex
#' networks: modules and universal roles, Journal of Statistical Mechanics:
#' Theory and Experiment, 02, P02001.

within.module.deg.z.score <- function(g, community.mem) {
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  Nv <- length(degs)
  vs <- unname(which(degs != 0))

  zs <- vector(length=Nv)
  zs2 <- vector('list', length=length(vs))

  edges <- E(g)
  zs2 <- foreach (i=seq_along(vs)) %dopar% {
    si <- which(community.mem==community.mem[vs[i]])
    K <- vapply(si, function(x) length(edges[x %--% si]), numeric(1))
    Ki <- length(edges[vs[i] %--% si])
    (Ki - mean(K)) / sd(K)
  }

  zs[vs] <- unlist(zs2)
  zs <- ifelse(is.na(zs), 0, zs)
  zs <- ifelse(is.infinite(zs), 0, zs)
  names(zs) <- V(g)$name
  return(zs)
}
