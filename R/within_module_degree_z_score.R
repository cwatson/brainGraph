#' Calculate vertex within-module degree z-score
#'
#' This function calculates the within-module degree z-score of each vertex in a
#' graph, based on some module membership. This is a measure of the connectivity
#' from a given vertex to other vertices in its module.
#'
#' The within-module degree z-score is:
#' \deqn{z_i = \frac{\kappa_i - \bar{\kappa}_{s_i}}{\sigma_{\kappa_{s_i}}}}
#' where \eqn{\kappa_i} is the number of edges from vertex \emph{i} to vertices
#' in the same module \eqn{s_i}, \eqn{\bar{\kappa}_{s_i}} is the average of
#' \eqn{\kappa} over all vertices in \eqn{s_i}, and \eqn{\sigma_{\kappa_{s_i}}}
#' is the standard deviation.
#'
#' @param g The graph
#' @param memb The community membership indices of each vertex
#' @export
#'
#' @return A vector of the within-module degree z-scores for each vertex of the
#' graph.
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Guimera, R. and Amaral, L.A.N. (2005) Cartography of complex
#' networks: modules and universal roles, Journal of Statistical Mechanics:
#' Theory and Experiment, 02, P02001.

within_module_deg_z_score <- function(g, memb) {
  i <- NULL
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  vs <- which(degs > 0)
  es <- E(g)
  z <- Ki <- rep(0, length(degs))

  Ki[vs] <- foreach(i=vs, .combine='c') %dopar% {
    length(es[i %--% which(memb == memb[i])])
  }
  di <- lapply(seq_len(max(memb)), function(x) Ki[memb == x])
  Ksi <- vapply(di, mean, numeric(1))
  sigKsi <- vapply(di, sd, numeric(1))

  z[vs] <- (Ki[vs] - Ksi[memb[vs]]) / sigKsi[memb[vs]]
  z <- ifelse(!is.finite(z), 0, z)
  return(z)
}
