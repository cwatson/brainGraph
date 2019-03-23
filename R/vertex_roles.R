#' Gateway coefficient, participation coefficient, and within-mod degree z-score
#'
#' \code{gateway_coeff} calculates the gateway coefficient of each vertex,
#' based on community membership.
#'
#' The gateway coefficient \eqn{G_i} of vertex \emph{i} is:
#' \deqn{G_i = 1 - \sum_{S=1}^{N_M} \left ( \frac{\kappa_{iS}}{\kappa_i} \right
#' )^2 (g_{iS})^2}
#' where \eqn{\kappa_{iS}} is the number of edges from vertex \emph{i} to
#' vertices in module \emph{S}, and \eqn{\kappa_i} is the degree of vertex
#' \emph{i}. \eqn{N_M} equals the number of modules. \eqn{g_{ii}} is a weight,
#' defined as:
#' \deqn{g_{iS} = 1 - \bar{\kappa_{iS}} \bar{c_{iS}}}
#' where
#' \deqn{\bar{\kappa_{iS}} = \frac{\kappa_{iS}}{\sum_j \kappa_{jS}}}
#' for all nodes \eqn{j} in node \eqn{i}'s module, and
#' \deqn{\bar{c_{iS}} = c_{iS} / max(c_n)}
#'
#' @param g An \code{igraph} graph object
#' @param memb A numeric vector of membership indices of each vertex
#' @param centr Character string; the type of centrality to use in calculating
#'   GC (default: \code{btwn.cent})
#' @export
#'
#' @return A vector of the participation coefficients, within-module degree
#'   z-scores, or gateway coefficients for each vertex of the graph.
#'
#' @name VertexRoles
#' @aliases gateway_coeff
#' @rdname vertex_roles
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Vargas, E.R. and Wahl, L.M. (2014) The gateway coefficient: a
#'   novel metric for identifying critical connections in modular networks.
#'   \emph{Eur Phys J B}, \bold{87}, 161--170.
#'   \url{https://dx.doi.org/10.1140/epjb/e2014-40800-7}

gateway_coeff <- function(g, memb, centr=c('btwn.cent', 'degree', 'strength')) {
  stopifnot(is_igraph(g))
  Ki <- check_degree(g)
  centr <- match.arg(centr)
  if (centr == 'btwn.cent') {
    if ('btwn.cent' %in% vertex_attr_names(g)) {
      cent <- V(g)$btwn.cent
    } else {
      cent <- centr_betw(g)$res
    }
  } else if (centr == 'degree') {
    cent <- Ki
  } else if (centr == 'strength') {
    cent <- check_strength(g)
  }
  N <- max(memb)
  if (N == 1) return(rep(0, length(memb)))
  Cn <- max(vapply(seq_len(N), function(x) sum(cent[which(memb == x)]), numeric(1)))

  A <- as_adj(g, sparse=FALSE, names=FALSE)
  Kis <- vapply(seq_len(N), function(x) colSums(A * (memb==x)), numeric(length(Ki)))
  M <- which(tabulate(memb) > 1)
  Kjs <- matrix(0, nrow=N, ncol=N)
  Kjs[M, M] <- vapply(M, function(x) colSums(Kis[which(memb==x), M, drop=FALSE]), numeric(length(M)))
  barKis <- Cis <- matrix(0, nrow=length(Ki), ncol=N)

  for (i in which(Ki > 0)) {
    barKis[i, ] <- Kis[i, ] / Kjs[, memb[i]]
    Vi <- which(A[, i] == 1)
    Cis[i, ] <- vapply(seq_len(N), function(x) sum(cent[Vi[memb[Vi] == x]]), numeric(1))
  }

  barCis <- Cis / Cn
  gis <- 1 - barKis * barCis
  G <- 1 - ((1 / Ki^2) * rowSums(Kis^2 * gis^2, na.rm=TRUE))

  return(G)
}

#' Participation coefficient
#'
#' \code{part_coeff} calculates the participation coefficient of each vertex,
#' based on community membership.
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
#' @export
#'
#' @aliases part_coeff
#' @rdname vertex_roles
#'
#' @references Guimera, R. and Amaral, L.A.N. (2005) Cartography of complex
#'   networks: modules and universal roles. \emph{Journal of Statistical
#'   Mechanics: Theory and Experiment}, \bold{02}, P02001.
#'   \url{https://dx.doi.org/10.1088/1742-5468/2005/02/P02001}

part_coeff <- function(g, memb) {
  stopifnot(is_igraph(g))
  Ki <- check_degree(g)
  N <- max(memb)
  A <- as_adj(g, sparse=FALSE, names=FALSE)
  Kis <- vapply(seq_len(N), function(x) colSums(A * (memb == x)), numeric(length(Ki)))
  PC <- 1 - ((1 / Ki^2) * rowSums(Kis^2))

  return(PC)
}

#' Calculate vertex within-module degree z-score
#'
#' \code{within_module_deg_z_score} is a measure of the connectivity from a
#' given vertex to other vertices in its module/community.
#'
#' The within-module degree z-score is:
#' \deqn{z_i = \frac{\kappa_i - \bar{\kappa}_{s_i}}{\sigma_{\kappa_{s_i}}}}
#' where \eqn{\kappa_i} is the number of edges from vertex \emph{i} to vertices
#' in the same module \eqn{s_i}, \eqn{\bar{\kappa}_{s_i}} is the average of
#' \eqn{\kappa} over all vertices in \eqn{s_i}, and \eqn{\sigma_{\kappa_{s_i}}}
#' is the standard deviation.
#'
#' @export
#'
#' @aliases within_module_deg_z_score
#' @rdname vertex_roles

within_module_deg_z_score <- function(g, memb) {
  stopifnot(is_igraph(g))
  N <- max(memb)
  A <- as_adj(g, sparse=FALSE, names=FALSE)
  z <- Ki <- rep(0, nrow(A))
  Ksi <- sigKsi <- rep(0, N)

  for (S in seq_len(N)) {
    x <- A[memb==S, ] %*% (memb==S)
    Ki[memb==S] <- x
    Ksi[S] <- mean(x)
    sigKsi[S] <- sd(x)
  }
  z <- (Ki - Ksi[memb]) / sigKsi[memb]
  z <- ifelse(!is.finite(z), 0, z)
  return(z)
}
