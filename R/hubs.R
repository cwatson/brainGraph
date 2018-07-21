#' Calculate vertex hubness
#'
#' \code{hubness} calculates the "hubness" (see reference) of the vertices in a
#' graph. These are vertices which meet at least two of the following four
#' criteria:
#' \enumerate{
#'   \item Have high degree/strength
#'   \item Have high betweenness centrality
#'   \item Have low clustering coefficient
#'   \item Have low average path length
#' }
#' For each criterion, "high" or "low" means "in the top 20\%" across all
#' vertices. Vertices meeting any of the criteria get a value of 1 for that
#' metric; these are summed to yield the hubness score which ranges from 0-4. As
#' in the reference article, vertices with a score of 2 or higher are to be
#' considered hubs, although that determination isn't made in this function.
#'
#' @inheritParams efficiency
#' @inheritParams xfm.weights
#' @export
#'
#' @return A numeric vector with the vertices' hubness score
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references van den Heuvel M.P., Mandl R.C.W., Stam C.J., Kahn R.S., Pol
#'   H.E.H. (2010) \emph{Aberrant frontal and temporal complex network structure
#'   in schizophrenia: a graph theoretical analysis}. The Journal of
#'   Neuroscience, 30(47):15915-15926.

hubness <- function(g, xfm.type=g$xfm.type, weights=NULL) {
  stopifnot(is_igraph(g))
  Nv <- vcount(g)
  cutoff <- round(0.2 * Nv)
  weights <- check_weights(g, weights)

  S <- strength(g, weights=weights)
  Cp <- transitivity(g, type='weighted', isolates='zero', weights=weights)
  btwn <- centr_betw(g)$res

  if (is_weighted(g)) g <- xfm.weights(g, xfm.type)
  Lpv <- distances(g, weights=weights)
  Lpv[is.infinite(Lpv)] <- NA
  Lp <- rowMeans(Lpv, na.rm=TRUE)

  M <- matrix(c(S, btwn, Cp, Lp), nrow=Nv, ncol=4, byrow=FALSE)
  H <- matrix(0, nrow=Nv, ncol=4)
  for (i in 1:2) {
    H[order(M[, i], decreasing=TRUE)[1:cutoff], i] <- 1
  }
  for (i in 3:4) {
    H[order(M[, i])[1:cutoff], i] <- 1
  }
  hubs <- rowSums(H)
  return(hubs)
}
