#' Calculate vertex hubness
#'
#' \code{hubness} calculates the \dQuote{hubness} (see reference) of the
#' vertices in a graph. These are vertices which meet at least two of the
#' following four criteria:
#' \enumerate{
#'   \item Have high degree/strength
#'   \item Have high betweenness centrality
#'   \item Have low clustering coefficient
#'   \item Have low average path length
#' }
#' For each criterion, \dQuote{high} or \dQuote{low} means \dQuote{in the top
#' 20\%} across all vertices. Vertices meeting any of the criteria get a value
#' of 1 for that metric; these are summed to yield the hubness score which
#' ranges from 0-4. As in the reference article, vertices with a score of 2 or
#' higher are to be considered hubs, although that determination isn't made in
#' this function.
#'
#' @param prop.keep Numeric (between 0 and 1) indicating the proportion of
#'   vertices to consider as having a high score. Default: 0.2 (20\%)
#' @inheritParams efficiency
#' @inheritParams xfm.weights
#' @export
#'
#' @return A numeric vector with the vertices' hubness score
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references van den Heuvel, M.P. and Mandl, R.C.W. and Stam, C.J. and Kahn,
#'   R.S. and Pol, H.E.H. (2010) Aberrant frontal and temporal complex network
#'   structure in schizophrenia: a graph theoretical analysis. \emph{The Journal
#'   of Neuroscience}, \bold{30(47)}, 15915--15926.
#'   \url{https://dx.doi.org/10.1523/JNEUROSCI.2874-10.2010}

hubness <- function(g, xfm.type=g$xfm.type, weights=NULL, prop.keep=0.2) {
  stopifnot(is_igraph(g), prop.keep <= 1, prop.keep >= 0)
  Nv <- vcount(g)
  cutoff <- round(prop.keep * Nv)
  weights <- check_weights(g, weights)

  S <- strength(g, weights=weights)
  Cp <- transitivity(g, type='weighted', isolates='zero', weights=weights)
  btwn <- centr_betw(g)$res

  vattr <- 'Lp'
  if (is_weighted(g)) {
    g <- xfm.weights(g, xfm.type)
    vattr <- 'Lp.wt'
  }
  if (vattr %in% vertex_attr_names(g)) {
    Lp <- vertex_attr(g, vattr)
  } else {
    Lp <- mean_distance_wt(g, 'vertex', weights=weights)
  }

  M <- matrix(c(-S, -btwn, Cp, Lp), nrow=Nv, ncol=4L)
  H <- matrix(0L, nrow=Nv, ncol=4L)
  for (i in 1L:4L) H[order(M[, i])[1L:cutoff], i] <- 1L
  hubs <- rowSums(H)
  return(hubs)
}
