#' Calculate the normalized rich club coefficient
#'
#' This function will generate a number of random graphs, calculate their rich
#' club coefficients (\eqn{\phi}), and return \eqn{\phi} of the graph of
#' interest divided by the mean across random graphs, i.e. \eqn{\phi_{norm}}.
#' If random graphs have already been generated, you can supply a list as an
#' argument (since graph generation is time consuming).
#'
#' @param g The igraph graph object of interest
#' @param N The number of random graphs to generate (default: 100)
#' @param rand A list of igraph graph objects, if random graphs have already
#' been generated (default: NULL)
#' @param ... Other parameters (passed to \emph{rich.club.coeff})
#' @export
#'
#' @return A list with two elements:
#' \item{phi.rand}{A matrix with \emph{N} rows and \emph{max(degree(g))}
#' columns, where each row contains the coefficients for a given random graph.}
#' \item{phi.norm}{A named vector of the normalized rich club coefficients.}
#'
#' @seealso \code{\link{rich.club.coeff}, \link{sim.rand.graph.par}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Colizza V., Flammini A., Serrano M.A., Vespignani A. (2006)
#' \emph{Detecting rich-club ordering in complex networks}. Nature Physics,
#' 2:110-115.

rich.club.norm <- function(g, N=1e2, rand=NULL, ...) {
  if (is.null(rand)) {
    rand <- sim.rand.graph.par(g, N, clustering=F)
  } else {
    if (!all(sapply(rand, is_igraph))) {
      stop('Argument "rand" must be a list of igraph graph objects!')
    }
  }

  phi <- plyr::laply(rand, function(x)
                     sapply(seq_len(max(V(g)$degree)), function(y)
                            rich.club.coeff(x, y, ...)$phi),
                     .parallel=T)
  phi.norm <- g$rich$phi / colMeans(phi)
  return(list(phi.rand=phi, phi.norm=phi.norm))
}
