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
#' \item{phi.orig}{A vector of the rich-club coefficients for the original
#' graph.}
#' \item{phi.norm}{A named vector of the normalized rich club coefficients.}
#' \item{p}{The p-value based on the \emph{N} random graphs generated.}
#'
#' @seealso \code{\link{rich.club.coeff}, \link{sim.rand.graph.par}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Colizza V., Flammini A., Serrano M.A., Vespignani A. (2006)
#' \emph{Detecting rich-club ordering in complex networks}. Nature Physics,
#' 2:110-115.

rich.club.norm <- function(g, N=1e2, rand=NULL, ...) {
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if (is.null(rand)) {
    rand <- sim.rand.graph.par(g, N, clustering=F)
  } else {
    if (!all(vapply(rand, is_igraph, logical(1)))) {
      stop('Argument "rand" must be a list of igraph graph objects!')
    }
    N <- length(rand)
  }

  max.deg <- max(V(g)$degree)
  if (all(vapply(rand, function(x) 'rich' %in% graph_attr_names(x), logical(1)))) {
    phi.rand <- aaply(rand, .margins=1, function(x) x$rich$phi)
  } else {
    phi.rand <- laply(rand, function(x)
                            vapply(seq_len(max.deg), function(y)
                                   rich.club.coeff(x, y, ...)$phi, numeric(1)),
                            .parallel=T)
  }
  phi.orig <- g$rich$phi
  phi.norm <- phi.orig / colMeans(phi.rand)
  p <- vapply(seq_len(max.deg), function(x)
              sum(phi.rand[, x] >= phi.orig[x]) / N,
              numeric(1))
  return(list(phi.rand=phi.rand, phi.orig=phi.orig, phi.norm=phi.norm, p=p))
}
