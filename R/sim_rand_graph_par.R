#' Simulate N random graphs w/ same clustering and degree sequence as the input.
#'
#' Simulates \code{N} simple random graphs with the same clustering (optional)
#' and degree sequence as the input. Essentially a wrapper for
#' \code{\link[igraph]{sample_degseq}} (or \code{\link{sim.rand.graph.clust}})
#' and \code{\link{set_brainGraph_attr}}. It uses \code{\link[foreach]{foreach}}
#' for parallel processing.
#'
#' If you do not want to match by clustering, then simple rewiring of the input
#' graph is performed (the number of rewire's equaling the larger of \code{1e4}
#' and \eqn{10 \times m}, where \eqn{m} is the graph's edge count).
#'
#' @param g An \code{igraph} graph object
#' @param N Integer; the number of random graphs to simulate (default: 100)
#' @param clustering Logical; whether or not to control for clustering (default:
#'   \code{FALSE})
#' @param ... Other parameters (passed to \code{\link{sim.rand.graph.clust}})
#' @export
#'
#' @return A \emph{list} of \emph{N} random graphs with vertex and graph
#'   attributes
#'
#' @family Null graph functions
#' @seealso \code{\link[igraph]{rewire}, \link[igraph]{sample_degseq},
#'   \link[igraph]{keeping_degseq}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' rand1 <- sim.rand.graph.par(g[[1]][[N]], N=1e3)
#' rand1.cl <- sim.rand.graph.par(g[[1]][[N]], N=1e2,
#'   clustering=T, max.iters=1e3)
#' }

sim.rand.graph.par <- function(g, N=100, clustering=FALSE, ...) {
  stopifnot(is_igraph(g))
  if (isTRUE(clustering)) {
    r <- foreach(i=seq_len(N), .packages=c('igraph', 'brainGraph')) %dopar% {
      tmp <- sim.rand.graph.clust(g, ...)
      tmp <- set_brainGraph_attr(tmp, rand=TRUE)
      tmp
    }
  } else {
    if (is_connected(g)) {
      V(g)$degree <- degree(g)
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- sample_degseq(V(g)$degree, method='vl')
        tmp <- set_brainGraph_attr(tmp, rand=TRUE)
        tmp
      }
    } else {
      iters <- max(10*ecount(g), 1e4)
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- rewire(g, keeping_degseq(loops=F, iters))
        tmp <- set_brainGraph_attr(tmp, rand=TRUE)
        tmp
      }
    }
  }
}
