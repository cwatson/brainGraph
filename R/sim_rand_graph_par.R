#' Simulate N random graphs w/ same clustering and degree sequence as the input.
#'
#' This function will simulate N simple random graphs with the same clustering
#' and degree sequence as the input. Essentially a wrapper for
#' \code{\link{sim.rand.graph.clust}} and
#' \code{\link{set.brainGraph.attributes}}. It uses
#' \code{\link[foreach]{foreach}} to speed it up. If you do not want to match by
#' clustering, then it will do a simple rewiring of the given graph (the number
#' of rewire's equaling the larger of 1e4 and 10 * number of edges).
#'
#' @param g A graph with the characteristics for simulation of random graphs
#' @param N The number of iterations
#' @param clustering Logical for whether or not to control for clustering
#' @param ... Other parameters (passed to \code{\link{sim.rand.graph.clust}})
#' @export
#'
#' @return A list of \emph{N} random graphs with vertex and graph attributes.
#'
#' @seealso \code{\link{sim.rand.graph.clust}, \link[igraph]{rewire}}
#' @examples
#' \dontrun{
#' rand1 <- sim.rand.graph.par(g1[[N]], N=1e3, clustering=F)
#' rand1.cl <- sim.rand.graph.par(g1[[N]], N=1e2, max.iters=1e3)
#' }

sim.rand.graph.par <- function(g, N, clustering=TRUE, ...) {
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if (clustering == TRUE) {
    r <- foreach(i=seq_len(N),
                  .packages=c('igraph', 'brainGraph')) %dopar% {
      tmp <- sim.rand.graph.clust(g, ...)
      tmp$g <- set.brainGraph.attributes(tmp$g, rand=TRUE)
      tmp
    }
  } else {
    if (is_connected(g)) {
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- sample_degseq(V(g)$degree, method='vl')
        tmp <- set.brainGraph.attributes(tmp, rand=TRUE)
        tmp
      }
    } else {
      iters <- max(10*ecount(g), 1e4)
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- rewire(g, keeping_degseq(loops=F, iters))
        tmp <- set.brainGraph.attributes(tmp, rand=TRUE)
        tmp
      }
    }
  }
}
