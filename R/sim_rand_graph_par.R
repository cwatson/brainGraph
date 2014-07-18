#' Simulate N random graphs w/ same clustering and degree sequence as the input.
#'
#' This function will simulate N simple random graphs with the same clustering
#' and degree sequence as the input. Essentially a wrapper for
#' \code{\link{sim.rand.graph.clust}} and
#' \code{\link{set.brainGraph.attributes}}. It uses \code{\link{foreach}} to
#' speed it up. If you do not want to match by clustering, then it will do a
#' simple rewiring of the given graph (1e4 times).
#'
#' @param g A graph with the characteristics for simulation of random graphs
#' @param N The number of iterations
#' @param clustering Logical for whether or not clustering should be controlled
#' for
#' @export
#'
#' @return A random graph with vertex and graph attributes.
#'
#' @seealso \code{\link{sim.rand.graph.clust}, \link{rewire}}

sim.rand.graph.par <- function(g, N, clustering=TRUE) {
  if (clustering == TRUE) {
    r <- foreach(i=1:N,
                  .packages=c('igraph', 'brainGraph')) %dopar% {
      tmp <- sim.rand.graph.clust(g, V(g)$degree, g$Cp)
      tmp <- set.brainGraph.attributes(tmp, rand=TRUE)
      tmp
    }
  } else {
    r <- foreach(i=1:N,
                 .packages=c('igraph', 'brainGraph')) %dopar% {
      tmp <- rewire(g, 'simple', 1e4)
      tmp <- set.brainGraph.attributes(tmp, rand=TRUE)
      tmp
    }
  }
}
