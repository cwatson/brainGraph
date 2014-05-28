#' Simulate N random graphs w/ same clustering and degree sequence as the input.
#'
#' This function will simulate N simple random graphs with the same clustering
#' and degree sequence as the input. It will also calculate (and return) avg.
#' path length, transitivity, mean degree, and modularity. Essentially a wrapper
#' for \code{\link{sim.rand.graph.clust}} and
#' \code{\link{set.brainGraph.attributes}}. It uses \code{\link{foreach}} to
#' speed it up.
#'
#' @param g A graph with the characteristics for simulation of random graphs
#' @param N The number of iterations
#' @export
#'
#' @return A random graph with vertex and graph attributes.

sim.rand.graph.par <- function(g, N) {
  r <- foreach(i=1:N,
                .packages=c('igraph', 'brainGraph')) %dopar% {
    tmp <- sim.rand.graph.clust(V(g)$degree, g$Cp)
    tmp <- set.brainGraph.attributes(tmp, rand=TRUE)
    tmp
  }

}
