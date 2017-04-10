#' Simulate a random graph with given degree sequence and clustering.
#'
#' This function will simulate a random graph with a given degree sequence and
#' clustering coefficient. This function calls \link{choose.edges} to decide
#' which edges will be re-wired.
#'
#' @param g An \code{igraph} graph object from which null graphs are simulated
#' @param cl The clustering measure (default: transitivity)
#' @param max.iters The maximum number of iterations to perform (default: 100)
#' @export
#'
#' @return An \code{igraph} graph object
#' @seealso \code{\link{choose.edges}, \link[igraph]{rewire},
#' \link[igraph]{transitivity}, \link[igraph]{keeping_degseq}}
#'
#' @references Bansal S., Khandelwal S., Meyers L.A. (2009) \emph{Exploring
#' biological network structure with clustered random networks}. BMC
#' Bioinformatics, 10:405-421.

sim.rand.graph.clust <- function(g, cl=g$transitivity, max.iters=100) {
  g <- rewire(g, keeping_degseq(loops=F, 1e4))

  cur.iter <- 0
  while ((transitivity(g) < cl) & (cur.iter < max.iters)) {
    repeat {
      g.cand <- g
      A <- as_adj(g.cand, sparse=FALSE, names=FALSE)
      degs <- colSums(A)
      degs.large <- which(degs > 1)

      # If E(y1, y2) and E(z1, z2) don't exist, rewire 2 edges
      repeat {
        e <- choose.edges(A, degs, degs.large)
        if ( (A[e$y1, e$y2] == 0) && (A[e$z1, e$z2] == 0) &&
             (e$y1 != e$y2) && (e$z1 != e$z2) ) {
          break
        }
      }

      A[e$y1, e$z1] <- A[e$z1, e$y1] <- A[e$y2, e$z2] <- A[e$z2, e$y2] <- 0
      A[e$y1, e$y2] <- A[e$y2, e$y1] <- A[e$z1, e$z2] <- A[e$z2, e$z1] <- 1
      g.cand <- graph_from_adjacency_matrix(A, mode='undirected')

      if (transitivity(g.cand) > transitivity(g)) break
    }

    cur.iter <- cur.iter + 1
    g <- g.cand
  }

  return(g)
}
