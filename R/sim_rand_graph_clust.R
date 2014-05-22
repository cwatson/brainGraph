#' Simulate a random graph with given degree sequence and clustering.
#'
#' This function will simulate a random graph with a given degree sequence and
#' clustering coefficient. This function calls \link{choose.edges} to decide
#' which edges will be re-wired.
#'
#' @param d The degree sequence
#' @param c The clustering measure (clustering coeff, transitivity, etc.)
#' @export
#'
#' @return A random graph

sim.rand.graph.clust <- function(d, c) {
  g <- degree.sequence.game(d, method='simple.no.multiple')

  while (transitivity(g) < c) {

    repeat {
      g.cand <- g
  
      # If E(y1, y2) and E(z1, z2) don't exist, rewire 2 edges
      repeat {
        e <- choose.edges(g.cand)
        if ( (length(E(g.cand)[e$y1 %--% e$y2]) == 0) &&
             (length(E(g.cand)[e$z1 %--% e$z2]) == 0)) {
          break
        }
      }

      g.cand <- g.cand - E(g.cand, P=c(e$y1, e$z1))
      g.cand <- g.cand - E(g.cand, P=c(e$y2, e$z2))
      g.cand <- g.cand + edges(c(e$y1, e$y2))
      g.cand <- g.cand + edges(c(e$z1, e$z2))
    
      if (transitivity(g.cand) > transitivity(g)) {
        break
      }
    }

    g <- g.cand
  }

  return(g)
}
