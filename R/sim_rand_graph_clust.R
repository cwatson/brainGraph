#' Simulate a random graph with given degree sequence and clustering.
#'
#' This function will simulate a random graph with a given degree sequence and
#' clustering coefficient. This function calls \link{choose.edges} to decide
#' which edges will be re-wired.
#'
#' @param graph The graph from which null graphs are simulated
#' @param d The degree sequence
#' @param cl The clustering measure (clustering coeff, transitivity, etc.)
#' @export
#'
#' @return A random graph
#'
#' @seealso \code{\link{degree.sequence.game}, \link{choose.edges},
#' \link{rewire}, \link{transitivity}, \link{keeping_degseq}}

sim.rand.graph.clust <- function(graph, d, cl) {
  #g <- degree.sequence.game(d, method='simple.no.multiple')
  g <- rewire(graph, keeping_degseq(loops=F, 1e4))

  #g.all <- vector('list')

  while (transitivity(g) < cl) {
    m <- 1
    repeat {
      g.cand <- g
  
      n <- 1
      # If E(y1, y2) and E(z1, z2) don't exist, rewire 2 edges
      repeat {
        n <- n + 1
        if (n >= 100) {
          g.cand <- rewire(graph, keeping_degseq(loops=F, 1e4))
          n <- 1
        }
        e <- choose.edges(g.cand)
        if ( (!are.connected(g.cand, e$y1, e$y2)) &&
             (!are.connected(g.cand, e$z1, e$z2)) &&
             (e$y1 != e$y2) && (e$z1 != e$z2) ) {
          break
        }
      }

      g.cand <- g.cand - E(g.cand, P=c(e$y1, e$z1))
      g.cand <- g.cand - E(g.cand, P=c(e$y2, e$z2))
      g.cand <- g.cand + edges(c(e$y1, e$y2)) + edges(c(e$z1, e$z2))

      #g.all[[length(g.all)+1]] <- g.cand

      m <- m + 1
      if (m %% 100 == 0) {
        g.cand <- rewire(graph, keeping_degseq(loops=F, 1e4))
        m <- 1
        n <- 1
      }
      if (transitivity(g.cand) > transitivity(g)) {
        break
      }
    }

    g <- g.cand
  }

  return(g)#, g.all))
}
