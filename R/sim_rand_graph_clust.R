#' Simulate a random graph with given degree sequence and clustering.
#'
#' This function will simulate a random graph with a given degree sequence and
#' clustering coefficient. This function calls \link{choose.edges} to decide
#' which edges will be re-wired.
#'
#' @param graph The graph from which null graphs are simulated
#' @param cl The clustering measure (default: transitivity)
#' @param max.iters The maximum number of iterations to perform (default: 100)
#' @export
#'
#' @return A list with components:
#' \item{g}{The random graph that was generated}
#' \item{iters}{The total number of iterations performed}
#' \item{cl}{The clustering coefficient at each step}
#'
#' @seealso \code{\link{choose.edges}, \link[igraph]{rewire},
#' \link[igraph]{transitivity}, \link[igraph]{keeping_degseq}}

sim.rand.graph.clust <- function(graph, cl=graph$transitivity, max.iters=100) {
  g <- rewire(graph, keeping_degseq(loops=F, 1e4))

  #g.all <- vector('list')

  cur.iter <- 0
  cl.all <- loop.time <- vector(length=length(max.iters))
  while ((transitivity(g) < cl) & (cur.iter < max.iters)) {
    start.time <- Sys.time()
    repeat {
      g.cand <- g
  
      # If E(y1, y2) and E(z1, z2) don't exist, rewire 2 edges
      repeat {
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

      if (transitivity(g.cand) > transitivity(g)) {
        break
      } else {
        # Undo the changes that were just made, and try again
        g.cand <- g.cand - E(g.cand, P=c(e$y1, e$y2))
        g.cand <- g.cand - E(g.cand, P=c(e$z1, e$z2))
        g.cand <- g.cand + edges(c(e$y1, e$z1)) + edges(c(e$y2, e$z2))
      }
    }

    cur.iter <- cur.iter + 1
    cl.all[cur.iter] <- transitivity(g.cand)
    g <- g.cand
    loop.time[cur.iter] <- as.numeric(Sys.time() - start.time)
  }

  return(list(g=g, iters=cur.iter, cl=cl.all, time=loop.time))#, g.all))
}
