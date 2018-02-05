#' Simulate N random graphs w/ same clustering and degree sequence as the input.
#'
#' \code{sim.rand.graph.par} simulates \code{N} simple random graphs with the
#' same clustering (optional) and degree sequence as the input. Essentially a
#' wrapper for \code{\link[igraph]{sample_degseq}} (or, if you want to match by
#' clustering, \code{\link{sim.rand.graph.clust}}) and
#' \code{\link{set_brainGraph_attr}}. It uses \code{\link[foreach]{foreach}} for
#' parallel processing.
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
#' @return A \emph{list} of \emph{N} random graphs with some additional vertex
#'   and graph attributes
#'
#' @name RandomGraphs
#' @aliases sim.rand.graph.par
#' @rdname random_graphs
#'
#' @family Random graph functions
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
  iters <- max(10*ecount(g), 1e4)
  if (isTRUE(clustering)) {
    r <- foreach(i=seq_len(N), .packages=c('igraph', 'brainGraph')) %dopar% {
      tmp <- sim.rand.graph.clust(g, rewire.iters=iters, ...)
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
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- rewire(g, keeping_degseq(loops=FALSE, iters))
        tmp <- set_brainGraph_attr(tmp, rand=TRUE)
        tmp
      }
    }
  }
}

#' Simulate a random graph with given degree sequence and clustering.
#'
#' \code{sim.rand.graph.clust} simulates a random graph with a given degree
#' sequence \emph{and} clustering coefficient. Increasing the \code{max.iters}
#' value will result in a closer match of clustering with the observed graph.
#'
#' @param rewire.iters Integer; number of rewiring iterations for the initial
#'   graph randomization (default: 1e4)
#' @param cl The clustering measure (default: transitivity)
#' @param max.iters The maximum number of iterations to perform; choosing a
#'   lower number may result in clustering that is further away from the
#'   observed graph's (default: 100)
#' @export
#'
#' @return An \code{igraph} graph object
#'
#' @aliases sim.rand.graph.clust
#' @rdname random_graphs
#' @family Random graph functions
#'
#' @seealso \code{\link[igraph]{rewire}, \link[igraph]{transitivity},
#'   \link[igraph]{keeping_degseq}}
#'
#' @references Bansal S., Khandelwal S., Meyers L.A. (2009) \emph{Exploring
#' biological network structure with clustered random networks}. BMC
#' Bioinformatics, 10:405-421.

sim.rand.graph.clust <- function(g, rewire.iters=1e4, cl=g$transitivity, max.iters=100) {
  g <- rewire(g, keeping_degseq(loops=FALSE, rewire.iters))

  cur.iter <- 0
  while ((transitivity(g) < cl) & (cur.iter < max.iters)) {
    repeat {
      g.cand <- g
      A <- as_adj(g.cand, sparse=FALSE, names=FALSE)
      degs <- colSums(A)
      degs.large <- which(degs > 1)

      # If E(y1, y2) and E(z1, z2) don't exist, rewire 2 edges
      repeat {
        e <- choose.edges(A, degs.large)
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

#' Select edges for re-wiring.
#'
#' This function selects edges to be re-wired when simulating random graphs
#' while controlling for \emph{clustering}. It is based on the algorithm by
#' Bansal et al. (2009), BMC Bioinformatics.
#'
#' @param A Numeric (adjacency) matrix
#' @param degs.large Integer vector of vertex numbers with degree greater than
#'   one
#'
#' @return A data frame with four elements; two edges will be removed and two
#'   edges will be added between the four vertices.
#'
#' @keywords internal
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Bansal S., Khandelwal S., Meyers L.A. (2009) \emph{Exploring
#'   biological network structure with clustered random networks}. BMC
#'   Bioinformatics, 10:405-421.

choose.edges <- function(A, degs.large) {

  # Uniformly select random node with degree > 1
  #=============================================================================
  get.node <- function(A, degs.large) {
    repeat {
      x <- sample(degs.large, 1)
      neighb <- intersect(which(A[x, ] == 1), degs.large)
      if (length(neighb) >= 2) return(list(x, neighb))
    }
  }
  # Uniformly select 2 random neighbors with degree > 1
  #=============================================================================
  get.neighbors <- function(nbrhood) {
    y <- sample(nbrhood, 2)
    return(data.frame(y1=y[1], y2=y[2]))
  }
  #=============================================================================
  # Uniformly select a random neighbor from the y's with degree > 1
  #=============================================================================
  get.neighbors.z1 <- function(A, node, y1) {
    y1.neighb <- which(A[y1, ] == 1)
    choices <- setdiff(y1.neighb, node)
    if (length(choices) <= 1) {
      return(choices)
    } else {
      return(sample(choices, 1))
    }
  }
  #=============================================================================
  get.neighbors.z2 <- function(A, node, degs.large, y2, z1) {
    y2.neighb <- which(A[y2, ] == 1)
    repeat {
      choices <- setdiff(y2.neighb, c(node, z1))
      if (length(choices) == 1) {
        return(choices)
      } else if (length(choices) > 1) {
        return(sample(choices, 1))
      }
      tmp <- get.node(A, degs.large)
      node <- tmp[[1]]
      neighb <- tmp[[2]]
      n <- get.neighbors(neighb)
      y1 <- n$y1
      y2 <- n$y2

      z1 <- get.neighbors.z1(A, node, y1)
      get.neighbors.z2(A, node, degs.large, y2, z1)
    }
  }
  #=============================================================================
  #=============================================================================
  tmp <- get.node(A, degs.large)
  x <- tmp[[1]]
  neighb <- tmp[[2]]

  n <- get.neighbors(neighb)
  y1 <- n$y1
  y2 <- n$y2
  z1 <- get.neighbors.z1(A, x, y1)
  z2 <- get.neighbors.z2(A, x, degs.large, y2, z1)

  return(data.frame(y1, y2, z1, z2))
}
