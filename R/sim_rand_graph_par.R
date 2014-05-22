#' Simulate N random graphs w/ same clustering and degree sequence as the input.
#'
#' This function will simulate N simple random graphs with the same clustering
#' and degree sequence as the input. It will also calculate (and return) avg.
#' path length, transitivity, mean degree, and modularity. Essentially a wrapper
#' for \code{\link{sim.rand.graph.clust}}. It uses \code{\link{foreach}} to
#' speed it up.
#'
#' @param g A graph with the characteristics for simulation of random graphs
#' @param N The number of iterations
#' @export
#'
#' @return A list with the following components:
#' \item{Lp.rand}{A vector of the characteristic path lengths for the N random
#' graphs that were created.}
#' \item{Cp.rand}{A vector of the clustering coefficients for the N random
#' graphs that were created.}
#' \item{g.eff}{A vector of the global efficiencies for the N random
#' graphs that were created.}
#' \item{l.eff}{A vector of the mean local efficiencies for the N random
#' graphs that were created.}
#' \item{rich}{A vector of the rich club coefficients for the N random
#' graphs that were created.}
#' \item{mod}{A vector of the modularity for the N random graphs.}
#' \item{Ek}{The number of edges for each graph.}
#' \item{triangles}{The number of triangles for each graph.}

sim.rand.graph.par <- function(g, N) {
  deg.dist <- V(g)$degree
  deg.thresh <- length(deg.dist) - ceiling(0.1 * length(deg.dist))
  cl <- g$cl.coeff

  Lp <- vector(length=N)
  Cp <- vector(length=N)
  g.eff <- vector(length=N)
  l.eff <- vector(length=N)
  rich <- vector(length=N)
  mod <- vector(length=N)

  x <- matrix(nrow=6, ncol=N)
  x <- foreach (i=1:N, .combine=cbind,
                .packages=c('igraph', 'brainGraph')) %dopar% {
    tmp <- sim.rand.graph.clust(deg.dist, cl)
    Lp <- average.path.length(tmp)
    Cp <- transitivity(tmp)
    g.eff <- global.eff(tmp)
    l.eff <- mean(local.eff(tmp))
    rich <- rich.club.coeff(tmp, sort(degree(tmp))[deg.thresh])$coeff
    mod <- max(fastgreedy.community(tmp)$modularity)
    Ek <- ecount(tmp)
    num.tri <- graph.motifs(tmp)[4]

    c(Lp, Cp, g.eff, l.eff, rich, mod, Ek, num.tri)
  }

  list(Lp.rand=x[1, ], Cp.rand=x[2, ], g.eff=x[3, ], l.eff=x[4, ], rich=x[5, ],
       mod=x[6, ], Ek=x[7, ], triangles=x[8, ])
}
