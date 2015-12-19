#' Calculate graph small-worldness
#'
#' This function will calculate the characteristic path length and clustering
#' coefficient, which are used to calculate small-worldness.
#'
#' @param g The graph (or list of graphs) of interest
#' @param rand List of (lists of) equivalent random graphs (output from
#' \code{\link{sim.rand.graph.par}})
#' @export
#'
#' @return A data frame with the following components:
#' \item{density}{The range of density thresholds used.}
#' \item{N}{The number of random graphs that were generated.}
#' \item{Lp}{The characteristic path length.}
#' \item{Cp}{The clustering coefficient.}
#' \item{Lp.rand}{The mean characteristic path length of the random graphs with
#' the same degree distribution as g.}
#' \item{Cp.rand}{The mean clustering coefficient of the random graphs with
#' the same degree distribution as g.}
#' \item{Lp.norm}{The normalized characteristic path length.}
#' \item{Cp.norm}{The normalized clustering coefficient.}
#' \item{sigma}{The small-world measure of the graph.}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Watts D.J., Strogatz S.H. (1998) \emph{Collective dynamics of
#' 'small-world' networks}. Nature, 393:440-442.

small.world <- function(g, rand) {
  if (is.igraph(g)) g <- list(g)  # Single graph at a single density
  Lp <- vapply(g, function(x) graph_attr(x, 'Lp'), numeric(1))
  Cp <- vapply(g, function(x) graph_attr(x, 'Cp'), numeric(1))
  densities <- vapply(g, function(x) graph_attr(x, 'density'), numeric(1))

  if (is.igraph(rand[[1]])) {
    Lp.rand <- mean(vapply(rand, function(x) graph_attr(x, 'Lp'), numeric(1)))
    Cp.rand <- mean(vapply(rand, function(x) graph_attr(x, 'Cp'), numeric(1)))
    N <- length(rand)
  } else {
    if (length(rand[[1]]) == 1) {  # If there's 1 rand graph for each density
      Lp.rand <- vapply(rand,
        function(x) vapply(x,
            function(y) graph_attr(y, 'Lp'), numeric(1)),
                               numeric(1))
      Cp.rand <- vapply(rand,
        function(x) vapply(x,
            function(y) graph_attr(y, 'Cp'), numeric(1)),
                               numeric(1))
      N <- 1
    } else {
      Lp.rand <- colMeans(vapply(rand,
          function(x) vapply(x,
              function(y) graph_attr(y, 'Lp'), numeric(1)),
                                 numeric(length(rand[[1]]))))
      Cp.rand <- colMeans(vapply(rand,
          function(x) vapply(x,
              function(y) graph_attr(y, 'Cp'), numeric(1)),
                                 numeric(length(rand[[1]]))))
      N <- lengths(rand)
    }
  }
  Cp.norm <- Cp / Cp.rand
  Lp.norm <- Lp / Lp.rand
  sigma <- Cp.norm / Lp.norm
  return(data.table(density=densities, N, Lp, Cp, Lp.rand, Cp.rand, Lp.norm,
                    Cp.norm, sigma))
}
