#' Get measures to calculate small-worldness.
#'
#' This function will calculate the characteristic path length and clustering
#' coefficient, which are used to calculate small-worldness.
#'
#' @param g The graph (or list of graphs) of interest
#' @param rand List of equivalent random graphs (output from
#' \code{\link{sim.rand.graph.par}}
#' @export
#'
#' @return A data frame with the following components:
#' \item{density}{The range of density thresholds used.}
#' \item{Lp}{The characteristic path length.}
#' \item{Cp}{The clustering coefficient.}
#' \item{Lp.rand}{The mean characteristic path length of the random graphs with
#' the same degree distribution as g.}
#' \item{Cp.rand}{The mean clustering coefficient of the random graphs with
#' the same degree distribution as g.}
#' \item{sigma}{The small-world measure of the graph.}
#'
#' @references Watts D.J., Strogatz S.H. (1998) \emph{Collective dynamics of
#' 'small-world' networks}. Nature, 393:440-442.

small.world <- function(g, rand) {
  Lp <- vapply(g, function(x) x$Lp, numeric(1))
  Cp <- vapply(g, function(x) x$Cp, numeric(1))
  density <- vapply(g, function(x) x$density, numeric(1))

  Lp.rand <- colMeans(vapply(rand,
      function(x) vapply(x,
          function(y) y$Lp, numeric(1)), numeric(length(rand[[1]]))))
  Cp.rand <- colMeans(vapply(rand,
      function(x) vapply(x,
          function(y) y$Cp, numeric(1)), numeric(length(rand[[1]]))))

  sigma <- (Cp / Cp.rand) / (Lp / Lp.rand)
  return(data.frame(density=density, Lp=Lp, Cp=Cp, Lp.rand=Lp.rand,
                    Cp.rand=Cp.rand, sigma=sigma))
}
