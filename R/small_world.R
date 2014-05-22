#' Get measures to calculate small-worldness.
#'
#' This function will calculate the characteristic path length and clustering
#' coefficient, which are used to calculate small-worldness.
#'
#' @param g the adjacency graph (or a list of graphs)
#' @param rand list of output from sim.random.graph
#' @param threshes Vector of thresholds (either by cost or correlation coeff)
#' @export
#'
#' @return A data frame with the following components:
#' \item{cost}{The range of cost thresholds used.}
#' \item{Lp}{The characteristic path length.}
#' \item{Cp}{The clustering coefficient.}
#' \item{Lp.rand}{The mean characteristic path length of the random graphs with
#' the same degree distribution as g.}
#' \item{Cp.rand}{The mean clustering coefficient of the random graphs with
#' the same degree distribution as g.}
#' \item{sigma}{The small-world measure of the graph.}

small.world <- function(g, rand, threshes) {
  Lp <- vapply(g, average.path.length, numeric(1))
  Cp <- vapply(g, transitivity, numeric(1))

  Lp.rand <- vapply(rand, function(x) mean(x$Lp.rand), numeric(1))
  Cp.rand <- vapply(rand, function(x) mean(x$Cp.rand), numeric(1))

  sigma <- (Cp / Cp.rand) / (Lp / Lp.rand)
  return(data.frame(cost=threshes, Lp=Lp, Cp=Cp, Lp.rand=Lp.rand,
                    Cp.rand=Cp.rand, sigma=sigma))
}
