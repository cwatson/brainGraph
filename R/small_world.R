#' Calculate graph small-worldness
#'
#' \code{small.world} calculates the normalizaed characteristic path length and
#' clustering coefficient based on observed and random graphs. These are used to
#' calculate the small-world coefficient \eqn{\sigma}.
#'
#' @param g.list A \code{brainGraphList} object or list of graphs
#' @param rand List of (lists of) equivalent random graphs (output from
#' \code{\link{sim.rand.graph.par}})
#' @export
#'
#' @return A \code{data.table} with the following components:
#' \item{density}{The range of density thresholds used.}
#' \item{N}{The number of random graphs that were generated.}
#' \item{Lp,Lp.rand,Lp.norm}{The observed, average random, and normalized
#'   characteristic path length.}
#' \item{Cp,Cp.rand,Cp.norm}{The observed, average random, and normalized
#'   clustering coefficient.}
#' \item{sigma}{The small-world measure of the graph.}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Watts, D.J. and Strogatz S.H. (1998) Collective dynamics of
#'   "small-world" networks. \emph{Nature}, \bold{393}, 440--442.
#'   \url{https://dx.doi.org/10.1515/9781400841356.301}

small.world <- function(g.list, rand) {
  Study.ID <- NULL
  if (is_igraph(g.list) || is.brainGraph(g.list)) {
    g.list <- list(g.list)
  } else if (inherits(g.list, 'brainGraphList')) {
    g.list <- g.list$graphs
    rand <- rand$graphs
  } else if (!is.list(g.list)) {
    stop('Input should be a list, igraph, brainGraph, or brainGraphList object.')
  }

  Lp <- vapply(g.list, graph_attr, numeric(1), 'Lp')
  Cp <- vapply(g.list, graph_attr, numeric(1), 'Cp')
  densities <- vapply(g.list, graph_attr, numeric(1), 'density')

  N <- lengths(rand)
  Lp.rand <- colMeans(vapply(rand, vapply, numeric(N[1]), graph_attr, numeric(1), 'Lp'))
  Cp.rand <- colMeans(vapply(rand, vapply, numeric(N[1]), graph_attr, numeric(1), 'Cp'))

  Cp.norm <- Cp / Cp.rand
  Lp.norm <- Lp / Lp.rand
  sigma <- Cp.norm / Lp.norm
  DT <- data.table(density=densities, N, Lp, Cp, Lp.rand, Cp.rand,
                   Lp.norm, Cp.norm, sigma)
  if (!is.null(names(g.list))) DT[, Study.ID := names(g.list)]
  attrs <- c('Group', 'threshold')
  for (x in attrs) {
    if (x %in% graph_attr_names(g.list[[1]])) {
      DT[, eval(x) := sapply(g.list, graph_attr, x)]
    }
  }
  return(DT)
}
