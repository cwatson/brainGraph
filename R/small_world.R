#' Calculate graph small-worldness
#'
#' \code{small.world} calculates the normalized characteristic path length and
#' clustering coefficient based on observed and random graphs, used to calculate
#' the small-world coefficient \eqn{\sigma}.
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
#'   \doi{10.1038/30918}

small.world <- function(g.list, rand) {
  level <- NULL
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')

  if (is_igraph(g.list) || is.brainGraph(g.list)) {
    g.list <- list(g.list)
  } else if (is.brainGraphList(g.list)) {
    level <- g.list$level
    g.list <- g.list$graphs
    rand <- rand$graphs
  } else if (!is.list(g.list)) {
    stop('Input should be a list, igraph, brainGraph, or brainGraphList object.')
  }

  Lp <- vapply(g.list, graph_attr, numeric(1L), 'Lp')
  Cp <- vapply(g.list, graph_attr, numeric(1L), 'Cp')
  densities <- vapply(g.list, graph_attr, numeric(1L), 'density')

  N <- lengths(rand)
  Lp.rand <- colMeans(vapply(rand, vapply, numeric(N[1L]), graph_attr, numeric(1L), 'Lp'))
  Cp.rand <- colMeans(vapply(rand, vapply, numeric(N[1L]), graph_attr, numeric(1L), 'Cp'))

  Cp.norm <- Cp / Cp.rand
  Lp.norm <- Lp / Lp.rand
  sigma <- Cp.norm / Lp.norm
  DT <- data.table(density=densities, N, Lp, Cp, Lp.rand, Cp.rand,
                   Lp.norm, Cp.norm, sigma)
  if (!is.null(names(g.list))) {
    DT[, eval(sID) := names(g.list)]
  }
  attrs <- c('Group', 'threshold')
  for (x in attrs) {
    if (x %in% graph_attr_names(g.list[[1L]])) {
      DT[, eval(x) := sapply(g.list, graph_attr, x)]
    }
  }
  if (hasName(DT, 'Group')) setnames(DT, 'Group', gID)
  if (level == 'group' && DT[, all(get(sID) == get(gID))]) DT[, eval(sID) := NULL]
  return(DT)
}
