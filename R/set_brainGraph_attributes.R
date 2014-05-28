#' Set a number of graph and vertex attributes useful in MRI analyses
#'
#' This function will set a number of graph and vertex attributes of a given
#' igraph object. It is required that a community structure has already been
#' determined, and \code{\link{part.coeff}} and
#' \code{\link{within.module.deg.z.score}} have been calculated.
#'
#' @param g A graph
#' @export
#'
#' @return g A copy of the same graph, with the attributes added.

set.brainGraph.attributes <- function(g, vcols, ecols, PC, z, coords=NULL,
                                      rand=FALSE) {
  if (length(coords) > 0) {
    V(g)$x <- coords[, 1]
    V(g)$y <- coords[, 2]
    V(g)$z <- coords[, 3]
  }

  Nv <- vcount(g)
  g$density <- round(2 * ecount(g) / (Nv * (Nv - 1)), digits=3)
  g$conn.comp <- max(clusters(g)$csize)
  g$num.tri <- graph.motifs(g)[4]
  V(g)$degree <- degree(g)

  V(g)$btwn.cent <- centralization.betweenness(g)$res
  V(g)$ev.cent <- centralization.evcent(g)$vector

  t <- transitivity(g, type='local')
  V(g)$transitivity <- ifelse(is.nan(t), 0, t)
  g$Cp <- transitivity(g)

  g$Lp <- average.path.length(g)

  g$assortativity <- assortativity.degree(g)

  g$g.eff <- global.eff(g)
  V(g)$l.eff <- local.eff(g)

  if (rand == TRUE) {
    deg.thresh <- Nv - ceiling(0.1 * Nv)
    g$mod <- max(fastgreedy.community(g)$modularity)
    g$rich <- rich.club.coeff(g, sort(degree(g))[deg.thresh])$coeff
  } else {
    V(g)$color <- vcols
    E(g)$color <- ecols
    V(g)$PC <- PC
    V(g)$z.score <- z
  }

  g
}
