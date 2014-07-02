#' Set a number of graph and vertex attributes useful in MRI analyses
#'
#' This function will set a number of graph and vertex attributes of a given
#' igraph object.
#'
#' @param g A graph
#' @param coords A matrix of the spatial coordinates
#' @param rand Logical indicating whether the graph is random or not
#'
#' @export
#'
#' @return g A copy of the same graph, with the following attributes:
#' \item{Graph-level}{Density, largest connected component size, number of
#' triangles, transitivity, average path length, assortativity, global
#' efficiency, modularity, and rich-club coefficient}
#' \item{Vertex-level}{Degree, betweenness centrality, eigenvector centrality,
#' transitivity (local), local efficiency, color, community membership,
#' participation coefficient, within-module degree z-score, and coordinates (x,
#' y, and z)}
#' \item{Edge-level}{Color}
#'
#' @seealso \code{\link{clusters}, \link{graph.motifs},
#' \link{centralization.betweenness}, \link{centralization.evcent},
#' \link{transitivity}, \link{average.path.length}, \link{assortativity.degree},
#' \link{global.eff}, \link{local.eff}, \link{rich.club.coeff},
#' \link{edge.betweenness.community}, \link{community.measures},
#' \link{color.edges}, \link{part.coeff}, \link{within.module.degree.z.score}}

set.brainGraph.attributes <- function(g, coords=NULL, rand=FALSE) {
  if (length(coords) > 0) {
    V(g)$x <- coords[, 1]
    V(g)$y <- coords[, 2]
    V(g)$z <- coords[, 3]
    V(g)$name <- rownames(coords)
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

  # Currently, takes the top 10% of nodes for rich club
  deg.thresh <- Nv - ceiling(0.1 * Nv)
  g$rich <- rich.club.coeff(g, sort(degree(g))[deg.thresh])$coeff

  if (rand == TRUE) {
    g$mod <- max(edge.betweenness.community(g)$modularity)
  } else {
    # Community stuff
    comm <- community.measures(g)
    V(g)$color <- comm$vcolors[comm$community$membership]
    E(g)$color <- color.edges(g, comm)
    V(g)$comm <- comm[[1]]$membership

    A <- get.adjacency(g, sparse=F)
    V(g)$PC <- part.coeff(A, g, comm[[1]]$membership)
    V(g)$z.score <- within.module.deg.z.score(A, g, comm[[1]]$membership)
    g$mod <- max(comm[[1]]$modularity)
  }

  g
}
