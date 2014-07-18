#' Set a number of graph and vertex attributes useful in MRI analyses
#'
#' This function will set a number of graph and vertex attributes of a given
#' igraph object.
#'
#' @param g A graph
#' @param atlas A character vector indicating which atlas was used for the nodes
#' @param coords A matrix of the spatial coordinates
#' @param rand Logical indicating whether the graph is random or not
#'
#' @export
#'
#' @return g A copy of the same graph, with the following attributes:
#' \item{Graph-level}{Atlas, density, connected component sizes, diameter,
#' number of triangles, transitivity, average path length, assortativity,
#' global efficiency, modularity, and rich-club coefficient}
#' \item{Vertex-level}{Degree, betweenness centrality, eigenvector centrality,
#' transitivity (local), coreness, local efficiency, color, community membership,
#' participation coefficient, within-module degree z-score, and coordinates (x,
#' y, and z)}
#' \item{Edge-level}{Color, lobe color}
#'
#' @seealso \code{\link{clusters}, \link{graph.motifs}, \link{diameter},
#' \link{centralization.betweenness}, \link{centralization.evcent},
#' \link{subgraph.centrality}, \link{hub.score}, \link{authority.score},
#' \link{transitivity}, \link{average.path.length}, \link{assortativity.degree},
#' \link{global.eff}, \link{local.eff}, \link{rich.club.coeff},
#' \link{edge.betweenness.community}, \link{color.edges}, \link{part.coeff},
#' \link{within.module.deg.z.score},\link{graph.coreness}}

set.brainGraph.attributes <- function(g, atlas=NULL, coords=NULL, rand=FALSE) {
  Nv <- vcount(g)

  # Graph-level attributes
  #-----------------------------------------------------------------------------
  g$density <- round(2 * ecount(g) / (Nv * (Nv - 1)), digits=3)
  comps <- rev(table(clusters(g)$csize))
  g$conn.comp <- data.frame(size=as.numeric(names(comps)),
                            number=unclass(unname(comps)))
  g$num.tri <- graph.motifs(g)[4]
  g$diameter <- diameter(g)

  t <- transitivity(g, type='local')
  g$Cp <- transitivity(g)
  g$Lp <- average.path.length(g)
  g$assortativity <- assortativity.degree(g)
  g$g.eff <- global.eff(g)

  # Get the rich club coeff for all possible degree values
  R <- lapply(1:max(degree(g)), function(x) rich.club.coeff(g, x))
  coef <- sapply(R, function(x) x$coeff)
  Nk <- sapply(R, function(x) x$Nk)
  Ek <- sapply(R, function(x) x$Ek)
  g$rich <- data.frame(R=round(coef, 4), Nk=Nk, Ek=Ek)

  # Vertex-level attributes
  #-----------------------------------------------------------------------------
  # Give each vertex a 'lobe' attribute, and colors for each lobe
  if (!is.null(atlas)) {
    lobe.cols <- c('red', 'green', 'blue', 'magenta', 'yellow')
    g$atlas <- atlas
    V(g)$lobe.color <- color.vertices(atlas=atlas)

    if (atlas == 'dkt') {
      lobes <- vector('integer', length=Nv)
      lobes[dkt$frontal] <- 1
      lobes[dkt$parietal] <- 2
      lobes[dkt$temporal] <- 3
      lobes[dkt$occipital] <- 4
      lobes[dkt$insula] <- 5
      V(g)$lobe <- lobes

      V(g)$circle.layout <- c(dkt$frontal.lh, dkt$parietal.lh, dkt$temporal.lh,
                              dkt$occipital.lh, dkt$insula.lh, dkt$frontal.rh,
                              dkt$parietal.rh, dkt$temporal.rh, dkt$occipital.rh,
                              dkt$insula.rh)

    } else if (atlas == 'dk') {
      lobes <- vector('integer', length=Nv)
      lobes[dk$frontal] <- 1
      lobes[dk$parietal] <- 2
      lobes[dk$temporal] <- 3
      lobes[dk$occipital] <- 4
      lobes[dk$insula] <- 5
      V(g)$lobe <- lobes

      V(g)$circle.layout <- c(dk$frontal.lh, dk$parietal.lh, dk$temporal.lh,
                              dk$occipital.lh, dk$insula.lh, dk$frontal.rh,
                              dk$parietal.rh, dk$temporal.rh, dk$occipital.rh,
                              dk$insula.rh)
    }
    E(g)$lobe.color <- color.edges(g, lobes=V(g)$lobe, lobe.cols=lobe.cols)
  }

  # Add the spatial coordinates for plotting over the brain
  if (length(coords) > 0) {
    V(g)$x <- coords[, 1]
    V(g)$y <- coords[, 2]
    V(g)$z <- coords[, 3]
    V(g)$name <- rownames(coords)
  }

  V(g)$degree <- degree(g)
  V(g)$btwn.cent <- centralization.betweenness(g)$res
  V(g)$ev.cent <- centralization.evcent(g)$vector
  V(g)$subgraph.cent <- subgraph.centrality(g)
  # Calculate the coreness of each vertex
  V(g)$coreness <- graph.coreness(g)

  V(g)$transitivity <- ifelse(is.nan(t), 0, t)
  V(g)$l.eff <- local.eff(g)

  # Calculate both the hub and authority scores for each vertex
  # Check if 'g' is directed; if not, both scores are the same
  V(g)$hub.score <- hub.score(g)$vector
  if (is.directed(g)) {
    V(g)$authority.score <- authority.score(g)$vector
  }

  if (rand == TRUE) {
    g$mod <- max(edge.betweenness.community(g)$modularity)
  } else {
    # Community stuff
    comm <- edge.betweenness.community(g)
    V(g)$comm <- comm$membership
    vcolors <- color.vertices(comm)
    V(g)$color <- vcolors[V(g)$comm]
    E(g)$color <- color.edges(g, V(g)$comm)

    V(g)$circle.layout.comm <- order(V(g)$comm,
                                     V(g)$degree)

    V(g)$PC <- part.coeff(g, V(g)$comm)
    V(g)$z.score <- within.module.deg.z.score(g, V(g)$comm)
    g$mod <- max(comm$modularity)
  }

  g
}
