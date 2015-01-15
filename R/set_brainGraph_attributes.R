#' Set a number of graph and vertex attributes useful in MRI analyses
#'
#' This function will set a number of graph, vertex, and edge attributes of a
#' given igraph object.
#'
#' @param g An igraph object
#' @param atlas A character vector indicating which atlas was used for the nodes
#' @param coords A matrix of the spatial coordinates
#' @param rand Logical indicating whether the graph is random or not
#'
#' @export
#'
#' @return g A copy of the same graph, with the following attributes:
#' \item{Graph-level}{Atlas, density, connected component sizes, diameter,
#' number of triangles, transitivity, average path length, assortativity,
#' clique sizes, global & local efficiency, modularity, hub score, and rich-club
#' coefficient}
#' \item{Vertex-level}{Degree, strength, betweenness centrality, eigenvector
#' centrality, transitivity (local), coreness, local & nodal efficiency, color
#' (community), color (lobe), color (component), membership (community),
#' membership (component), participation coefficient, within-module degree
#' z-score, and coordinates (x, y, and z)}
#' \item{Edge-level}{Color (community), color (lobe), color (component), edge
#' betweenness}
#'
#' @seealso \code{\link{components}, \link{graph.motifs}, \link{diameter},
#' \link{cliques}, \link{centralization.betweenness}, \link{edge.betweenness},
#' \link{centralization.evcent}, \link{subgraph.centrality}, \link{hub.score},
#' \link{authority.score}, \link{transitivity}, \link{average.path.length},
#' \link{assortativity.degree}, \link{graph.efficiency}, \link{rich.club.coeff},
#' \link{edge.betweenness.community}, \link{color.edges}, \link{part.coeff},
#' \link{within.module.deg.z.score},\link{graph.coreness}}

set.brainGraph.attributes <- function(g, atlas=NULL, coords=NULL, rand=FALSE) {
  V(g)$degree <- degree(g)
  # Graph-level attributes
  #-----------------------------------------------------------------------------
  g$density <- round(graph.density(g), digits=3)

  # Connected components
  clusts <- components(g)
  comps <- rev(table(clusts$csize))
  g$conn.comp <- data.frame(size=as.numeric(names(comps)),
                            number=unclass(unname(comps)))
  #g$cliques <- table(sapply(cliques(g), length))
  g$num.tri <- graph.motifs(g)[4]
  g$diameter <- diameter(g)

  g$transitivity <- transitivity(g)
  g$Cp <- transitivity(g, type='localaverage')
  g$Lp <- average.path.length(g)
  g$assortativity <- assortativity.degree(g)
  g$E.global <- graph.efficiency(g, type='global')

  # Get the rich club coeff for all possible degree values
  R <- lapply(1:max(V(g)$degree), function(x) rich.club.coeff(g, x))
  coef <- sapply(R, function(x) x$coeff)
  Nk <- sapply(R, function(x) x$Nk)
  Ek <- sapply(R, function(x) x$Ek)
  g$rich <- data.frame(R=round(coef, 4), Nk=Nk, Ek=Ek)

  if (is.directed(g)) {
    hubs <- hub.score(g)
    g$hub.score <- hubs$value
    authorities <- authority.score(g)
    g$authority.score <- authorities$value
    V(g)$hub.score <- hubs$vector
    V(g)$authority.score <- authorities$vector
  }

  # Vertex-level attributes
  #-----------------------------------------------------------------------------
  # Give each vertex a 'lobe' attribute, and colors for each lobe
  if (!is.null(atlas)) {
    lobe.cols <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orange',
                   'lightgreen')
    g$atlas <- atlas
    V(g)$color.lobe <- color.vertices(atlas=atlas)

    atlas.list <- eval(parse(text=atlas))
    lobes <- lobe.hemi <- vector('integer', length=vcount(g))

    lobes[atlas.list$frontal] <- 1
    lobes[atlas.list$parietal] <- 2
    lobes[atlas.list$temporal] <- 3
    lobes[atlas.list$occipital] <- 4
    lobes[atlas.list$insula] <- 5
    lobes[atlas.list$cingulate] <- 6

    lobe.hemi[atlas.list$frontal.lh] <- 1
    lobe.hemi[atlas.list$frontal.rh] <- 2
    lobe.hemi[atlas.list$parietal.lh] <- 3
    lobe.hemi[atlas.list$parietal.rh] <- 4
    lobe.hemi[atlas.list$temporal.lh] <- 5
    lobe.hemi[atlas.list$temporal.rh] <- 6
    lobe.hemi[atlas.list$occipital.lh] <- 7
    lobe.hemi[atlas.list$occipital.rh] <- 8
    lobe.hemi[atlas.list$insula.lh] <- 9
    lobe.hemi[atlas.list$insula.rh] <- 10
    lobe.hemi[atlas.list$cingulate.lh] <- 11
    lobe.hemi[atlas.list$cingulate.rh] <- 12

    if (atlas %in% c('dkt', 'dk')) {
      V(g)$circle.layout <- with(atlas.list,
                                 c(frontal.lh, insula.lh, cingulate.lh, 
                                   temporal.lh, parietal.lh, occipital.lh, 
                                   occipital.rh, parietal.rh, temporal.rh,
                                   cingulate.rh, insula.rh, frontal.rh))

    } else if (atlas == 'aal90') {
      lobes[atlas.list$limbic] <- 6
      lobes[atlas.list$scgm] <- 7
      V(g)$circle.layout <- with(atlas.list,
                                 c(frontal.lh, insula.lh, limbic.lh, scgm.lh,
                                   temporal.lh, parietal.lh, occipital.lh,
                                   occipital.rh, parietal.rh, temporal.rh,
                                   scgm.rh, limbic.rh, insula.rh, frontal.rh))
      lobe.hemi[atlas.list$limbic.lh] <- 11
      lobe.hemi[atlas.list$limbic.rh] <- 12
      lobe.hemi[atlas.list$scgm.lh] <- 13
      lobe.hemi[atlas.list$scgm.rh] <- 14

    } else if (atlas %in% c('lpba40', 'hoa112', 'brainsuite')) {
      lobes[atlas.list$scgm] <- 7
      V(g)$circle.layout <- with(atlas.list,
                                 c(frontal.lh, insula.lh, cingulate.lh, scgm.lh,
                                   temporal.lh, parietal.lh, occipital.lh, 
                                   occipital.rh, parietal.rh, temporal.rh,
                                   scgm.rh, cingulate.rh, insula.rh, frontal.rh))
      lobe.hemi[atlas.list$scgm.lh] <- 13
      lobe.hemi[atlas.list$scgm.rh] <- 14
    }

    V(g)$lobe <- lobes
    V(g)$lobe.hemi <- lobe.hemi

    g$assortativity.lobe <- assortativity(g, V(g)$lobe)
    g$assortativity.lobe.hemi <- assortativity(g, V(g)$lobe.hemi)
    E(g)$color.lobe <- color.edges(g, lobes=V(g)$lobe, lobe.cols=lobe.cols)
  }

  # Add the spatial coordinates for plotting over the brain
  if (length(coords) > 0) {
    V(g)$x <- coords[, 1]
    V(g)$y <- coords[, 2]
    V(g)$z <- coords[, 3]
    V(g)$name <- rownames(coords)
  }

  if (is.weighted(g)) {
    V(g)$strength <- graph.strength(g)
  }
  V(g)$btwn.cent <- centralization.betweenness(g)$res
  V(g)$ev.cent <- centralization.evcent(g)$vector
  V(g)$subgraph.cent <- subgraph.centrality(g)
  V(g)$coreness <- graph.coreness(g)

  V(g)$transitivity <- transitivity(g, type='local', isolates='zero')

  V(g)$E.local <- graph.efficiency(g, type='local')
  V(g)$E.nodal <- graph.efficiency(g, type='nodal')
  g$E.local <- mean(V(g)$E.local)

  if (rand == TRUE) {
    g$mod <- max(multilevel.community(g)$modularity)
  } else {
    # Community stuff
    comm <- multilevel.community(g)
    V(g)$comm <- comm$membership
    vcolors <- color.vertices(V(g)$comm)
    V(g)$color.comm <- vcolors[V(g)$comm]
    E(g)$color.comm <- color.edges(g, V(g)$comm)

    V(g)$circle.layout.comm <- order(V(g)$comm,
                                     V(g)$degree)

    V(g)$PC <- part.coeff(g, V(g)$comm)
    V(g)$z.score <- within.module.deg.z.score(g, V(g)$comm)
    g$mod <- max(comm$modularity)
  }

  V(g)$comp <- clusts$membership
  vcolors <- color.vertices(clusts$membership)
  V(g)$color.comp <- vcolors[clusts$membership]

  # Edge attributes
  #-----------------------------------------------------------------------------
  E(g)$color.comp <- color.edges(g, clusts$membership)
  E(g)$btwn <- edge.betweenness(g)

  g
}
