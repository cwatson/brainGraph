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
  # For AAL90: frontal, parietal, temporal, occipital, limbic, insula, SCGM
  if (!is.null(atlas)) {
    lobe.cols <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orange',
                   'lightgreen')
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
    } else if (atlas == 'aal90') {
      lobes <- vector('integer', length=Nv)
      lobes[aal90$frontal] <- 1
      lobes[aal90$parietal] <- 2
      lobes[aal90$temporal] <- 3
      lobes[aal90$occipital] <- 4
      lobes[aal90$insula] <- 5
      lobes[aal90$limbic] <- 6
      lobes[aal90$scgm] <- 7
      V(g)$lobe <- lobes

      V(g)$circle.layout <- c(aal90$frontal.lh, aal90$parietal.lh, aal90$temporal.lh,
                              aal90$occipital.lh, aal90$insula.lh, aal90$limbic.lh,
                              aal90$scgm.lh,
                              aal90$frontal.rh, aal90$parietal.rh, aal90$temporal.rh,
                              aal90$occipital.rh, aal90$insula.rh, aal90$limbic.rh,
                              aal90$scgm.rh)

    } else if (atlas == 'lpba40') {
      lobes <- vector('integer', length=Nv)
      lobes[lpba40$frontal] <- 1
      lobes[lpba40$parietal] <- 2
      lobes[lpba40$temporal] <- 3
      lobes[lpba40$occipital] <- 4
      lobes[lpba40$insula] <- 5
      lobes[lpba40$cingulate] <- 6
      lobes[lpba40$scgm] <- 7
      V(g)$lobe <- lobes

      V(g)$circle.layout <- c(lpba40$frontal.lh, lpba40$parietal.lh, lpba40$temporal.lh,
                              lpba40$occipital.lh, lpba40$insula.lh, lpba40$cingulate.lh,
                              lpba40$scgm.lh,
                              lpba40$frontal.rh, lpba40$parietal.rh, lpba40$temporal.rh,
                              lpba40$occipital.rh, lpba40$insula.rh, lpba40$cingulate.rh,
                              lpba40$scgm.rh)
    } else if (atlas == 'hoa112') {
      lobes <- vector('integer', length=Nv)
      lobes[hoa112$frontal] <- 1
      lobes[hoa112$parietal] <- 2
      lobes[hoa112$temporal] <- 3
      lobes[hoa112$occipital] <- 4
      lobes[hoa112$insula] <- 5
      lobes[hoa112$cingulate] <- 6
      lobes[hoa112$scgm] <- 7
      V(g)$lobe <- lobes

      V(g)$circle.layout <- c(hoa112$frontal.lh, hoa112$parietal.lh, hoa112$temporal.lh,
                              hoa112$occipital.lh, hoa112$insula.lh, hoa112$cingulate.lh,
                              hoa112$scgm.lh,
                              hoa112$frontal.rh, hoa112$parietal.rh, hoa112$temporal.rh,
                              hoa112$occipital.rh, hoa112$insula.rh, hoa112$cingulate.rh,
                              hoa112$scgm.rh)
    } else if (atlas == 'brainsuite') {
      lobes <- vector('integer', length=Nv)
      lobes[brainsuite$frontal] <- 1
      lobes[brainsuite$parietal] <- 2
      lobes[brainsuite$temporal] <- 3
      lobes[brainsuite$occipital] <- 4
      lobes[brainsuite$insula] <- 5
      lobes[brainsuite$cingulate] <- 6
      lobes[brainsuite$scgm] <- 7
      V(g)$lobe <- lobes

      V(g)$circle.layout <- c(brainsuite$frontal.lh, brainsuite$parietal.lh,
                              brainsuite$temporal.lh, brainsuite$occipital.lh,
                              brainsuite$insula.lh, brainsuite$cingulate.lh,
                              brainsuite$scgm.lh,
                              brainsuite$frontal.rh, brainsuite$parietal.rh,
                              brainsuite$temporal.rh, brainsuite$occipital.rh,
                              brainsuite$insula.rh, brainsuite$cingulate.rh,
                              brainsuite$scgm.rh)
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
