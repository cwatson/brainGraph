#' Set graph, vertex, and edge attributes common in MRI analyses
#'
#' This function sets a number of graph, vertex, and edge attributes for a
#' given \code{igraph} graph object. These are all measures that are common in
#' MRI analyses of brain networks.
#'
#' @param g An \code{igraph} graph object
#' @param atlas Character vector indicating which atlas was used (default:
#'   \code{NULL})
#' @param modality Character vector indicating imaging modality (e.g. 'dti')
#'   (default: \code{NULL})
#' @param subject Character vector indicating subject ID (default: \code{NULL})
#' @param group Character vector indicating group membership (default: NULL)
#' @param rand Logical indicating if the graph is random or not (default: FALSE)
#' @param use.parallel Logical indicating whether or not to use \emph{foreach}
#'   (default: TRUE)
#' @param A Numeric matrix; the (weighted) adjacency matrix, which can be used
#'   for faster calculation of local efficiency (default: \code{NULL})
#' @export
#'
#' @return g An \code{igraph} graph object with the following attributes:
#' \item{Graph-level}{Package version, atlas, density, connected component sizes,
#' diameter, \# of triangles, transitivity, average path length, assortativity,
#' clique number, global & local efficiency, modularity, vulnerability, hub score,
#' rich-club coefficient, \# of hubs, edge asymmetry, and modality}
#' \item{Vertex-level}{Degree, strength, betweenness/eigenvector and leverage
#' centralities, hubs, transitivity (local), coreness, local & nodal efficiency,
#' color (community), color (lobe), color (component), membership (community),
#' membership (component), gateway and participation coefficients, within-module
#' degree z-score, vulnerability, and coordinates (x, y, and z)}
#' \item{Edge-level}{Color (community), color (lobe), color (component), edge
#' betweenness, Euclidean distance (in mm)}
#'
#' @seealso \code{\link[igraph]{components}, \link[igraph]{diameter},
#' \link[igraph]{clique_num}, \link[igraph]{centr_betw}, \link{part_coeff},
#' \link[igraph]{edge.betweenness}, \link[igraph]{centr_eigen},
#' \link{gateway_coeff}, \link[igraph]{hub.score},
#' \link[igraph]{authority.score}, \link[igraph]{transitivity},
#' \link[igraph]{mean_distance}, \link[igraph]{assortativity.degree},
#' \link[igraph]{cluster_louvain}, \link{efficiency},
#' \link{set_edge_color}, \link{rich_club_coeff},
#' \link{within_module_deg_z_score}, \link[igraph]{coreness},
#' \link{edge_spatial_dist}, \link{vulnerability}, \link{centr_lev},
#' \link{edge_asymmetry}, \link[igraph]{graph.knn}, \link{vertex_spatial_dist}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

set_brainGraph_attr <- function(g, atlas=NULL, modality=NULL,
                                subject=NULL, group=NULL, rand=FALSE,
                                use.parallel=TRUE, A=NULL) {
  name <- NULL
  stopifnot(is_igraph(g))

  g$version <- packageVersion('brainGraph')
  if (!'degree' %in% vertex_attr_names(g)) V(g)$degree <- degree(g)
  g$Cp <- transitivity(g, type='localaverage')
  g$Lp <- mean_distance(g)
  # Get the rich club coeff for all possible degree values
  R <- lapply(1:max(V(g)$degree), function(x) rich_club_coeff(g, x))
  phi <- vapply(R, with, numeric(1), phi)
  Nk <- vapply(R, with, numeric(1), Nk)
  Ek <- vapply(R, with, numeric(1), Ek)
  g$rich <- data.frame(phi=round(phi, 4), Nk=Nk, Ek=Ek)
  g$E.global <- efficiency(g, 'global', weights=NA)
  comm <- cluster_louvain(g, weights=NA)
  g$mod <- max(comm$modularity)

  if (!isTRUE(rand)) {
    if (!is.null(group)) g$Group <- group
    if (!is.null(subject)) g$name <- subject
    if (!is.null(modality)) g$modality <- modality

    # Graph-level attributes
    #-----------------------------------------------------------------------------
    g$density <- round(graph.density(g), digits=3)

    # Connected components
    clusts <- components(g)
    comps <- rev(table(clusts$csize))
    g$conn.comp <- data.frame(size=as.integer(names(comps)),
                              number=as.integer(comps))
    g$max.comp <- g$conn.comp[1, 1]
    g$clique.num <- clique_num(g)
    g$num.tri <- sum(count_triangles(g)) / 3
    g$diameter <- diameter(g, weights=NA)
    g$transitivity <- transitivity(g)
    g$assortativity <- assortativity_degree(g)

    if (is.weighted(g)) {
      V(g)$strength <- graph.strength(g)
      g$strength <- mean(V(g)$strength)
      V(g)$knn.wt <- graph.knn(g)$knn
      R <- lapply(1:max(V(g)$degree),
                  function(x) rich_club_coeff(g, x, weighted=TRUE))
      phi <- vapply(R, with, numeric(1), phi)
      Nk <- vapply(R, with, numeric(1), Nk)
      Ek <- vapply(R, with, numeric(1), Ek)
      g$rich.wt <- data.frame(phi=round(phi, 4), Nk=Nk, Ek=Ek)
      comm.wt <- cluster_louvain(g)
      g$mod.wt <- max(comm.wt$modularity)
      x <- comm.wt$membership
      V(g)$comm.wt <- match(x, order(table(x), decreasing=TRUE))
      V(g)$color.comm.wt <- set_vertex_color(V(g)$comm.wt)[V(g)$comm.wt]
      E(g)$color.comm.wt <- set_edge_color(g, V(g)$comm.wt)
      V(g)$GC.wt <- gateway_coeff(g, V(g)$comm.wt)
      V(g)$PC.wt <- part_coeff(g, V(g)$comm.wt)
      V(g)$z.score.wt <- within_module_deg_z_score(g, V(g)$comm.wt)
      V(g)$transitivity.wt <- transitivity(g, type='weighted')

      # Need to convert weights for distance measures
      E(g)$weight <- 1 / E(g)$weight
      V(g)$E.local.wt <- efficiency(g, type='local',
                                          use.parallel=use.parallel, A=A)
      g$E.local.wt <- mean(V(g)$E.local.wt)
      V(g)$E.nodal.wt <- efficiency(g, 'nodal')
      g$E.global.wt <- mean(V(g)$E.nodal.wt)
      g$diameter.wt <- diameter(g)
      Lpv.wt <- distances(g)
      Lpv.wt[is.infinite(Lpv.wt)] <- NA
      V(g)$Lp.wt <- rowMeans(Lpv.wt, na.rm=TRUE)

      # Convert back to connection strength
      E(g)$weight <- 1 / E(g)$weight
    }

    if (is_directed(g)) {
      hubs <- hub_score(g)
      g$hub.score <- hubs$value
      authorities <- authority_score(g)
      g$authority.score <- authorities$value
      V(g)$hub.score <- hubs$vector
      V(g)$authority.score <- authorities$vector
    }

    # Vertex-level attributes
    #-----------------------------------------------------------------------------
    # 'lobe', 'hemi', 'lobe.hemi' attributes, and colors for each lobe
    if (!is.null(atlas)) {
      g$atlas <- atlas
      atlas.dt <- eval(parse(text=atlas))
      if (!is_named(g)) V(g)$name <- atlas.dt[, name]

      g <- assign_lobes(g)
      g$assortativity.lobe <- assortativity_nominal(g, as.integer(factor(V(g)$lobe)))
      g$assortativity.lobe.hemi <- assortativity_nominal(g, V(g)$lobe.hemi)

      g$asymm <- edge_asymmetry(g)$asymm
      V(g)$asymm <- edge_asymmetry(g, 'vertex')$asymm

      E(g)$dist <- edge_spatial_dist(g)
      g$spatial.dist <- mean(E(g)$dist)
      V(g)$dist <- vertex_spatial_dist(g)
      V(g)$dist.strength <- V(g)$dist * V(g)$degree

      if (atlas %in% c('destrieux', 'destrieux.scgm')) {
        g$assortativity.class <- assortativity_nominal(g, V(g)$class)
      }
      if (atlas == 'dosenbach160') {
        g$assortativity.network <- assortativity_nominal(g, as.integer(factor(V(g)$network)))
      }
    }

    V(g)$knn <- graph.knn(g, weights=NA)$knn

    Lpv <- distances(g, weights=NA)
    Lpv[is.infinite(Lpv)] <- NA
    V(g)$Lp <- rowMeans(Lpv, na.rm=TRUE)

    E(g)$btwn <- edge.betweenness(g)
    V(g)$btwn.cent <- centr_betw(g)$res
    V(g)$hubs <- 0  # I define hubs as vertices w/ btwn.cent > mean + sd
    V(g)$hubs[which(V(g)$btwn.cent > mean(V(g)$btwn.cent) + sd(V(g)$btwn.cent))] <- 1
    g$num.hubs <- sum(V(g)$hubs)
    V(g)$ev.cent <- centr_eigen(g)$vector
    V(g)$lev.cent <- centr_lev(g)
    V(g)$coreness <- coreness(g)
    V(g)$transitivity <- transitivity(g, type='local', isolates='zero')
    V(g)$E.local <- efficiency(g, type='local', weights=NA,
                                     use.parallel=use.parallel, A=A)
    V(g)$E.nodal <- efficiency(g, type='nodal', weights=NA)
    g$E.local <- mean(V(g)$E.local)
    V(g)$vulnerability <- vulnerability(g, use.parallel=use.parallel)
    g$vulnerability <- max(V(g)$vulnerability)
    V(g)$eccentricity <- eccentricity(g)

    # Community stuff
    x <- comm$membership
    V(g)$comm <- match(x, order(table(x), decreasing=TRUE))
    V(g)$color.comm <- set_vertex_color(V(g)$comm)[V(g)$comm]
    E(g)$color.comm <- set_edge_color(g, V(g)$comm)

    x <- clusts$membership
    V(g)$comp <- match(x, order(table(x), decreasing=TRUE))
    V(g)$color.comp <- set_vertex_color(V(g)$comp)[V(g)$comp]
    E(g)$color.comp <- set_edge_color(g, V(g)$comp)

    V(g)$circle.layout.comm <- order(V(g)$comm, V(g)$degree)

    V(g)$GC <- gateway_coeff(g, V(g)$comm)
    V(g)$PC <- part_coeff(g, V(g)$comm)
    V(g)$z.score <- within_module_deg_z_score(g, V(g)$comm)
  }

  g
}

#' @inheritParams set_brainGraph_attr
#' @export
#' @rdname set_brainGraph_attr

set.brainGraph.attributes <- function(g, atlas=NULL, modality=NULL,
                                      subject=NULL, group=NULL, rand=FALSE,
                                      use.parallel=TRUE) {
  .Deprecated('set_brainGraph_attr')
  set_brainGraph_attr(g, atlas, modality, subject, group, rand, use.parallel)
}
