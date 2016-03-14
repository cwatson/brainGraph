#' Set a number of graph and vertex attributes useful in MRI analyses
#'
#' This function will set a number of graph, vertex, and edge attributes of a
#' given igraph object.
#'
#' @param g An igraph object
#' @param atlas A character vector indicating which atlas was used for the nodes
#' @param modality A character vector indicating imaging modality (e.g. 'dti')
#' @param subject A character vector indicating subject ID (default: NULL)
#' @param group A character vector indicating group membership (default: NULL)
#' @param rand Logical indicating if the graph is random or not (default: FALSE)
#' @export
#'
#' @return g A copy of the same graph, with the following attributes:
#' \item{Graph-level}{Package version, atlas, density, connected component sizes,
#' diameter, \# of triangles, transitivity, average path length, assortativity,
#' clique number, global & local efficiency, modularity, vulnerability, hub score,
#' rich-club coefficient, \# of hubs, edge asymmetry, and modality}
#' \item{Vertex-level}{Degree, strength, betweenness/eigenvector and leverage
#' centralities, hubs, transitivity (local), coreness, local & nodal efficiency,
#' color (community), color (lobe), color (component), membership (community),
#' membership (component), participation coefficient, within-module degree
#' z-score, vulnerability, and coordinates (x, y, and z)}
#' \item{Edge-level}{Color (community), color (lobe), color (component), edge
#' betweenness, Euclidean distance (in mm)}
#'
#' @seealso \code{\link[igraph]{components}, \link[igraph]{diameter},
#' \link[igraph]{clique_num}, \link[igraph]{centr_betw}, \link{part.coeff},
#' \link[igraph]{edge.betweenness}, \link[igraph]{centr_eigen},
#' \link[igraph]{hub.score},
#' \link[igraph]{authority.score}, \link[igraph]{transitivity},
#' \link[igraph]{mean_distance}, \link[igraph]{assortativity.degree},
#' \link[igraph]{cluster_louvain}, \link{graph.efficiency}, \link{color.edges},
#' \link{rich.club.coeff}, \link{within_module_deg_z_score},
#' \link[igraph]{coreness}, \link{edge_spatial_dist}, \link{vulnerability},
#' \link{centr_lev}, \link{edge_asymmetry}, \link[igraph]{graph.knn},
#' \link{vertex_spatial_dist}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

set.brainGraph.attributes <- function(g, atlas=NULL, modality=NULL,
                                      subject=NULL, group=NULL, rand=FALSE) {
  name <- NULL
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }

  g$version <- packageVersion('brainGraph')
  if (!'degree' %in% graph_attr_names(g)) V(g)$degree <- degree(g)
  g$Cp <- transitivity(g, type='localaverage')
  g$Lp <- mean_distance(g)
  # Get the rich club coeff for all possible degree values
  R <- lapply(1:max(V(g)$degree), function(x) rich.club.coeff(g, x))
  phi <- vapply(R, with, numeric(1), phi)
  Nk <- vapply(R, with, numeric(1), Nk)
  Ek <- vapply(R, with, numeric(1), Ek)
  g$rich <- data.frame(phi=round(phi, 4), Nk=Nk, Ek=Ek)
  g$E.global <- graph.efficiency(g, 'global', weights=NA)
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
                              number=unname(comps))
    g$max.comp <- g$conn.comp[1, 1]
    g$clique.num <- clique_num(g)
    g$num.tri <- sum(count_triangles(g)) / 3
    g$diameter <- diameter(g, weights=NA)
    g$transitivity <- transitivity(g)
    g$assortativity <- assortativity_degree(g)

    if (is.weighted(g)) {
      V(g)$strength <- graph.strength(g)
      V(g)$knn.wt <- graph.knn(g)$knn
      V(g)$E.local.wt <- graph.efficiency(g, type='local')
      g$E.local.wt <- mean(V(g)$E.local.wt)
      V(g)$E.nodal.wt <- graph.efficiency(g, 'nodal')
      g$E.global.wt <- mean(V(g)$E.nodal.wt)
      g$diameter.wt <- diameter(g)
      R <- lapply(1:max(V(g)$degree),
                  function(x) rich.club.coeff(g, x, weighted=TRUE))
      phi <- vapply(R, with, numeric(1), phi)
      Nk <- vapply(R, with, numeric(1), Nk)
      Ek <- vapply(R, with, numeric(1), Ek)
      g$rich.wt <- data.frame(phi=round(phi, 4), Nk=Nk, Ek=Ek)
      comm.wt <- cluster_louvain(g)
      g$mod.wt <- max(comm.wt$modularity)

      x <- comm.wt$membership
      V(g)$comm.wt <- match(x, order(table(x), decreasing=TRUE))
      V(g)$color.comm.wt <- color.vertices(V(g)$comm.wt)[V(g)$comm.wt]
      E(g)$color.comm.wt <- color.edges(g, V(g)$comm.wt)

      V(g)$PC.wt <- part.coeff(g, V(g)$comm.wt)
      V(g)$z.score.wt <- within_module_deg_z_score(g, V(g)$comm.wt)

      V(g)$transitivity.wt <- transitivity(g, type='weighted')
      Lpv.wt <- distances(g)
      Lpv.wt[is.infinite(Lpv.wt)] <- NA
      V(g)$Lp.wt <- rowMeans(Lpv.wt, na.rm=TRUE)
    }

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
    # 'lobe', 'hemi', 'lobe.hemi' attributes, and colors for each lobe
    if (!is.null(atlas)) {
      g$atlas <- atlas
      atlas.dt <- eval(parse(text=atlas))
      if (!'name' %in% vertex_attr_names(g)) V(g)$name <- atlas.dt[, name]

      g <- assign_lobes(g)
      V(g)$color.lobe <- group.cols[V(g)$lobe]
      E(g)$color.lobe <- color.edges(g, V(g)$lobe)
      g$assortativity.lobe <- assortativity_nominal(g, V(g)$lobe)
      g$assortativity.lobe.hemi <- assortativity_nominal(g, V(g)$lobe.hemi)

      g$asymm <- edge_asymmetry(g)$asymm
      V(g)$asymm <- edge_asymmetry(g, 'vertex')$asymm

      if (atlas == 'destrieux') {
        V(g)$color.class <- group.cols[V(g)$class]
        g$assortativity.class <- assortativity_nominal(g, V(g)$class)
        E(g)$color.class <- color.edges(g, V(g)$class)
      }

      # Add the spatial coordinates for plotting over the brain
      if ('name' %in% vertex_attr_names(g)) {
        x <- y <- z <- x.mni <- y.mni <- z.mni <- NULL
        vorder <- match(V(g)$name, atlas.dt$name)
        V(g)$x <- atlas.dt[vorder, x]
        V(g)$y <- atlas.dt[vorder, y]
        V(g)$z <- atlas.dt[vorder, z]
        V(g)$x.mni <- atlas.dt[vorder, x.mni]
        V(g)$y.mni <- atlas.dt[vorder, y.mni]
        V(g)$z.mni <- atlas.dt[vorder, z.mni]
      }
    }

    V(g)$knn <- graph.knn(g, weights=NA)$knn
    V(g)$btwn.cent <- centr_betw(g)$res
    V(g)$hubs <- 0  # I define hubs as vertices w/ btwn.cent > mean + sd
    V(g)$hubs[which(V(g)$btwn.cent > mean(V(g)$btwn.cent) + sd(V(g)$btwn.cent))] <- 1
    g$num.hubs <- sum(V(g)$hubs)
    Lpv <- distances(g, weights=NA)
    Lpv[is.infinite(Lpv)] <- NA
    V(g)$Lp <- rowMeans(Lpv, na.rm=TRUE)
    V(g)$ev.cent <- centr_eigen(g)$vector
    V(g)$lev.cent <- centr_lev(g)
    V(g)$coreness <- coreness(g)
    V(g)$transitivity <- transitivity(g, type='local', isolates='zero')
    V(g)$E.local <- graph.efficiency(g, type='local', weights=NA)
    V(g)$E.nodal <- graph.efficiency(g, type='nodal', weights=NA)
    g$E.local <- mean(V(g)$E.local)
    V(g)$vulnerability <- vulnerability(g)
    g$vulnerability <- max(V(g)$vulnerability)
    V(g)$eccentricity <- eccentricity(g)

    # Community stuff
    x <- comm$membership
    V(g)$comm <- match(x, order(table(x), decreasing=TRUE))
    V(g)$color.comm <- color.vertices(V(g)$comm)[V(g)$comm]
    E(g)$color.comm <- color.edges(g, V(g)$comm)

    x <- clusts$membership
    V(g)$comp <- match(x, order(table(x), decreasing=TRUE))
    V(g)$color.comp <- color.vertices(V(g)$comp)[V(g)$comp]
    E(g)$color.comp <- color.edges(g, V(g)$comp)

    V(g)$circle.layout.comm <- order(V(g)$comm, V(g)$degree)

    V(g)$PC <- part.coeff(g, V(g)$comm)
    V(g)$z.score <- within_module_deg_z_score(g, V(g)$comm)

    # Edge attributes
    #-----------------------------------------------------------------------------
    E(g)$btwn <- edge.betweenness(g)
    E(g)$dist <- edge_spatial_dist(g)

    g$spatial.dist <- mean(E(g)$dist)
    V(g)$dist <- vertex_spatial_dist(g)
    V(g)$dist.strength <- V(g)$dist * V(g)$degree
  }

  g
}
