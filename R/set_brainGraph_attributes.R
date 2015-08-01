#' Set a number of graph and vertex attributes useful in MRI analyses
#'
#' This function will set a number of graph, vertex, and edge attributes of a
#' given igraph object.
#'
#' @param g An igraph object
#' @param atlas A character vector indicating which atlas was used for the nodes
#' @param rand Logical indicating if the graph is random or not (default: FALSE)
#'
#' @export
#'
#' @return g A copy of the same graph, with the following attributes:
#' \item{Graph-level}{Package version, atlas, density, connected component sizes,
#' diameter, \# of triangles, transitivity, average path length, assortativity,
#' clique number, global & local efficiency, modularity, vulnerability, hub score,
#' rich-club coefficient, \# of hubs, and edge asymmetry}
#' \item{Vertex-level}{Degree, strength, betweenness/eigenvector/subgraph and
#' leverage centralities, hubs, transitivity (local), coreness, local & nodal
#' efficiency, color (community), color (lobe), color (component), membership
#' (community), membership (component), participation coefficient, within-module
#' degree z-score, vulnerability, and coordinates (x, y, and z)}
#' \item{Edge-level}{Color (community), color (lobe), color (component), edge
#' betweenness, Euclidean distance (in mm)}
#'
#' @seealso \code{\link[igraph]{components}, \link[igraph]{diameter},
#' \link[igraph]{clique_num}, \link[igraph]{centr_betw}, \link{part.coeff},
#' \link[igraph]{edge.betweenness}, \link[igraph]{centr_eigen},
#' \link[igraph]{subgraph.centrality}, \link[igraph]{hub.score},
#' \link[igraph]{authority.score}, \link[igraph]{transitivity},
#' \link[igraph]{mean_distance}, \link[igraph]{assortativity.degree},
#' \link[igraph]{cluster_louvain}, \link{graph.efficiency}, \link{color.edges},
#' \link{rich.club.coeff}, \link{within_module_deg_z_score},
#' \link[igraph]{coreness}, \link{spatial.dist}, \link{vulnerability},
#' \link{centr_lev}, \link{edge_asymmetry}, \link[igraph]{graph.knn}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

set.brainGraph.attributes <- function(g, atlas=NULL, rand=FALSE) {

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
  g$E.global <- graph.efficiency(g, 'global')
  comm <- cluster_louvain(g)
  g$mod <- max(comm$modularity)

  if (!isTRUE(rand)) {

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
    g$diameter <- diameter(g)
    g$transitivity <- transitivity(g)
    g$assortativity <- assortativity.degree(g)

    if (is.weighted(g)) {
      V(g)$strength <- graph.strength(g)
      R <- lapply(1:max(V(g)$degree),
                  function(x) rich.club.coeff(g, x, weighted=T))
      phi <- vapply(R, with, numeric(1), phi)
      Nk <- vapply(R, with, numeric(1), Nk)
      Ek <- vapply(R, with, numeric(1), Ek)
      g$rich.wt <- data.frame(phi=round(phi, 4), Nk=Nk, Ek=Ek)
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
      atlas.dt <- eval(parse(text=data(list=atlas)))

      g <- assign_lobes(g, atlas.dt)
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

    V(g)$knn <- graph.knn(g)$knn
    V(g)$btwn.cent <- centr_betw(g)$res
    V(g)$hubs <- 0  # I define hubs as vertices w/ btwn.cent > mean + sd
    V(g)$hubs[which(V(g)$btwn.cent > mean(V(g)$btwn.cent) + sd(V(g)$btwn.cent))] <- 1
    g$num.hubs <- sum(V(g)$hubs)
    V(g)$ev.cent <- centr_eigen(g)$vector
    V(g)$subgraph.cent <- subgraph.centrality(g)
    V(g)$lev.cent <- centr_lev(g)
    V(g)$coreness <- graph.coreness(g)
    V(g)$transitivity <- transitivity(g, type='local', isolates='zero')
    V(g)$E.local <- graph.efficiency(g, type='local')
    V(g)$E.nodal <- graph.efficiency(g, type='nodal')
    g$E.local <- mean(V(g)$E.local)
    V(g)$vulnerability <- vulnerability(g)
    g$vulnerability <- max(V(g)$vulnerability)

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
    E(g)$dist <- spatial.dist(g)
  }

  g
}
