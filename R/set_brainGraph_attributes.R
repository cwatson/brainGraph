#' Set graph, vertex, and edge attributes common in MRI analyses
#'
#' This function sets a number of graph, vertex, and edge attributes for a
#' given \code{igraph} graph object. These are all measures that are common in
#' MRI analyses of brain networks.
#'
#' \code{xfm.type} allows you to choose from 3 options for transforming edge
#' weights when calculating distance-based metrics (e.g., shortest paths). There
#' is no "best-practice" for choosing one over the other, but the reciprocal is
#' probably most common.
#' \itemize{
#'   \item \code{1/w}: reciprocal (default)
#'   \item \code{-log(w)}: the negative (natural) logarithm
#'   \item \code{1-w}: subtract weights from 1
#' }
#'
#' \code{clust.method} allows you to choose from any of the clustering
#' (community detection) functions available in \code{igraph}. These functions
#' all begin with \code{clust_}; the function argument should not include this
#' leading character string. The default value is \code{louvain}, which calls
#' \code{\link[igraph]{cluster_louvain}}. If there are any negative edge
#' weights, and the selected method is anything other than \code{spinglass} or
#' \code{walktrap}, then \code{walktrap} is used (calling
#' \code{\link[igraph]{cluster_walktrap}}). If \code{edge_betweenness} is
#' selected and the graph is weighted, then the edges are first transformed (via
#' \code{\link{xfm.weights}}), because the algorithm considers edges as
#' \emph{distances}.
#'
#' Since \code{v2.4.0}, hubs are calculated by the new function
#' \code{\link{hubness}}. It is calculated using edge weights in addition to the
#' unweighted version of the graph.
#'
#' @param g An \code{igraph} graph object
#' @param atlas Character vector indicating which atlas was used (default:
#'   \code{NULL})
#' @param rand Logical indicating if the graph is random or not (default:
#'   \code{FALSE})
#' @param use.parallel Logical indicating whether or not to use \emph{foreach}
#'   (default: \code{TRUE})
#' @param A Numeric matrix; the (weighted) adjacency matrix, which can be used
#'   for faster calculation of local efficiency (default: \code{NULL})
#' @param xfm.type Character string indicating how to transform edge weights
#'   (default: \code{1/w} [reciprocal])
#' @param clust.method Character string indicating which method to use for
#'   community detection. Default: \code{'louvain'}
#' @param ... Other arguments passed to \code{\link{make_brainGraph}}
#' @export
#'
#' @return g An \code{igraph} graph object with the following attributes:
#'   \item{Graph-level}{Density, connected component sizes, diameter, \# of
#'     triangles, transitivity, average path length, assortativity, global &
#'     local efficiency, modularity, vulnerability, hub score, rich-club
#'     coefficient, \# of hubs, edge asymmetry, and modality}
#'   \item{Vertex-level}{Degree, strength; betweenness, eigenvector, and
#'     leverage centralities; hubs; transitivity (local); k-core, s-core; local
#'     & nodal efficiency; color (community, lobe, component); membership
#'     (community, lobe, component); gateway and participation coefficients,
#'     within-module degree z-score; vulnerability; and coordinates (x, y, and
#'     z)}
#'   \item{Edge-level}{Color (community, lobe, component), edge betweenness,
#'     Euclidean distance (in mm), weight (if weighted)}
#'
#' @seealso \code{\link[igraph]{components}}, \code{\link[igraph]{diameter}},
#' \code{\link[igraph]{clique_num}}, \code{\link[igraph]{centr_betw}},
#' \code{\link{part_coeff}}, \code{\link[igraph]{edge.betweenness}},
#' \code{\link[igraph]{centr_eigen}}, \code{\link{gateway_coeff}},
#' \code{\link[igraph]{transitivity}}, \code{\link[igraph]{mean_distance}},
#' \code{\link[igraph]{assortativity_degree}}, \code{\link{efficiency}},
#' \code{\link[igraph]{assortativity_nominal}}, \code{\link[igraph]{coreness}},
#' \code{\link[igraph]{communities}}, \code{\link{set_edge_color}},
#' \code{\link{rich_club_coeff}}, \code{\link{s_core}}, \code{\link{centr_lev}},
#' \code{\link{within_module_deg_z_score}}, \code{\link{edge_spatial_dist}},
#' \code{\link{vulnerability}}, \code{\link{edge_asymmetry}},
#' \code{\link[igraph]{graph.knn}}, \code{\link{vertex_spatial_dist}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

set_brainGraph_attr <- function(g, atlas=NULL, rand=FALSE, use.parallel=TRUE, A=NULL,
                                xfm.type=c('1/w', '-log(w)', '1-w'),
                                clust.method='louvain', ...) {
  name <- NULL
  stopifnot(is_igraph(g))
  clust.funs <- ls('package:igraph')[grep('cluster_', ls('package:igraph'))]
  clust.funs <- gsub('cluster_', '', clust.funs)
  if (!clust.method %in% clust.funs) {
    stop('Invalid clustering method! You must choose from the following:\n', paste(clust.funs, collapse='\n'))
  }

  if (!'degree' %in% vertex_attr_names(g)) V(g)$degree <- degree(g)
  g$Cp <- transitivity(g, type='localaverage')
  g$Lp <- mean_distance(g)
  g$rich <- rich_club_all(g)
  g$E.global <- efficiency(g, 'global', weights=NA)

  # Handle different cases for different community detection methods
  if (clust.method == 'spinglass' & !is.connected(g)) {
    warning('Invalid clustering method for an unconnected graph; using "louvain".')
    clust.method <- 'louvain'
  }
  g$clust.method <- clust.method
  if (clust.method %in% c('edge_betweenness', 'fast_greedy', 'walktrap')) {
    comm <- eval(parse(text=paste0('cluster_', clust.method, '(g, weights=NULL)')))
  } else if (clust.method == 'infomap') {
    comm <- cluster_infomap(g, e.weights=NULL)
  } else {
    comm <- eval(parse(text=paste0('cluster_', clust.method, '(g, weights=NA)')))
  }
  g$mod <- max(comm$modularity)

  if (!isTRUE(rand)) {
    # Graph-level attributes
    #-----------------------------------------------------------------------------
    g$density <- graph.density(g)

    # Connected components
    clusts <- components(g)
    comps <- rev(table(clusts$csize))
    g$conn.comp <- data.frame(size=as.integer(names(comps)),
                              number=as.integer(comps))
    g$max.comp <- g$conn.comp[1, 1]
    g$num.tri <- sum(count_triangles(g)) / 3
    g$diameter <- diameter(g, weights=NA)
    g$transitivity <- transitivity(g)
    g$assortativity <- assortativity_degree(g)

    if (is_weighted(g)) {
      xfm.type <- match.arg(xfm.type)
      # Handle different cases for different community detection methods
      if (any(E(g)$weight < 0) & !clust.method %in% c('spinglass', 'walktrap')) {
        warning('Invalid clustering method for negative edge weights; using "walktrap".')
        clust.method <- 'walktrap'
      }
      g$clust.method.wt <- clust.method
      if (clust.method == 'edge_betweenness') {
        g <- xfm.weights(g, xfm.type)
        comm.wt <- cluster_edge_betweenness(g)
        g <- xfm.weights(g, xfm.type, invert=TRUE)
      } else {
        comm.wt <- eval(parse(text=paste0('cluster_', clust.method, '(g)')))
      }
      V(g)$strength <- graph.strength(g)
      g$strength <- mean(V(g)$strength)
      V(g)$knn.wt <- graph.knn(g)$knn
      V(g)$s.core <- s_core(g, A)
      g$rich.wt <- rich_club_all(g, weighted=TRUE)
      g$mod.wt <- max(comm.wt$modularity)
      x <- comm.wt$membership
      V(g)$comm.wt <- match(x, order(table(x), decreasing=TRUE))
      g <- set_vertex_color(g, 'color.comm.wt', V(g)$comm.wt)
      g <- set_edge_color(g, 'color.comm.wt', V(g)$comm.wt)
      V(g)$GC.wt <- gateway_coeff(g, V(g)$comm.wt)
      V(g)$PC.wt <- part_coeff(g, V(g)$comm.wt)
      V(g)$z.score.wt <- within_module_deg_z_score(g, V(g)$comm.wt)
      V(g)$transitivity.wt <- transitivity(g, type='weighted')

      # Need to convert weights for distance measures
      g <- xfm.weights(g, xfm.type)
      V(g)$E.local.wt <- efficiency(g, type='local',
                                    use.parallel=use.parallel, A=A)
      g$E.local.wt <- mean(V(g)$E.local.wt)
      V(g)$E.nodal.wt <- efficiency(g, 'nodal')
      g$E.global.wt <- mean(V(g)$E.nodal.wt)
      g$diameter.wt <- diameter(g)
      Lpv.wt <- distances(g)
      Lpv.wt[is.infinite(Lpv.wt)] <- NA
      V(g)$Lp.wt <- rowMeans(Lpv.wt, na.rm=TRUE)
      g$Lp.wt <- mean(Lpv.wt[upper.tri(Lpv.wt)], na.rm=T)

      # Convert back to connection strength
      g <- xfm.weights(g, xfm.type, invert=TRUE)
      V(g)$hubs.wt <- hubness(g)
      g$num.hubs.wt <- sum(V(g)$hubs.wt >= 2)
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
      atlas.dt <- get(atlas)
      if (!is_named(g)) V(g)$name <- atlas.dt[, name]

      g <- make_brainGraph(g, atlas, ...)
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
    V(g)$hubs <- hubness(g, weights=NA)
    g$num.hubs <- sum(V(g)$hubs >= 2)
    V(g)$ev.cent <- centr_eigen(g)$vector
    V(g)$lev.cent <- centr_lev(g)
    V(g)$k.core <- coreness(g)
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
    g <- set_vertex_color(g, 'color.comm', V(g)$comm)
    g <- set_edge_color(g, 'color.comm', V(g)$comm)

    x <- clusts$membership
    V(g)$comp <- match(x, order(table(x), decreasing=TRUE))
    g <- set_vertex_color(g, 'color.comp', V(g)$comp)
    g <- set_edge_color(g, 'color.comp', V(g)$comp)

    V(g)$circle.layout.comm <- order(V(g)$comm, V(g)$degree)

    V(g)$GC <- gateway_coeff(g, V(g)$comm)
    V(g)$PC <- part_coeff(g, V(g)$comm)
    V(g)$z.score <- within_module_deg_z_score(g, V(g)$comm)
  }

  return(g)
}
