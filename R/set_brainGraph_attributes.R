#' Set graph, vertex, and edge attributes common in MRI analyses
#'
#' \code{set_brainGraph_attr} sets a number of graph, vertex, and edge
#' attributes for a given graph object. Specifically, it calculates measures
#' that are common in MRI analyses of brain networks.
#'
#' Including \code{type='random'} in the function call will reduce the number of
#' attributes calculated.
#'
#' @section Community detection:
#' \code{clust.method} allows you to choose from any of the clustering
#' (community detection) functions available in \code{igraph}. These functions
#' all begin with \code{clust_}; the function argument should not include this
#' leading character string. There are a few possibilities, depending on the
#' value and the type of input graph:
#' \enumerate{
#'   \item By default, \code{louvain} is used, calling
#'     \code{\link[igraph]{cluster_louvain}}
#'   \item Uses \code{spinglass} if there are any negative edges and/or the
#'     selected method is \code{spinglass}
#'   \item Uses \code{walktrap} if there are any negative edge weights and any
#'     other method (besided \code{spinglass}) is selected
#'   \item Automatically transforms the edge weights if \code{edge_betweenness}
#'     is selected and the graph is weighted, because the algorithm considers
#'     edges as \emph{distances}
#' }
#'
#' @section Transforming edge weights:
#' For distance-based measures, it is important to transform the edge weights so
#' that the \emph{strongest} connections are re-mapped to having the
#' \emph{lowest} weights. Then you may calculate e.g., the \emph{shortest path
#' length} which will include the strongest connections.
#'
#' \code{xfm.type} allows you to choose from 5 options for transforming edge
#' weights when calculating distance-based metrics (e.g., shortest paths). There
#' is no \dQuote{best-practice} for choosing one over the other, but the reciprocal is
#' probably most common.
#' \describe{
#'   \item{\code{1/w}}{reciprocal (default)}
#'   \item{\code{-log(w)}}{the negative (natural) logarithm}
#'   \item{\code{1-w}}{subtract weights from 1}
#'   \item{\code{-log10(w/max(w))}}{negative (base-10) log of normalized
#'   weights}
#'   \item{\code{-log10(w/max(w)+1)}}{same as above, but add 1 before taking
#'     the log}
#' }
#'
#' To transform the weights back to original values, specify \code{invert=TRUE}.
#'
#' @param g A graph object
#' @param use.parallel Logical indicating whether to use \emph{foreach}.
#'   Default: \code{TRUE}
#' @param A Numeric matrix; the (weighted) adjacency matrix, which can be used
#'   for faster calculation of local efficiency. Default: \code{NULL}
#' @param clust.method Character string indicating which method to use for
#'   community detection. Default: \code{'louvain'}
#' @inheritParams CreateGraphs
#' @export
#'
#' @return A graph object with the following attributes:
#'   \item{Graph-level}{Density, connected component sizes, diameter, # of
#'     triangles, transitivity, average path length, assortativity, global &
#'     local efficiency, modularity, vulnerability, hub score, rich-club
#'     coefficient, # of hubs, edge asymmetry}
#'   \item{Vertex-level}{Degree, strength; betweenness, eigenvector, and
#'     leverage centralities; hubs; transitivity (local); k-core, s-core; local
#'     & nodal efficiency; color (community, lobe, component); membership
#'     (community, lobe, component); gateway and participation coefficients,
#'     within-module degree z-score; vulnerability; and coordinates (x, y, and
#'     z)}
#'   \item{Edge-level}{Color (community, lobe, component), edge betweenness,
#'     Euclidean distance (in mm), weight (if weighted)}
#'
#TODO: change these to *destinations* and use in "@return"
#' @seealso \code{\link[igraph]{components}}, \code{\link[igraph]{diameter}},
#' \code{\link[igraph]{centr_betw}},
#' \code{\link{VertexRoles}}, \code{\link[igraph]{edge.betweenness}},
#' \code{\link[igraph]{centr_eigen}},
#' \code{\link[igraph]{transitivity}}, \code{\link[igraph]{mean_distance}},
#' \code{\link[igraph]{assortativity_degree}}, \code{\link{efficiency}},
#' \code{\link[igraph]{assortativity_nominal}}, \code{\link[igraph]{coreness}},
#' \code{\link[igraph]{communities}},
#' \code{\link{rich_club_coeff}}, \code{\link{s_core}}, \code{\link{centr_lev}},
#' \code{\link{edge_spatial_dist}},
#' \code{\link{vulnerability}}, \code{\link{edge_asymmetry}},
#' \code{\link[igraph]{graph.knn}}, \code{\link{vertex_spatial_dist}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @name Attributes
#' @rdname attributes

set_brainGraph_attr <- function(g, type=c('observed', 'random'),
                                use.parallel=TRUE, A=NULL,
                                xfm.type=c('1/w', '-log(w)', '1-w',
                                           '-log10(w/max(w))',
                                           '-log10(w/max(w)+1)'),
                                clust.method='louvain') {
  if (!is_igraph(g) && !is.brainGraph(g)) {
    stop('Input graph must have class either "brainGraph" or "igraph"')
  }
  clust.funs <- ls('package:igraph')[grep('cluster_', ls('package:igraph'))]
  clust.funs <- gsub('cluster_', '', clust.funs)
  if (!clust.method %in% clust.funs) {
    stop('Invalid clustering method! You must choose from the following:\n',
         paste(clust.funs, collapse='\n'))
  }

  V(g)$degree <- degree(g)
  g$Cp <- transitivity(g, type='localaverage')
  g$Lp <- mean_distance(g)
  g$rich <- rich_club_all(g)
  g$E.global <- efficiency(g, 'global', weights=NA)

  # Handle different cases for different community detection methods
  if (clust.method == 'spinglass' && !is_connected(g)) {
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

  type <- match.arg(type)
  if (type == 'observed') {
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
      if (any(E(g)$weight < 0) && !clust.method %in% c('spinglass', 'walktrap')) {
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
      g <- set_graph_colors(g, 'color.comm.wt', V(g)$comm.wt)
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
    if (!is.null(g$atlas)) {
      g$assortativity.lobe <- assortativity_nominal(g, as.integer(factor(V(g)$lobe)))
      g$assortativity.lobe.hemi <- assortativity_nominal(g, V(g)$lobe.hemi)

      g$asymm <- edge_asymmetry(g)$asymm
      V(g)$asymm <- edge_asymmetry(g, 'vertex')$asymm

      E(g)$dist <- edge_spatial_dist(g)
      g$spatial.dist <- mean(E(g)$dist)
      V(g)$dist <- vertex_spatial_dist(g)
      V(g)$dist.strength <- V(g)$dist * V(g)$degree

      if (g$atlas %in% c('destrieux', 'destrieux.scgm')) {
        g$assortativity.class <- assortativity_nominal(g, V(g)$class)
      }
      if (g$atlas == 'dosenbach160') {
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
    g <- set_graph_colors(g, 'color.comm', V(g)$comm)

    x <- clusts$membership
    V(g)$comp <- match(x, order(table(x), decreasing=TRUE))
    g <- set_graph_colors(g, 'color.comp', V(g)$comp)

    V(g)$circle.layout.comm <- order(V(g)$comm, V(g)$degree)

    V(g)$GC <- gateway_coeff(g, V(g)$comm)
    V(g)$PC <- part_coeff(g, V(g)$comm)
    V(g)$z.score <- within_module_deg_z_score(g, V(g)$comm)
  }

  return(g)
}

#' Delete all attributes of a graph
#'
#' Deletes all graph-, vertex-, and edge-level attributes of an \code{igraph}
#' graph object.
#'
#' @param g An \code{igraph} graph object
#' @param keep.names Logical indicating whether to keep the \code{name} vertex
#'   attribute (default: \code{FALSE})
#'
#' @keywords internal
#' @return An \code{igraph} graph object
#' @seealso \code{\link[igraph]{delete_graph_attr},
#'   \link[igraph]{delete_vertex_attr}, \link[igraph]{delete_edge_attr}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

delete_all_attr <- function(g, keep.names=FALSE) {
  for (att in graph_attr_names(g)) g <- delete_graph_attr(g, att)
  for (att in edge_attr_names(g)) g <- delete_edge_attr(g, att)
  if (isTRUE(keep.names)) {
    vattrs <- setdiff(vertex_attr_names(g), 'name')
  } else {
    vattrs <- vertex_attr_names(g)
  }
  for (att in vattrs) g <- delete_vertex_attr(g, att)

  return(g)
}

#' Color graph vertices and edges
#'
#' \code{set_graph_colors} takes an integer vector representing membership of
#' some grouping (e.g., a community or connected component) and creates a
#' character vector of colors for each grouping. Isolated vertices will be
#' colored \emph{gray}. Edges are assigned the same color if connected to
#' vertices in the same group, and assigned \emph{gray} otherwise.
#'
#' @param g An \code{igraph} graph object
#' @param name Character string of the name of the attribute to add
#' @param memb An integer vector representing membership of e.g. a community
#'
#' @return The same graph with additional vertex and edge attributes
#' @keywords internal

set_graph_colors <- function(g, name, memb) {
  stopifnot(length(memb) == vcount(g))
  memb <- as.integer(factor(memb))
  big.groups <- which(table(memb) > 1)

  # Vertex colors
  group.cols.memb <- rep('gray', length=max(memb))
  group.cols.memb[big.groups] <- group.cols[big.groups]

  # Edge colors
  newcols <- rep('gray50', length=ecount(g))
  tmp <- vector('list', length=max(big.groups))
  for (i in big.groups) {
    x <- which(memb == i)
    tmp[[i]] <- as.vector(E(g)[x %--% x])
    if (!is.null(tmp[[i]])) newcols[tmp[[i]]] <- group.cols[i]
  }

  g <- set_vertex_attr(g, name, value=group.cols.memb[memb])
  g <- set_edge_attr(g, name, value=newcols)
  return(g)
}

#' Transform edge weights
#'
#' @param xfm.type Character string specifying how to transform the weights
#'   (default: \code{1/w})
#' @param invert Logical indicating whether or not to invert the transformation
#'   (default: \code{FALSE})
#' @export
#' @return \code{xfm.weights} returns the same graph object, with transformed
#'   edge weights plus a graph attribute (\code{xfm.type}) recording the method
#'   of transformation
#' @rdname attributes

xfm.weights <- function(g, xfm.type=c('1/w', '-log(w)', '1-w',
                                      '-log10(w/max(w))', '-log10(w/max(w)+1)'),
                        invert=FALSE) {
  stopifnot(is_igraph(g), is_weighted(g))
  xfm.type <- match.arg(xfm.type)
  if (xfm.type == '1/w') {
    E(g)$weight <- 1 / E(g)$weight

  } else if (xfm.type == '-log(w)') {
    if (isTRUE(invert)) {
      E(g)$weight <- exp(-E(g)$weight)
    } else {
      E(g)$weight <- -log(E(g)$weight)
    }

  } else if (xfm.type == '1-w') {
    E(g)$weight <- 1 - E(g)$weight

  } else if (xfm.type == '-log10(w/max(w))') {
    if (isTRUE(invert)) {
      E(g)$weight <- g$max.weight / 10^(E(g)$weight)
    } else {
      g$max.weight <- max(E(g)$weight)
      E(g)$weight <- -log10(E(g)$weight / g$max.weight)
    }

  } else if (xfm.type == '-log10(w/max(w)+1)') {
    if (isTRUE(invert)) {
      E(g)$weight <- g$max.weight * (1 / 10^(E(g)$weight) - 1)
    } else {
      g$max.weight <- max(E(g)$weight)
      E(g)$weight <- -log10(E(g)$weight / g$max.weight + 1)
    }
  }

  g$xfm.type <- xfm.type
  return(g)
}
