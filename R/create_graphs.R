################################################################################
# MAIN CREATION FUNCTIONS
################################################################################

#' Create a brainGraph object
#'
#' \code{make_brainGraph} is the main creation function for creating a
#' \code{brainGraph} graph object. This is simply an \code{igraph} graph
#' object with additional attributes (at all levels). Several of the graph-level
#' attributes serve the purpose of providing metadata on how the connectivity
#' matrices/networks were created.
#'
#' @section Graph-level attributes:
#' Graph-level attributes added are:
#' \describe{
#'   \item{version}{The R, \code{brainGraph}, and \code{igraph} package versions
#'     used to create the graph}
#'   \item{date}{The creation date, from \code{\link{as.POSIXct}}}
#'   \item{atlas}{Character string denoting the brain atlas used}
#'   \item{type}{Character string specifying whether this is an \emph{observed}
#'     or \emph{random} graph}
#'   \item{modality}{The imaging modality; you can choose anything you like,
#'     but the \code{summary.brainGraph} knows about \code{dti}, \code{fmri},
#'     \code{thickness}, \code{area}, and \code{volume}}
#'   \item{weighting}{What edge weights represent; you can choose anything you
#'     like, but \code{summary.brainGraph} knows about \code{fa}, \code{sld}
#'     (streamline density, tractography), \code{pearson}, \code{spearman},
#'     \code{kendall}, and \code{partial} (partial correlation coefficient)}
#'   \item{threshold}{Numeric indicating the threshold used to create the final
#'     connectivity matrix (if any)}
#'   \item{name}{Character string specifying the study ID or group/contrast
#'     name, depending on the \code{level} argument}
#'   \item{Group}{Character string specifying the experimental group the given
#'     subject belongs to, or if it is a group-level graph}
#' }
#'
#' @section Vertex attributes:
#' Vertex-level attributes added are:
#' \describe{
#'   \item{name}{The names of the brain regions in the network}
#'   \item{lobe}{The names of the major brain lobes for each vertex}
#'   \item{hemi}{The names of the hemisphere for each vertex}
#'   \item{lobe.hemi}{The lobe-hemisphere combination}
#'   \item{class}{The tissue class (if applicable)}
#'   \item{network}{The network (if the atlas is \code{dosenbach160})}
#'   \item{x,y,z}{The coordinates of the (centers-of-mass) brain regions in MNI
#'     space}
#'   \item{color.lobe,color.class,color.network}{Colors for vertices of their
#'     respective membership}
#'   \item{circle.layout}{Integer vector indicating the order (going
#'     counter-clockwise from the top) for circular layouts}
#' }
#'
#' @section Edge attributes:
#' Edge-level attributes added are:
#' \describe{
#'   \item{color.lobe,color.class,color.network}{Correspond(s) to the vertex
#'   attribute(s) of the same name. Inter-group edges will be colored
#'   \emph{gray}}
#' }
#'
#' @param x An \code{igraph} graph object or numeric matrix
#' @param atlas Character string specifying the brain atlas
#' @param type Character string indicating the type of graphs. Default:
#'   \code{observed}
#' @param level Character string indicating whether the graphs are subject-,
#'   group-, or contrast-specific. Default: \code{'subject'}
#' @param set.attrs Logical indicating whether to assign all graph-, vertex-,
#'   and edge-level attributes (via \code{\link{set_brainGraph_attr}}). Default:
#'   \code{TRUE}
#' @param modality Character string indicating imaging modality (e.g. 'dti').
#'   Default: \code{NULL}
#' @param weighting Character string indicating how the edges are weighted
#'   (e.g., 'fa', 'pearson', etc.). Default: \code{NULL}
#' @param threshold Numeric indicating the threshold used when sparsifying the
#'   connectivity matrix (if any). Default: \code{NULL}
#' @param name Character string indicating subject ID or group/contrast name,
#'   depending on the \code{level}. Default: \code{NULL}
#' @param Group Character string indicating group membership. Default:
#'   \code{NULL}
#' @param ... Arguments passed to \code{\link{set_brainGraph_attr}}
#' @export
#'
#' @return A \code{brainGraph} graph object with additional attributes:
#'   \item{version}{(graph) The current versions of R, \code{brainGraph}, and
#'     \code{igraph}}
#'   \item{date}{(graph)}
#'   \item{atlas}{(graph)}
#'   \item{type}{(graph)}
#'   \item{modality}{(graph)}
#'   \item{weighting}{(graph)}
#'   \item{threshold}{(graph)}
#'   \item{name}{(graph) The subject ID, group name, or contrast name (depending
#'     on the value of \code{level})}
#'   \item{Group}{(graph) only if \code{group} is specified}
#'   \item{lobe}{(vertex) Character vector of lobe names}
#'   \item{hemi}{(vertex) Character vector of hemispheres (\code{'L'},
#'     \code{'R'}, or \code{'B'})}
#'   \item{lobe.hemi}{(vertex) Integer vector indicating the lobe and
#'     hemisphere}
#'   \item{class}{(vertex) Character vector of class names (if applicable)}
#'   \item{network}{(vertex) Character vector of network names (if
#'     applicable)}
#'   \item{x, y, z, x.mni, y.mni, z.mni}{(vertex) Spatial coordinates}
#'   \item{color.lobe}{(vertex and edge) Colors based on \emph{lobe}}
#'   \item{color.class,color.network}{(vertex and edge) If applicable}
#'   \item{circle.layout}{(vertex) Integer vector for ordering the vertices for
#'     plots with circular layout}
#' @name CreateGraphs
#' @rdname make_brainGraph
#' @family Graph creation functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

make_brainGraph <- function(x, atlas, type=c('observed', 'random'),
                            level=c('subject', 'group', 'contrast'), set.attrs=TRUE,
                            modality=NULL, weighting=NULL, threshold=NULL,
                            name=NULL, Group=NULL, ...) {
  UseMethod('make_brainGraph')
}

#' Create a brainGraph object from an igraph graph
#'
#' @export
#' @method make_brainGraph igraph
#' @rdname make_brainGraph

make_brainGraph.igraph <- function(x, atlas, type=c('observed', 'random'),
                                   level=c('subject', 'group', 'contrast'),
                                   set.attrs=TRUE, modality=NULL, weighting=NULL,
                                   threshold=NULL, name=NULL, Group=NULL, ...) {
  lobe <- hemi <- index <- class <- network <- x.mni <- y.mni <- z.mni <- NULL

  x <- get_metadata(x)
  x$atlas <- if (missing(atlas)) guess_atlas(x) else atlas
  DT <- get(atlas)
  if (!is_named(x)) {
    V(x)$name <- DT$name
  } else {
    nonmatches <- !V(x)$name %in% DT[, name]
    if (any(nonmatches)) {
      stop(paste('Check the following vertex names: ',
                 paste(V(x)$name[nonmatches], collapse=' ')))
    }
  }

  vorder <- match(V(x)$name, DT$name)
  V(x)$lobe <- DT[vorder, as.character(lobe)]
  V(x)$lobe.hemi <- as.numeric(DT[vorder, interaction(lobe, hemi)])
  V(x)$hemi <- DT[vorder, as.character(hemi)]

  if (isTRUE(grepl('destr', x$atlas))) V(x)$class <- DT[vorder, as.numeric(class)]
  if (x$atlas == 'dosenbach160') V(x)$network <- DT[vorder, as.character(network)]

  level <- match.arg(level)
  type <- match.arg(type)
  x$level <- level
  x$type <- type
  attrs <- c('modality', 'weighting', 'threshold', 'name', 'Group')
  for (a in attrs) {
    if (!is.null(get(a))) x <- set_graph_attr(x, a, get(a))
  }
  if (level == 'group' && !is.null(Group)) x$name <- x$Group

  if (type == 'observed') {
    l.cir <- vector('integer')
    lobes <- DT[, levels(lobe)]
    V(x)$x <- V(x)$x.mni <- DT[vorder, x.mni]
    V(x)$y <- V(x)$y.mni <- DT[vorder, y.mni]
    V(x)$z <- V(x)$z.mni <- DT[vorder, z.mni]
    x <- set_graph_colors(x, 'color.lobe', DT[vorder, as.numeric(lobe)])
    if (x$atlas %in% c('destrieux', 'destrieux.scgm')) {
      x <- set_graph_colors(x, 'color.class', V(x)$class)
    } else if (x$atlas == 'dosenbach160') {
      x <- set_graph_colors(x, 'color.network', DT[vorder, as.numeric(network)])
      l.cir <- c(l.cir, which(V(x)$hemi == 'B'))
    }

    lobeorder <- list('Frontal', c('Insula', 'Central'), c('Limbic', 'Cingulate'),
                      'SCGM', 'Temporal', 'Parietal', 'Occipital', 'Cerebellum', 'Brainstem')
    if (!'Brainstem' %in% lobes) lobeorder <- lobeorder[-9]
    if (!'Cerebellum' %in% lobes) lobeorder <- lobeorder[-8]
    if (!'SCGM' %in% lobes) lobeorder <- lobeorder[-4]
    for (i in seq_along(lobeorder)) {
      l.cir <- c(l.cir, DT[lobe %in% lobeorder[[i]] & !hemi %in% c('B', 'R'),
                           .SD[order(-y.mni, x.mni), index]])
    }
    for (i in seq_along(lobeorder)) {
      l.cir <- c(l.cir, DT[lobe %in% rev(lobeorder)[[i]] & hemi == 'R',
                           .SD[order(y.mni, x.mni), index]])
    }
    V(x)$circle.layout <- l.cir
  }

  # Set a bunch of attributes
  if (isTRUE(set.attrs) && ecount(x) > 1) x <- set_brainGraph_attr(x, type, ...)

  class(x) <- c('brainGraph', class(x))
  return(x)
}

#' Create a brainGraph object from an adjacency matrix
#'
#' @param mode Character string defining how the matrix should be interpreted.
#'   Default: \code{'undirected'}
#' @param weighted Logical specifying whether to create a weighted network
#' @param diag Logical indicating whether to include the diagonal of the
#'   connectivity matrix. Default: \code{FALSE}
#' @export
#' @method make_brainGraph matrix
#'
#' @rdname make_brainGraph
#' @examples
#' \dontrun{
#' bg <- make_brainGraph(A, 'dkt', modality='dti', weighting='fa',
#'   mode='undirected', diag=FALSE, weighted=TRUE)
#' }

make_brainGraph.matrix <- function(x, atlas, type=c('observed', 'random'),
                                   level=c('subject', 'group', 'contrast'),
                                   set.attrs=TRUE, modality=NULL, weighting=NULL,
                                   threshold=NULL, name=NULL, Group=NULL,
                                   mode='undirected', weighted=NULL, diag=FALSE, ...) {
  g <- graph_from_adjacency_matrix(x, mode, weighted, diag)
  type <- match.arg(type)
  level <- match.arg(level)
  atlas <- if (missing(atlas)) guess_atlas(x) else atlas
  g <- make_brainGraph(g, atlas, type, level, set.attrs, modality, weighting,
                       threshold, name, Group, A=x, ...)
  return(g)
}

#' Create a graph with mediation-specific attributes
#'
#' This function only creates a graph for \emph{vertex}-level analyses.
#'
#' @param x A \code{bg_mediate} object
#' @param ... Other arguments passed to \code{\link{make_brainGraph}}
#' @inheritParams CreateGraphs
#' @export
#'
#' @return A \code{brainGraph_mediate} graph object with attributes:
#'   \item{Graph}{\emph{mediator}, \emph{treat}, \emph{outcome}, \emph{nobs}}
#'   \item{Vertex}{\emph{b?.acme, p?.acme}, \emph{b?.ade, p?.ade},
#'     \emph{b?.prop, p?.prop}, \emph{b.tot, p.tot}}
#' @family Graph creation functions

make_brainGraph.bg_mediate <- function(x, atlas=x$atlas, type='observed',
                                       level='contrast', set.attrs=FALSE,
                                       modality=NULL, weighting=NULL,
                                       threshold=NULL, ...) {
  stopifnot(inherits(x, 'bg_mediate'), x$level == 'vertex')
  med.sum <- summary(x)$DT
  g.med <- make_empty_brainGraph(atlas, ...)
  for (a in c('mediator', 'treat', 'outcome', 'nobs')) {
    g.med <- set_graph_attr(g.med, a, x[[a]])
  }
  for (a in c('b0.acme', 'b0.ade', 'b.tot', 'b0.prop')) {
    g.med <- set_vertex_attr(g.med, a, value=med.sum[[a]])
  }
  for (a in c('p0.acme', 'p0.ade', 'p.tot', 'p0.prop')) {
    g.med <- set_vertex_attr(g.med, a, value=1 - med.sum[[a]])
  }
  if (isTRUE(x$INT)) {
    for (a in c('b1.acme', 'b1.ade', 'b1.prop', 'b.avg.acme', 'b.avg.ade', 'b.avg.prop')) {
      g.med <- set_vertex_attr(g.med, a, value=med.sum[[a]])
    }
    for (a in c('p1.acme', 'p1.ade', 'p1.prop', 'p.avg.acme', 'p.avg.ade', 'p.avg.prop')) {
      g.med <- set_vertex_attr(g.med, a, value=1 - med.sum[[a]])
    }
  }
  class(g.med) <- c('brainGraph_mediate', class(g.med))
  return(g.med)
}

################################################################################
# OTHER METHODS
################################################################################

#' Determine whether the input is a brainGraph object
#'
#' @export
#' @rdname make_brainGraph
is.brainGraph <- function(x) inherits(x, 'brainGraph')

#' Print a summary of a brainGraph object
#'
#' @param object A \code{brainGraph} object
#' @param print.attrs Character string indicating whether or not to list the
#'   object's attributes (default: \code{all})
#' @export
#' @method summary brainGraph
#' @rdname make_brainGraph

summary.brainGraph <- function(object, print.attrs=c('all', 'graph', 'vertex', 'edge', 'none'), ...) {
  if (!is.brainGraph(object)) NextMethod(generic='summary', object=object)

  df <- print_bg_summary(object)
  if (object$level != 'subject') df <- df[-14, ]

  print.attrs <- match.arg(print.attrs)
  if (print.attrs == 'all') {
    attrtypes <- c('graph', 'vertex', 'edge')
  } else if (print.attrs == 'none') {
    attrtypes <- NULL
  } else {
    attrtypes <- print.attrs
  }
  attrs.l <- sapply(attrtypes, function(x) NULL)
  for (atype in attrtypes) {
    attrs <- switch(atype,
                    graph=graph_attr_names(object),
                    vertex=vertex_attr_names(object),
                    edge=edge_attr_names(object))
    if (length(attrs) > 0) attrs.l[[atype]] <- print_text_vector(attrs, 3)
  }
  out <- list(object=object, df=df, attrs=attrs.l, print.attrs=print.attrs)
  class(out) <- c('summary.brainGraph', class(out))
  out
}

#' @aliases summary.brainGraph
#' @method print summary.brainGraph

print.summary.brainGraph <- function(x, ...) {
  print_title_summary(paste0('Summary for *', x$object$type, '* ', x$object$level, '-level graph: ', x$object$name))
  print(x$df, right=FALSE, row.names=FALSE)
  cat('\n')

  if (x$print.attrs != 'none') {
    for (atype in names(x$attrs)) {
      if (is.null(x$attrs[[atype]])) {
        cat('No', tolower(atype), 'attributes!')
      } else {
        title <- paste(tools::toTitleCase(atype), 'attributes')
        width <- getOption('width') - nchar(title) - 1
        message(title, paste(rep('-', width / 2), collapse=''))
        print(x$attrs[[atype]], right=FALSE, row.names=FALSE)
      }
      cat('\n')
    }
  }
  invisible(x)
}

################################################################################
# OTHER CREATION FUNCTIONS
################################################################################

#' Create an empty graph with attributes for brainGraph
#'
#' \code{make_empty_brainGraph} creates an empty undirected \code{brainGraph}
#' object with vertex count equal to the atlas specified; i.e., it creates a
#' graph with 0 edges. Typically used to present results from an analysis in
#' which edges don't make sense (e.g., GLM comparing differences in a
#' vertex-level attribute).
#'
#' @export
#' @return \code{make_empty_brainGraph} -- An empty \code{brainGraph} graph
#'   object
#' @rdname make_brainGraph

make_empty_brainGraph <- function(atlas, type=c('observed', 'random'),
                                  level=c('subject', 'group', 'contrast'),
                                  modality=NULL, weighting=NULL, threshold=NULL,
                                  name=NULL, Group=NULL, ...) {
  n <- nrow(get(atlas))
  A <- matrix(0, nrow=n, ncol=n)
  type <- match.arg(type)
  level <- match.arg(level)
  g <- make_brainGraph(A, atlas, type, level, set.attrs=FALSE, modality,
                       weighting, threshold, name, Group, ...)
  return(g)
}

#' Create a graph of the union of multiple vertex neighborhoods
#'
#' This function accepts multiple vertices, creates graphs of their
#' neighborhoods (of order 1), and returns the union of those graphs.
#'
#' @param g An \code{igraph} graph object
#' @param vs Either a character or integer vector (vertex names or indices,
#' respectively) for the vertices of interest
#' @export
#'
#' @return An \code{igraph} graph object containing the union of all edges and
#'   vertices in the neighborhoods of the input vertices; only the vertex
#'   attribute \emph{name} will be present
#' @family Graph creation functions
#' @seealso \code{\link[igraph]{ego}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' subg <- make_ego_brainGraph(g1[[N]], c(24, 58))
#' subg <- make_ego_brainGraph(g1[[N]], c('lPCUN', 'rPCUN'))
#' }

make_ego_brainGraph <- function(g, vs) {

  subgs <- make_ego_graph(g, order=1, nodes=vs)
  if (is.character(vs)) vs <- which(V(g)$name %in% vs)

  for (i in seq_along(vs)) {
    subgs[[i]] <- delete_all_attr(subgs[[i]], keep.names=TRUE)
  }

  combine_graphs <- function(x, y) {
    n <- length(x)
    if (n < 2) {
      res <- x[[1]] %u% y
    } else {
      y <- x[[n]] %u% y
      x <- x[-n]
      res <- combine_graphs(x, y)
    }
    return(res)
  }

  inds <- unique(c(vs, unlist(lapply(vs, function(x) neighbors(g, x)))))
  subg.all <- combine_graphs(subgs, make_empty_graph(directed=F) + vertices(V(g)$name[inds]))
  return(subg.all)
}

#' Create the intersection of graphs based on a logical condition
#'
#' @param ... Graph objects or lists of graph objects
#' @param subgraph Character string specifying an equation (logical condition)
#'   for the vertices to subset
#' @export
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @return An \code{igraph} graph object
#' @examples
#' \dontrun{
#' res.mtpc <- mtpc(g, covars, ...)
#' g.mtpc <- make_glm_brainGraph(res.mtpc, atlas)
#' g.mtpc.int <- make_intersection_brainGraph(g.mtpc,
#'   subgraph='sig == 1')
#' }

make_intersection_brainGraph <- function(..., subgraph) {
  g <- inds <- NULL
  graphs <- args_as_list(...)
  stopifnot(all(sapply(graphs, inherits, 'brainGraph')))
  Nv <- vcount(graphs[[1]])

  subs <- lapply(graphs, subset_graph, subgraph)
  graphs.sub <- lapply(subs, with, g)
  inds.sub <- lapply(subs, with, inds)
  graphs.valid <- graphs.sub[which(sapply(graphs.sub, function(x) !is.null(x)))]

  if (length(graphs.valid) == 0) {
    return(make_empty_brainGraph(graphs[[1]]$atlas))
  } else if (length(graphs.valid) == 1) {
    return(graphs.valid[[1]])
  } else {
    g.int <- do.call(intersection, c(graphs.valid, keep.all.vertices=FALSE))
    memb <- which(V(graphs[[1]])$name %in% V(g.int)$name)
    g.int <- delete_all_attr(g.int)
    V(g.int)$name <- V(graphs[[1]])$name[memb]
    g.int <- graphs[[1]] %s% g.int
    g.int <- graphs[[1]] - vertices(setdiff(seq_len(Nv), memb))
    class(g.int) <- class(graphs[[1]])
    return(g.int)
  }
}
