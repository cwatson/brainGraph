################################################################################
# MAIN CREATION FUNCTIONS
################################################################################

#' Create a brainGraph object
#'
#' Create a \code{brainGraph} graph object, which is an \code{igraph} graph
#' object with additional attributes (at all levels). The values are dependent
#' on the specified brain atlas.
#'
#' For the \code{modality} argument, you can choose anything you like, but the
#' \code{summary.brainGraph} knows about \code{dti}, \code{fmri},
#' \code{thickness}, \code{area}, and \code{volume}.
#'
#' For the \code{weighting} argument, you can choose anything you like, but
#' \code{summary.brainGraph} knows about \code{fa}, \code{sld} (streamline
#' density, tractography), \code{pearson}, \code{spearman}, \code{kendall}, and
#' \code{partial} (partial correlation coefficient).
#'
#' @param g An \emph{igraph} graph object.
#' @param atlas Character string specifying the brain atlas
#' @param rand A character string indicating whether this function is being run
#'   for a random graph. Default: \code{FALSE}
#' @param modality Character vector indicating imaging modality (e.g. 'dti').
#'   Default: \code{NULL}
#' @param weighting Character string indicating how the edges are weighted
#'   (e.g., 'fa', 'pearson', etc.). Default: \code{NULL}
#' @param threshold Numeric indicating the level at which the matrices were
#'   thresholded (if at all). Default: \code{NULL}
#' @param subject Character vector indicating subject ID. Default: \code{NULL}
#' @param group Character vector indicating group membership. Default:
#'   \code{NULL}
#' @export
#'
#' @return A \code{brainGraph} graph object with additional attributes:
#'   \item{version}{(graph-level) The current version of \code{brainGraph}}
#'   \item{atlas}{(graph-level)}
#'   \item{lobe}{(vertex-leve) Character vector of lobe names}
#'   \item{hemi}{(vertex-leve) Character vector of hemispheres (\code{'L'},
#'     \code{'R'}, or \code{'B'})}
#'   \item{lobe.hemi}{Integer vector indicating the lobe and hemisphere}
#'   \item{class}{(vertex-leve) Character vector of class names (if applicable)}
#'   \item{network}{(vertex-leve) Character vector of network names (if
#'     applicable)}
#'   \item{modality}{(graph-level)}
#'   \item{weighting}{(graph-level)}
#'   \item{threshold}{(graph-level)}
#'   \item{name}{(graph-level) The subject ID (if specified by \code{subject})}
#'   \item{Group}{(graph-level) only if \code{group} is specified}
#'   \item{x, y, z, x.mni, y.mni, z.mni}{Spatial coordinates}
#'   \item{color.lobe}{(vertex- and edge-level) Colors based on \emph{lobe}}
#'   \item{color.class,color.network}{(vertex- and edge-level) If applicable}
#'   \item{circle.layout}{Integer vector for ordering the vertices for circle
#'     plots}
#' @family Graph creation functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

make_brainGraph <- function(obj, atlas, rand=FALSE, modality=NULL,
                            weighting=NULL, threshold=NULL, subject=NULL,
                            group=NULL, ...) {
  UseMethod('make_brainGraph', object=obj)
}

#' Create a brainGraph object from an igraph graph
#'
#' @param obj An \code{igraph} graph
#' @export
#' @method make_brainGraph igraph

make_brainGraph.igraph <- function(obj, atlas, rand=FALSE, modality=NULL,
                                   weighting=NULL, threshold=NULL, subject=NULL,
                                   group=NULL) {
  lobe <- hemi <- name <- index <- N <- class <- network <- x <- y <- z <-
    x.mni <- y.mni <- z.mni <- NULL

  obj$version <- packageVersion('brainGraph')
  obj$atlas <- atlas
  DT <- get(atlas)
  if (!is_named(obj)) {
    V(obj)$name <- DT$name
  } else {
    nonmatches <- !V(obj)$name %in% DT[, name]
    if (any(nonmatches)) {
      stop(paste('Check the following vertex names: ',
                 paste(V(obj)$name[nonmatches], collapse=' ')))
    }
  }

  vorder <- match(V(obj)$name, DT$name)
  lobe.nums <- DT[vorder, as.numeric(lobe)]
  V(obj)$lobe <- DT[vorder, as.character(lobe)]
  V(obj)$lobe.hemi <- as.numeric(DT[vorder, interaction(lobe, hemi)])
  V(obj)$hemi <- DT[vorder, as.character(hemi)]

  if (isTRUE(grepl('destr', obj$atlas))) V(obj)$class <- DT[vorder, as.numeric(class)]
  if (obj$atlas == 'dosenbach160') V(obj)$network <- DT[vorder, as.character(network)]

  if (!isTRUE(rand)) {
    # First add some "bookkeeping" attributes
    if (!is.null(modality)) obj$modality <- modality
    if (!is.null(weighting)) obj$weighting <- weighting
    if (!is.null(threshold)) obj$threshold <- threshold
    if (!is.null(subject)) obj$name <- subject
    if (!is.null(group)) obj$Group <- group

    l.cir <- vector('integer')
    lobes <- DT[, levels(lobe)]
    V(obj)$x <- V(obj)$x.mni <- DT[vorder, x.mni]
    V(obj)$y <- V(obj)$y.mni <- DT[vorder, y.mni]
    V(obj)$z <- V(obj)$z.mni <- DT[vorder, z.mni]
    V(obj)$color.lobe <- group.cols[lobe.nums]
    obj <- set_edge_color(obj, 'color.lobe', lobe.nums)
    if (obj$atlas %in% c('destrieux', 'destrieux.scgm')) {
      V(obj)$color.class <- group.cols[V(obj)$class]
      obj <- set_edge_color(obj, 'color.class', V(obj)$class)
    } else if (obj$atlas == 'dosenbach160') {
      V(obj)$color.network <- group.cols[DT[vorder, as.numeric(network)]]
      obj <- set_edge_color(obj, 'color.network', DT[vorder, as.numeric(network)])
      l.cir <- c(l.cir, which(V(obj)$hemi == 'B'))
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

    V(obj)$circle.layout <- l.cir
  }

  class(obj) <- c('brainGraph', class(obj))
  return(obj)
}

#' Create a brainGraph object from an adjacency matrix
#'
#' \code{make_brainGraph.matrix} creates a \code{brainGraph} object from an
#' adjacency matrix through \code{\link[igraph]{graph_from_adjacency_matrix}}.
#'
#' @param obj A numeric matrix
#' @param ... Arguments passed to
#'   \code{\link[igraph]{graph_from_adjacency_matrix}}
#' @export
#' @method make_brainGraph matrix
#'
#' @examples
#' \dontrun{
#' bg <- make_brainGraph(A, 'dkt', modality='dti', weighting='fa',
#'   mode='undirected', diag=FALSE, weighted=TRUE)
#' }

make_brainGraph.matrix <- function(obj, atlas, rand=FALSE, modality=NULL,
                                   weighting=NULL, threshold=NULL, subject=NULL,
                                   group=NULL, ...) {
  obj <- graph_from_adjacency_matrix(obj, ...)
  obj <- make_brainGraph(obj, atlas, rand, modality, weighting, threshold, subject, group)
  return(obj)
}

################################################################################
# OTHER METHODS
################################################################################

#' Determine whether the input is a brainGraph object
#'
#' @param obj An object to test
#' @keywords internal
#' @export
#' @rdname make_brainGraph
is.brainGraph <- function(obj) inherits(obj, 'brainGraph')

#' Print a summary of a brainGraph object
#'
#' @param object A \code{brainGraph} object
#' @param print.attrs Character string indicating whether or not to list the
#'   object's attributes (default: \code{all})
#' @param ... Unused
#' @export
#' @method summary brainGraph
#' @rdname make_brainGraph

summary.brainGraph <- function(object, print.attrs=c('all', 'none'), ...) {
  if (!is.brainGraph(object)) {
    NextMethod(generic='summary', object=object)
    return(invisible(object))
  }
  ver <- weighting <- name <- Group <- modality <- clustmethod <- 'N/A'

  if ('version' %in% graph_attr_names(object)) ver <- as.character(object$version)
  atlasfull <-
    switch(object$atlas,
           aal116='AAL-116', aal2.120=,aal2.94='AAL2', aal90='AAL-90',
           brainsuite='Brainsuite', craddock200='Craddock-200',
           destrieux='Destrieux', destrieux.scgm='Destrieux + SCGM',
           dk='Desikan-Killiany', dk.scgm='Desikan-Killiany + SCGM',
           dkt='Desikan-Killiany-Tourville', dkt.scgm='Desikan-Killiany-Tourville + SCGM',
           dosenbach160='Dosenbach-160', hoa112='Harvard-Oxford cortical and subcortical',
           lpba40='LONI probabilistic brain atlas', object$atlas)
  if ('modality' %in% graph_attr_names(object)) {
    modality <-
      switch(object$modality, dti='DTI', fmri='fMRI', thickness='Cortical thickness',
             area='Cortical surface area', volume='Cortical/subcortical volume', object$modality)
  }
  if (is_weighted(object)) {
    if ('weighting' %in% graph_attr_names(object)) {
      weighting <-
        switch(object$weighting, fa='FA (fractional anisotropy)',
               sld='Streamline density', pearson='Pearson correlation',
               spearman='Spearman\'s rank correlation',
               kendall='Kendall\'s rank correlation', partial='Partial correlation',
               object$weighting)
    }
  } else {
    weighting <- 'Unweighted'
  }
  if ('clust.method' %in% graph_attr_names(object)) {
    clustmethod <-
      switch(object$clust.method, edge_betweenness='Edge betweenness',
             fast_greedy='Greedy optimization (hierarchical agglomeration)',
             infomap='Infomap', label_prop='Label propagation',
             leading_eigen='Leading eigenvector',
             louvain='Louvain (multi-level modularity optimization)', optimal='Optimal',
             spinglass='Potts spin glass model', walktrap='Walktrap algorithm',
             object$clust.method)
  }
  dens.pct <- sprintf('%1.2f%s', 100 * graph.density(object), '%')
  if ('name' %in% graph_attr_names(object)) name <- object$name
  if ('Group' %in% graph_attr_names(object)) Group <- object$Group

  df <- data.frame(A=c('brainGraph version: ', 'Brain atlas used: ', 'Imaging modality: ',
                       'Edge weighting: ', 'Clustering method: ', 'Graph density: ',
                       'Subject ID: ', 'Group: '),
                   B=c(ver, atlasfull, modality, weighting, clustmethod, dens.pct, name, Group))
  dimnames(df)[[2]] <- rep('', 2)

  attrtypes <- c('graph', 'vertex', 'edge')
  attrs.l <- sapply(attrtypes, function(x) NULL)
  for (type in attrtypes) {
    attrs <- switch(type,
                    graph=graph_attr_names(object),
                    vertex=vertex_attr_names(object),
                    edge=edge_attr_names(object))
    len <- length(attrs)
    if (len > 0) {
      div <- seq_len(len)
      factors <- div[len %% div == 0L]
      if (len > 6 && length(factors) == 2L) {
        attrs <- c(attrs, '')
        len <- length(attrs)
        div <- seq_len(len)
        factors <- div[len %% div == 0L]
      }
      mod <- max(factors[which(factors <= 6)])
      attrs.df <- as.data.frame(split(attrs, ceiling(seq_along(attrs) / (len %/% mod))))
      dimnames(attrs.df)[[2]] <- rep('', ncol(attrs.df))
      attrs.l[[type]] <- attrs.df
    }
  }
  print.attrs <- match.arg(print.attrs)
  out <- list(df=df, attrs=attrs.l, print.attrs=print.attrs)
  class(out) <- c('summary.brainGraph', class(out))
  out
}

#' @aliases summary.brainGraph
#' @method print summary.brainGraph

print.summary.brainGraph <- function(x, ...) {
  print(x$df, right=FALSE, row.names=FALSE)
  cat('\n')

  if (x$print.attrs == 'all') {
    for (type in names(x$attrs)) {
      if (is.null(x$attrs[[type]])) {
        cat('No', tolower(type), 'attributes!')
      } else {
        title <- paste(type, 'attributes')
        width <- getOption('width') - nchar(title) - 1
        message(title, paste(rep('-', width / 2), collapse=''))
        print(x$attrs[[type]], right=FALSE, row.names=FALSE)
      }
      cat('\n')
    }
  }
  invisible(x)
}

#' Create an empty graph with attributes for brainGraph
#'
#' This function creates an empty undirected graph with vertex count equal to
#' the atlas specified, and includes some graph-, vertex-, and
#' edge-level attributes that are important for \code{brainGraph} functions.
#' Basically a wrapper for \code{\link[igraph]{make_empty_graph}}.
#'
#' @param atlas Character string of the atlas to create a graph from
#' @param ... Other arguments passed to \code{\link{make_brainGraph}}
#' @export
#'
#' @return \code{make_empty_brainGraph} -- An empty \code{brainGraph} graph
#'   object
#' @family Graph creation functions
#' @rdname make_brainGraph
#' @seealso \code{\link[igraph]{make_empty_graph}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

make_empty_brainGraph <- function(atlas, rand=FALSE, modality=NULL,
                                  weighting=NULL, threshold=NULL, subject=NULL,
                                  group=NULL,...) {
  n <- nrow(get(atlas))
  A <- matrix(0, nrow=n, ncol=n)
  g <- make_brainGraph(A, atlas, rand, modality, weighting, threshold, subject, group, ...)
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
#' @seealso \code{\link[igraph]{make_ego_graph}}
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

#' Create a graph with GLM-specific attributes
#'
#' \code{make_glm_brainGraph} will create graphs with attributes specific to the
#' results of \code{\link{brainGraph_GLM}} or \code{\link{mtpc}}. The function
#' returns a list, with one element for each specified contrast.
#'
#' This function only creates a graph for \emph{vertex}-level analyses.
#'
#' @param res.glm List as output by \code{\link{brainGraph_GLM}} or by
#' \code{\link{mtpc}}.
#' @param atlas Character string specifying the brain atlas to use
#' @param ... Other arguments passed to \code{\link{make_brainGraph}}
#' @export
#'
#' @return A list of \code{igraph} graph objects (length equal to the number of
#'   contrasts) with additional attributes:
#'   \item{Graph}{\emph{name} (contrast name), \emph{outcome} (the outcome
#'     variable), \emph{alpha} (the significance level); for MTPC:
#'     \emph{tau.mtpc}, \emph{S.mtpc}, \emph{S.crit}, \emph{A.crit}}
#'   \item{Vertex}{\emph{size2} (t-statistic), \emph{size} (the t-stat
#'     transformed for visualization purposes), \emph{p} (equal to \eqn{1-p}),
#'     \emph{p.fdr} (equal to \eqn{1-p_{FDR}}, the FDR-adjusted p-value),
#'     \emph{gamma} (the contrast of parameter estimaties, \emph{se} (the
#'     standard error of \emph{gamma}); \emph{A.mtpc}, \emph{sig} (binary
#'     indicating whether \code{A.mtpc > A.crit}) (for MTPC)}
#' @family Graph creation functions
#' @seealso \code{\link{brainGraph_GLM}, \link{mtpc}}

make_glm_brainGraph <- function(res.glm, atlas, ...) {
  contrast <- p <- p.fdr <- p.perm <- se <- stat <- A.mtpc <- region <- A.crit <- S.crit <- S.mtpc <- tau.mtpc <- NULL
  check.class <- inherits(res.glm, c('bg_GLM', 'mtpc'), which=TRUE)
  stopifnot(any(check.class == 1), res.glm$level == 'vertex')

  g.diffs <- vector('list', length=length(res.glm$con.name))
  for (i in seq_along(g.diffs)) {
    g.diffs[[i]] <- make_empty_brainGraph(atlas, ...)
    g.diffs[[i]]$name <- res.glm$con.name[i]
    g.diffs[[i]]$con.type <- res.glm$con.type
    g.diffs[[i]]$outcome <- res.glm$outcome
    g.diffs[[i]]$alt <- res.glm$alt

    if (check.class[1] == 1) {  # bg_GLM
      g.diffs[[i]]$alpha <- res.glm$alpha
      V(g.diffs[[i]])$p <- 1 - res.glm$DT[contrast == i, p]
      V(g.diffs[[i]])$p.fdr <- 1 - res.glm$DT[contrast == i, p.fdr]
      V(g.diffs[[i]])$gamma <- res.glm$DT[contrast == i, gamma]
      V(g.diffs[[i]])$se <- res.glm$DT[contrast == i, se]
      V(g.diffs[[i]])$size2 <- res.glm$DT[contrast == i, stat]
      V(g.diffs[[i]])$size <- vec.transform(V(g.diffs[[i]])$size2, 0, 20)
      if (isTRUE(res.glm$permute)) V(g.diffs[[i]])$p.perm <- 1 - res.glm$DT[contrast == i, p.perm]
      class(g.diffs[[i]]) <- c('brainGraph_GLM', class(g.diffs[[i]]))
    } else {  # mtpc
      g.diffs[[i]]$tau.mtpc <- res.glm$stats[contrast == i, tau.mtpc]
      g.diffs[[i]]$S.mtpc <- res.glm$stats[contrast == i, S.mtpc]
      g.diffs[[i]]$S.crit <- res.glm$stats[contrast == i, S.crit]
      g.diffs[[i]]$A.crit <- res.glm$stats[contrast == i, A.crit]
      V(g.diffs[[i]])$A.mtpc <- res.glm$DT[contrast == i, unique(A.mtpc), by=region]$V1
      V(g.diffs[[i]])$sig <- 0
      V(g.diffs[[i]])[res.glm$DT[contrast == i & A.mtpc > A.crit, unique(region)]]$sig <- 1
      class(g.diffs[[i]]) <- c('brainGraph_mtpc', class(g.diffs[[i]]))
    }
  }
  return(g.diffs)
}

#' Create a graph with NBS-specific attributes
#'
#' @param res.nbs List that is output by \code{\link{NBS}}
#' @param atlas Character string specifying the brain atlas to use
#' @param ... Other arguments passed to \code{\link{make_brainGraph}}
#' @export
#'
#' @return A list of \code{igraph} graph objects (length equal to the number of
#'   contrasts) with additional attributes:
#'   \item{Graph}{\emph{name} (contrast name)}
#'   \item{Vertex}{\emph{comp} (integer vector indicating connected component
#'     membership), \emph{p.nbs} (P-value for each component)}
#'   \item{Edge}{\emph{stat} (the test statistic for each connection), \emph{p}
#'     (the P-value)}
#' @family Graph creation functions

make_nbs_brainGraph <- function(res.nbs, atlas, ...) {
  contrast <- p.perm <- csize <- NULL
  stopifnot(inherits(res.nbs, 'NBS'))
  g.nbs <- vector('list', length=length(res.nbs$con.name))
  for (i in seq_along(g.nbs)) {
    g.nbs[[i]] <- graph_from_adjacency_matrix(res.nbs$T.mat[, , i], diag=F, mode='undirected', weighted=TRUE)
    g.nbs[[i]]$name <- res.nbs$con.name[i]
    g.nbs[[i]]$con.type <- res.nbs$con.type
    g.nbs[[i]]$alt <- res.nbs$alt
    if (ecount(g.nbs[[i]]) > 0) {
      E(g.nbs[[i]])$stat <- E(g.nbs[[i]])$weight
      E(g.nbs[[i]])$p <- 1 - E(graph_from_adjacency_matrix(res.nbs$p.mat[, , i], diag=F, mode='undirected', weighted=TRUE))$weight
      if (any(E(g.nbs[[i]])$weight < 0)) g.nbs[[i]] <- delete_edge_attr(g.nbs[[i]], 'weight')
      clusts <- components(g.nbs[[i]])
      comps <- sort(unique(clusts$csize), decreasing=TRUE)
      x <- clusts$membership
      x.tab <- table(x)
      x.tab.st <- sort(x.tab, decreasing=TRUE)
      V(g.nbs[[i]])$comp <- match(x, order(x.tab, decreasing=TRUE))
      V(g.nbs[[i]])$p.nbs <- 0
      xdt <- copy(res.nbs$components$observed)
      for (j in seq_along(comps)) {
        inds <- which(xdt[contrast == i, csize[j]] == x.tab.st)
        V(g.nbs[[i]])[V(g.nbs[[i]])$comp %in% inds]$p.nbs <- 1 - xdt[contrast == i, p.perm[j]]
      }
      if (ecount(g.nbs[[i]]) > 1) {
        g.nbs[[i]] <- set_brainGraph_attr(g.nbs[[i]], atlas=atlas, ...)
      } else {
        g.nbs[[i]] <- make_brainGraph(g.nbs[[i]], atlas=atlas, ...)
      }
    }
    class(g.nbs[[i]]) <- c('brainGraph_NBS', class(g.nbs[[i]]))
  }
  return(g.nbs)
}

#' Create a graph with mediation-specific attributes
#'
#' This function only creates a graph for \emph{vertex}-level analyses.
#'
#' @param res.med List object output by \code{\link{brainGraph_mediate}}
#' @param atlas Character string specifying the brain atlas to use
#' @param ... Other arguments passed to \code{\link{make_brainGraph}}
#' @export
#'
#' @return A \code{brainGraph_mediate} graph object with attributes:
#'   \item{Graph}{\emph{mediator}, \emph{treat}, \emph{outcome}, \emph{nobs}}
#'   \item{Vertex}{\emph{b?.acme, p?.acme}, \emph{b?.ade, p?.ade},
#'     \emph{b?.prop, p?.prop}, \emph{b.tot, p.tot}}
#' @family Graph creation functions

make_mediate_brainGraph <- function(res.med, atlas, ...) {
  stopifnot(inherits(res.med, 'bg_mediate'), res.med$level == 'vertex')
  med.sum <- summary(res.med)$DT
  g.med <- make_empty_brainGraph(atlas, ...)
  g.med$mediator <- res.med$mediator
  g.med$treat <- res.med$treat
  g.med$outcome <- res.med$outcome
  g.med$nobs <- res.med$nobs
  V(g.med)$b0.acme <- med.sum$b0.acme
  V(g.med)$p0.acme <- 1 - med.sum$p0.acme
  V(g.med)$b0.ade <- med.sum$b0.ade
  V(g.med)$p0.ade <- 1 - med.sum$p0.ade
  V(g.med)$b.tot <- med.sum$b.tot
  V(g.med)$p.tot <- 1 - med.sum$p.tot
  V(g.med)$b0.prop <- med.sum$b0.prop
  V(g.med)$p0.prop <- 1 - med.sum$p0.prop
  if (isTRUE(res.med$INT)) {
    V(g.med)$b1.acme <- med.sum$b1.acme
    V(g.med)$p1.acme <- 1 - med.sum$p1.acme
    V(g.med)$b1.ade <- med.sum$b1.ade
    V(g.med)$p1.ade <- 1 - med.sum$p1.ade
    V(g.med)$b1.prop <- med.sum$b1.prop
    V(g.med)$p1.prop <- 1 - med.sum$p1.prop
    V(g.med)$b.avg.acme <- med.sum$b.avg.acme
    V(g.med)$b.avg.acme <- med.sum$b.avg.acme
    V(g.med)$p.avg.acme <- med.sum$p.avg.acme
    V(g.med)$b.avg.ade <- med.sum$b.avg.ade
    V(g.med)$p.avg.ade <- med.sum$p.avg.ade
    V(g.med)$b.avg.prop <- med.sum$b.avg.prop
    V(g.med)$p.avg.prop <- med.sum$p.avg.prop
  }
  class(g.med) <- c('brainGraph_mediate', class(g.med))
  return(g.med)
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
