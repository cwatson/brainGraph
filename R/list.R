################################################################################
# MAIN CREATION FUNCTIONS
################################################################################

#' Create a list of brainGraph graphs
#'
#' \code{make_brainGraphList} creates a \code{brainGraphList} object, a list
#' containing a set of graphs for all subjects in a study at a specific
#' threshold (or density), in addition to some graph-level attributes common to
#' those graphs.
#'
#' In addition to creating the initial \code{igraph} graphs from the
#' connectivity matrices, \code{\link{set_brainGraph_attr}} will be called on
#' each graph if \code{set.attrs=TRUE}; other arguments will be passed to that
#' function. You may display a progress bar by setting \code{.progress=TRUE}.
#'
#' This object can be considered comparable to a 4-D \emph{NIfTI} file,
#' particularly that returned by FSL's \emph{TBSS} prestats step since that file
#' contains the FA volumes for all study subjects.
#'
#' @note If the input is a \code{corr_mats} object, and the extent of the 3-D
#' array is greater than 1, then only the first will be converted to a graph.
#'
#' @param x 3-D numeric array of all subjects' connectivity matrices (for a
#'   single threshold) or a \code{corr_mats} object
#' @param gnames Character vector of graph names (e.g., study IDs if
#'   \code{level='subject'}). Default: \code{NULL}
#' @param groups Character (or factor) vector of group names (default:
#'   \code{NULL})
#' @param .progress Logical indicating whether to print a progress bar. Default:
#'   \code{TRUE}
#' @param ... Other arguments passed to \code{\link{set_brainGraph_attr}}
#' @inheritParams CreateGraphs
#' @export
#'
#' @return \code{make_brainGraphList} returns an object of class
#'   \code{brainGraphList} with elements:
#'   \item{threshold}{The specified threshold/density}
#'   \item{version}{The version of \code{brainGraph} used when creating the
#'     graphs}
#'   \item{atlas}{The atlas common to all the graphs}
#'   \item{modality}{The imaging modality (if supplied)}
#'   \item{weighting}{A string indicating what edge weights represent (if
#'     applicable)}
#'   \item{graphs}{A \emph{named list} of \code{brainGraph} graphs; the names
#'     correspond to the individual graphs' Study IDs}
#'
#' @name brainGraphList
#' @rdname brainGraphList
#' @family Graph creation functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' # Create a list, one for each threshold
#' g <- vector('list', length(thresholds))
#' for (i in seq_along(thresholds)) {
#'   g[[i]] <- make_brainGraphList(A.norm.sub[[i]], thresholds[i], atlas,
#'       covars.dti$Study.ID, covars.dti$Group, modality='dti', weighting='fa')
#' }
#' }

make_brainGraphList <- function(x, atlas, type=c('observed', 'random'),
                                level=c('subject', 'group', 'contrast'),
                                set.attrs=TRUE, modality=NULL, weighting=NULL,
                                threshold=NULL, gnames=NULL, groups=NULL, ...) {
  UseMethod('make_brainGraphList')
}

#' @export
#' @method make_brainGraphList array
#' @rdname brainGraphList

make_brainGraphList.array <- function(x, atlas, type=c('observed', 'random'),
                                      level=c('subject', 'group', 'contrast'),
                                      set.attrs=TRUE, modality=NULL,
                                      weighting=NULL, threshold=NULL,
                                      gnames=NULL, groups=NULL,
                                      mode='undirected', weighted=NULL, diag=FALSE,
                                      .progress=TRUE, ...) {
  i <- NULL
  level <- match.arg(level)
  kNumGraphs <- dim(x)[3]
  if (!is.null(gnames)) {
    stopifnot(length(gnames) == kNumGraphs)
    if (is.factor(gnames)) gnames <- as.character(gnames)
  } else {
    gnames <- seq_len(kNumGraphs)
  }
  if (!is.null(groups)) {
    stopifnot(length(groups) == kNumGraphs)
    if (is.factor(groups)) groups <- as.character(groups)
  } else {
    if (level == 'group') {
      groups <- gnames
    } else {
      groups <- rep(1, kNumGraphs)
    }
  }

  type <- match.arg(type)
  attrs <- c('atlas', 'type', 'level', 'modality', 'weighting', 'threshold')
  out <- setNames(vector('list', length(attrs)), attrs)
  for (a in attrs) {
    if (!is.null(get(a))) out[[a]] <- get(a)
  }
  out <- get_metadata(out)

  # Show a progress bar so you aren't left in the dark
  #---------------------------------------------------------
  if (isTRUE(.progress)) {
    os <- .Platform$OS.type
    ncpu <- if (os == 'windows') as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) else detectCores()
    print(paste('Start time:', format(as.POSIXct(Sys.time()), '%Y-%m-%d %H:%M:%OS0')))
    progbar <- txtProgressBar(min=0, max=kNumGraphs, style=3)
    loopfun <- function(...) {
      curVal <- get('counter', envir=env) + ncpu
      assign('counter', curVal, envir=env)
      setTxtProgressBar(get('progbar', envir=env), curVal)
      flush.console()
      make_brainGraph(...)
    }
  } else {
    loopfun <- make_brainGraph
  }

  env <- environment()
  counter <- 0
  g <- foreach(i=seq_len(kNumGraphs)) %dopar% {
    res <- loopfun(x[, , i], atlas, type, level, set.attrs, modality, weighting, threshold,
                   name=gnames[i], Group=groups[i], mode='undirected',
                   diag=FALSE, weighted=TRUE, use.parallel=FALSE, ...)
  }
  if (isTRUE(.progress)) close(progbar)

  out$graphs <- g
  names(out$graphs) <- gnames
  class(out) <- c('brainGraphList', 'list')
  return(out)
}

#' @export
#' @method make_brainGraphList corr_mats
#' @rdname brainGraphList

make_brainGraphList.corr_mats <- function(x, atlas=x$atlas, type='observed',
                                          level='group',
                                          set.attrs=TRUE, modality=NULL,
                                          weighting=NULL, threshold=NULL,
                                          gnames=names(x$r.thresh), groups=gnames,
                                          mode='undirected', weighted=NULL,
                                          diag=FALSE, .progress=TRUE, ...) {
  A <- abind::abind(x$r.thresh)
  if (isTRUE(weighted)) A <- x$R * A
  out <- make_brainGraphList(A, atlas=atlas, type=type, level=level,
                             set.attrs=set.attrs, modality=modality, weighting=weighting,
                             threshold=threshold, mode=mode, weighted=weighted, diag=diag,
                             gnames=gnames, groups=groups, .progress=.progress, ...)
  return(out)
}

#' Create a graph list with GLM-specific attributes
#'
#' Create a \code{brainGraphList} with attributes specific to the results of
#' \code{\link{brainGraph_GLM}} or \code{\link{mtpc}}. The \code{graphs} element
#' of the returned object will contain one graph for each contrast.
#'
#' @note Only valid for \emph{vertex}-level analyses.
#'
#' @param x A \code{bg_GLM} or \code{mtpc} object
#' @param atlas Character string specifying the brain atlas to use
#' @inheritParams brainGraphList
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
#' @rdname make_brainGraphList.bg_GLM
#' @family Graph creation functions
#' @seealso \code{\link{brainGraph_GLM}, \link{mtpc}}

make_brainGraphList.bg_GLM <- function(x, atlas, type='observed', level='contrast',
                                       set.attrs=FALSE, modality=NULL, weighting=NULL,
                                       threshold=NULL, gnames=x$con.name, ...) {
  contrast <- p <- p.fdr <- p.perm <- se <- stat <- region <- NULL
  if (x$level == 'graph') stop('Not valid for graph-level results.')

  attrs <- c('atlas', 'type', 'level', 'modality', 'weighting', 'threshold')
  out <- setNames(vector('list', length(attrs)), attrs)
  for (a in attrs) {
    if (!is.null(get(a))) out[[a]] <- get(a)
  }
  out <- get_metadata(out)

  g.diffs <- vector('list', length=length(x$con.name))
  for (i in seq_along(g.diffs)) {
    g.diffs[[i]] <- make_empty_brainGraph(atlas, type='observed', level='contrast',
                                          name=x$con.name[i], ...)
    for (a in c('con.type', 'outcome', 'alt', 'alpha')) {
      g.diffs[[i]] <- set_graph_attr(g.diffs[[i]], a, x[[a]])
    }

    V(g.diffs[[i]])$p <- 1 - x$DT[contrast == i, p]
    V(g.diffs[[i]])$p.fdr <- 1 - x$DT[contrast == i, p.fdr]
    V(g.diffs[[i]])$gamma <- x$DT[contrast == i, gamma]
    V(g.diffs[[i]])$se <- x$DT[contrast == i, se]
    V(g.diffs[[i]])$size2 <- x$DT[contrast == i, stat]
    V(g.diffs[[i]])$size <- vec.transform(V(g.diffs[[i]])$size2, 0, 20)
    if (isTRUE(x$permute)) V(g.diffs[[i]])$p.perm <- 1 - res.glm$DT[contrast == i, p.perm]
    class(g.diffs[[i]]) <- c('brainGraph_GLM', class(g.diffs[[i]]))
  }
  out$graphs <- g.diffs
  names(out$graphs) <- gnames
  class(out) <- c('brainGraphList', 'list')
  return(out)
}

#' Create graph list for MTPC results
#'
#' @export
#' @rdname make_brainGraphList.bg_GLM

make_brainGraphList.mtpc <- function(x, atlas, type='observed', level='contrast',
                                     set.attrs=FALSE, modality=NULL, weighting=NULL,
                                     threshold=NULL, gnames=x$con.name, ...) {
  A.mtpc <- A.crit <- S.crit <- S.mtpc <- tau.mtpc <- NULL
  if (x$level == 'graph') stop('Not valid for graph-level results.')

  attrs <- c('atlas', 'type', 'level', 'modality', 'weighting', 'threshold')
  out <- setNames(vector('list', length(attrs)), attrs)
  for (a in attrs) {
    if (!is.null(get(a))) out[[a]] <- get(a)
  }
  out <- get_metadata(out)

  g.diffs <- vector('list', length=length(x$con.name))
  for (i in seq_along(g.diffs)) {
    g.diffs[[i]] <- make_empty_brainGraph(atlas, type, level, modality, weighting,
                                          threshold, name=gnames[i])
    for (a in c('con.type', 'outcome', 'alt', 'alpha')) {
      g.diffs[[i]] <- set_graph_attr(g.diffs[[i]], a, x[[a]])
    }
    for (a in c('tau.mtpc', 'S.mtpc', 'S.crit', 'A.crit')) {
      g.diffs[[i]] <- set_graph_attr(g.diffs[[i]], a, x$stats[contrast == i, get(a)])
    }
    V(g.diffs[[i]])$A.mtpc <- x$DT[contrast == i, unique(A.mtpc), by=region]$V1
    V(g.diffs[[i]])$sig <- 0
    sig.regions <- as.integer(x$DT[contrast == i & A.mtpc > A.crit, unique(region)])
    V(g.diffs[[i]])[sig.regions]$sig <- 1
    class(g.diffs[[i]]) <- c('brainGraph_mtpc', class(g.diffs[[i]]))
  }
  out$graphs <- g.diffs
  names(out$graphs) <- gnames
  class(out) <- c('brainGraphList', 'list')
  return(out)
}

#' Create a graph list for NBS results
#'
#' @param x A \code{NBS} object
#' @inheritParams brainGraphList
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

make_brainGraphList.NBS <- function(x, atlas, type='observed', level='contrast',
                                    set.attrs=TRUE, modality=NULL, weighting=NULL,
                                    threshold=NULL, gnames=x$con.name, groups=NULL,
                                    mode='undirected', weighted=TRUE, diag=FALSE,
                                    ...) {
  contrast <- p.perm <- csize <- NULL

  out <- make_brainGraphList(x$T.mat, atlas=atlas, type=type, level=level, set.attrs=set.attrs,
                             modality=modality, weighting=weighting, threshold=threshold,
                             mode=mode, weighted=weighted, diag=diag, gnames=gnames, ...)

  for (i in seq_along(out[])) {
    class(out$graphs[[i]]) <- c('brainGraph_NBS', class(out$graphs[[i]]))
    out$graphs[[i]]$con.type <- x$con.type
    out$graphs[[i]]$alt <- x$alt

    if (ecount(out$graphs[[i]]) > 0) {
      E(out$graphs[[i]])$stat <- E(out$graphs[[i]])$weight
      E(out$graphs[[i]])$p <- 1 - E(graph_from_adjacency_matrix(x$p.mat[, , i], diag=diag, mode=mode, weighted=TRUE))$weight
      if (any(E(out$graphs[[i]])$weight < 0)) out$graphs[[i]] <- delete_edge_attr(out$graphs[[i]], 'weight')

      clusts <- components(out$graphs[[i]])
      comps <- sort(unique(clusts$csize), decreasing=TRUE)
      memb <- clusts$membership
      x.tab <- table(memb)
      x.tab.st <- sort(x.tab, decreasing=TRUE)
      V(out$graphs[[i]])$comp <- match(memb, order(x.tab, decreasing=TRUE))
      V(out$graphs[[i]])$p.nbs <- 0
      xdt <- copy(x$components$observed)
      for (j in seq_along(comps)) {
        inds <- which(xdt[contrast == i, csize[j]] == x.tab.st)
        V(out$graphs[[i]])[V(out$graphs[[i]])$comp %in% inds]$p.nbs <- 1 - xdt[contrast == i, p.perm[j]]
      }
    }
  }
  return(out)
}

################################################################################
# OTHER METHODS
################################################################################

#' Index subjects and/or groups in a brainGraphList
#'
#' The \code{[} method will let you subset/slice the graphs for individual
#' subjects and/or \emph{groups}.
#'
#' The first index is for subsetting the individual graphs. The second index is
#' for subsetting by group membership and requires that the graphs have a
#' \code{Group} graph attribute. When both are included, the first index cannot
#' have length or numeric value greater than the number of \emph{remaining}
#' subjects \emph{after} subsetting by group.
#'
#' If the indexing vector(s) is (are) \code{character}, the vector(s) must
#' contain one (or more) of the subject or group names. If \code{logical}, its
#' length must equal the number of subjects or groups.
#'
#' @param i Integer, character, or logical vector for subsetting by subject, or
#' by group (if \code{x$level='group'})
#' @param g Integer, character, or logical vector for subsetting by group (if
#'   \code{x$level='subject'})
#' @param drop If \code{TRUE} (the default), then return only the list of
#'   graphs; otherwise, subset the graphs and return the entire object
#' @export
#' @method [ brainGraphList
#'
#' @return \code{[} -- A \code{brainGraphList} object (if \code{drop=FALSE}) or
#'   a list of graphs
#' @name Extract.brainGraphList
#' @rdname brainGraphList
#' @examples
#' \dontrun{
#' # Subset the first 10 subjects, irrespective of group
#' my.bgl[1:10]
#'
#' # Return object for only 'Control' subjects
#' my.bgl[, 'Control']
#'
#' # Return object for groups 1 and 3
#' my.bgl[, c(1, 3)]
#'
#' # Subset the first 10 subjects of group 2
#' my.bgl[1:10, 2]
#'}

`[.brainGraphList` <- function(x, i, g=NULL, drop=TRUE) {
  # Subset groups; doesn't make sense for contrast lists
  if (!is.null(g) && x$level != 'contrast') {
    group.vec <- vapply(x$graphs, graph_attr, character(1), 'Group')
    groups <- unique(group.vec)
    kNumGroups <- length(groups)
    if (is.logical(g) && length(g) != kNumGroups) {
      stop('Logical indexing vector must be of the same length as the number of groups.')
    }
    if (!is.logical(g) && length(g) > kNumGroups) {
      warning('Indexing vector cannot have length greater than the number of groups.')
      g <- g[1:kNumGroups]
    }
    if (is.numeric(g) || is.logical(g)) g <- groups[g]
    if (is.character(g)) g <- which(group.vec %in% g)
    x$graphs <- x$graphs[g]
  }

  # Subset subjects
  kNumSubjects <- length(x$graphs)
  if (missing(i)) {
    if (isTRUE(drop)) {
      return(x$graphs)
    } else {
      i <- seq_len(kNumSubjects)
    }
  }
  if (is.logical(i) && length(i) != kNumSubjects) {
    stop('Logical indexing vector must be of the same length as the number of subjects.')
  }
  if (!is.logical(i) && length(i) > kNumSubjects) {
    warning('Indexing vector cannot have length greater than the number of subjects..')
    i <- i[1:kNumSubjects]
  }

  if (is.character(i)) {
    graphs <- names(x$graphs)
    i <- which(graphs %in% i)
  }
  if (isTRUE(drop)) {
    if (length(i) == 1) {
      x$graphs <- x$graphs[[i]]
    } else {
      x$graphs <- x$graphs[i]
    }
    return(x$graphs)
  } else {
    return(x)
  }
}

#' @method print brainGraphList
#' @keywords internal

print.brainGraphList <- function(x, ...) {
  modality <- weighting <- thresh <- 'N/A'

  kNumGraphs <- length(x$graphs)
  gnames <- names(x$graphs)
  message(rep('=', getOption('width') / 1.5))
  message(paste0('A "brainGraphList" object of *', x$type, '* graphs'))
  message(paste0('containing ', kNumGraphs, ' ', x$level, 's.'))
  message(rep('=', getOption('width') / 1.5))

  #TODO: put this stuff into "method_helpers" or whatever
  ver <- sapply(x$version, as.character)
  date_created <- sub('T', ' ', x$date)
  atlasfull <-
    switch(x$atlas,
           aal116='AAL-116', aal2.120=,aal2.94='AAL2', aal90='AAL-90',
           brainsuite='Brainsuite', craddock200='Craddock-200',
           destrieux='Destrieux', destrieux.scgm='Destrieux + SCGM',
           dk='Desikan-Killiany', dk.scgm='Desikan-Killiany + SCGM',
           dkt='Desikan-Killiany-Tourville', dkt.scgm='Desikan-Killiany-Tourville + SCGM',
           dosenbach160='Dosenbach-160', hoa112='Harvard-Oxford cortical and subcortical',
           lpba40='LONI probabilistic brain atlas', x$atlas)
  if (!is.null(x$modality)) {
    modality <-
      switch(x$modality, dti='DTI', fmri='fMRI', thickness='Cortical thickness',
             area='Cortical surface area', volume='Cortical/subcortical volume', x$modality)
  }
  if (!is.null(x$weighting)) {
    weighting <-
        switch(x$weighting, fa='FA (fractional anisotropy)',
               sld='Streamline density', pearson='Pearson correlation',
               spearman='Spearman\'s rank correlation',
               kendall='Kendall\'s rank correlation', partial='Partial correlation',
               x$weighting)
  }
  if (!is.null(x$threshold)) thresh <- prettyNum(x$threshold, ',')

  df <- data.frame(A=c('Softare versions:',
                       '       R release:', '      brainGraph:', '          igraph:',
                       'Date created:',
                       'Brain atlas used:', 'Imaging modality:',
                       'Edge weighting:', 'Threshold:'),
                   B=c('', ver, date_created, atlasfull, modality, weighting, thresh))
  dimnames(df)[[2]] <- rep('', 2)
  print(df, right=FALSE, row.names=FALSE)
  cat('\n')

  # Print subject/group names
  cat(paste(tools::toTitleCase(x$level), 'names:'), '\n')
  if (kNumGraphs < 10) {
    print(gnames)
  } else {
    splits <- split(gnames, ceiling(seq_along(gnames) / (kNumGraphs %/% 4)))
    lens <- lengths(splits)
    nsplits <- length(splits)
    splits[[nsplits]] <- c(splits[[nsplits]], rep('', (kNumGraphs %/% 4) - lens[nsplits]))
    attrs.df <- as.data.frame(splits)
    dimnames(attrs.df)[[2]] <- rep('', ncol(attrs.df))
    print(attrs.df)
  }
  invisible(x)
}

#' Coerce list of graphs to a brainGraphList object
#'
#' \code{as_brainGraphList} will coerce a list of graphs to a
#' \code{brainGraphList} object. It is assumed that certain metadata attributes
#' -- threshold, package version, brain atlas, imaging modality, edge weighting,
#' and whether these are random graphs -- are identical for all graphs in the
#' list. For any that are not, the vector of values will be stored; this may be
#' an issue for downstream analysis.
#'
#' @param g.list List of graph objects
#' @export
#'
#' @rdname brainGraphList

as_brainGraphList <- function(g.list, type=c('observed', 'random'),
                              level=c('subject', 'group', 'contrast')) {
  if (!inherits(g.list, 'list')) g.list <- list(g.list)

  type <- match.arg(type)
  level <- match.arg(level)
  attrs <- c('atlas', 'modality', 'weighting', 'threshold', 'version', 'sys', 'date')
  out <- setNames(vector('list', length=8), c('type', 'level', attrs))
  out$type <- type
  out$level <- level
  if (type == 'observed') {
    stopifnot(all(vapply(g.list, inherits, logical(1), 'brainGraph')))
    g1 <- g.list[[1]]
    ids <- vapply(g.list, graph_attr, character(1), 'name')
  } else {
    stopifnot(all(vapply(g.list, function(x) is_igraph(x[[1]]), logical(1))))
    g1 <- g.list[[1]][[1]]
    ids <- vapply(g.list, function(x) graph_attr(x[[1]], 'name'), character(1))
  }

  for (x in attrs) {
    if (x %in% graph_attr_names(g1)) out[[x]] <- graph_attr(g1, x)
  }
  out$graphs <- g.list
  names(out$graphs) <- ids
  class(out) <- c('brainGraphList', 'list')
  return(out)
}

#' Plot a brainGraphList and write to PDF
#'
#' The \code{plot} method will write a PDF file containing plots for all graphs
#' in the given object.
#'
#' You can choose to highlight edge differences between subsequent list
#' elements; in this case, new/different edges are colored pink. This is useful
#' mostly for a list of group-level graphs.
#'
#' @param x A \code{brainGraphList} object
#' @param filename.base Character string specifying the base of the filename
#' @param diffs Logical, indicating whether edge differences should be
#'   highlighted. Default: \code{FALSE}
#' @param ... Other parameters (passed to \code{\link{plot.brainGraph}})
#' @export
#' @method plot brainGraphList
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot.brainGraphList <- function(x, filename.base, diffs=FALSE, ...) {
  pdf(file=sprintf('%s.pdf', filename.base), onefile=TRUE)

  for (i in seq_along(x$graphs)) {
    plot(x$graphs[[i]], main=x$graphs[[i]]$name, ...)
    if (isTRUE(diffs) && i > 1) {
      g.diff <- graph.difference(x$graphs[[i]], x$graphs[[i-1]])
      class(g.diff) <- c('brainGraph', class(g.diff))
      if (hasArg('edge.color')) {
        ecols <- list(...)$edge.color
      } else {
        ecols <- rep('deeppink', ecount(g.diff))
      }
      ewidth <- 5
      plot(g.diff, add=TRUE, vertex.label=NA, vertex.shape='none',
           edge.width=ewidth, edge.color=ecols, mni=FALSE, subt=NULL)
    }
  }
  dev.off()
}
