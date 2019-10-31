################################################################################
# MAIN CREATION FUNCTIONS
################################################################################

#' Create a list of brainGraph graphs
#'
#' \code{make_brainGraphList} creates a \code{brainGraphList} object, a list
#' containing a set of graphs for all subjects (or group-average graphs) in a
#' study at a specific threshold (or density), in addition to some graph-level
#' attributes common to those graphs.
#'
#' In addition to creating the initial \code{igraph} graphs from the
#' connectivity matrices, \code{\link{set_brainGraph_attr}} will be called on
#' each graph if \code{set.attrs=TRUE}; other arguments will be passed to that
#' function. You may display a progress bar by setting \code{.progress=TRUE}.
#'
#' This object can be considered comparable to a 4-D \emph{NIfTI} file,
#' particularly that returned by FSL's \emph{TBSS} \dQuote{prestats} step since
#' that file contains the FA volumes for all study subjects.
#'
#' @note If the input is a \code{corr_mats} object, and the extent of the 3-D
#' array is greater than 1, then only the first will be converted to a graph.
#'
#' @param x 3-D numeric array of all subjects' connectivity matrices (for a
#'   single threshold) or a \code{corr_mats} object
#' @param gnames Character vector of graph names (e.g., study IDs if
#'   \code{level='subject'}). Default: \code{NULL}
#' @param ... Other arguments passed to \code{\link{set_brainGraph_attr}}
#' @inheritParams Creating_Graphs
#' @export
#'
#' @return \code{make_brainGraphList} returns an object of class
#'   \code{brainGraphList} with elements:
#'   \item{threshold}{The specified threshold/density}
#'   \item{version}{The versions of \code{R}, \code{igraph}, and
#'     \code{brainGraph} used when creating the graphs}
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
                                threshold=NULL, gnames=NULL, ...) {
  UseMethod('make_brainGraphList')
}

#' @param grpNames Character (or factor) vector of group names. If \code{level ==
#'   'group'}, then you do not need to include this argument (the group names
#'   will be the same as \code{gnames}). Default: \code{NULL})
#' @param .progress Logical indicating whether to print a progress bar. Default:
#'   \code{getOption('bg.progress')}
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach getDoParRegistered
#' @rdname brainGraphList

make_brainGraphList.array <- function(x, atlas, type=c('observed', 'random'),
                                      level=c('subject', 'group', 'contrast'),
                                      set.attrs=TRUE, modality=NULL,
                                      weighting=NULL, threshold=NULL,
                                      gnames=NULL, grpNames=NULL, subnet=NULL,
                                      mode='undirected', weighted=NULL, diag=FALSE,
                                      .progress=getOption('bg.progress'), ...) {
  i <- NULL
  level <- match.arg(level)
  kNumGraphs <- dim(x)[3L]
  if (is.null(gnames)) gnames <- seq_len(kNumGraphs)
  stopifnot(length(gnames) == kNumGraphs)
  gnames <- as.character(gnames)
  if (is.null(grpNames)) {
    grpNames <- if (level == 'group') gnames else rep(1, kNumGraphs)
  }
  stopifnot(length(grpNames) == kNumGraphs)
  grpNames <- as.character(grpNames)

  type <- match.arg(type)
  attrs <- c('atlas', 'type', 'level', 'modality', 'weighting', 'threshold')
  out <- setNames(vector('list', length(attrs)), attrs)
  for (a in attrs) {
    if (!is.null(get(a))) out[[a]] <- get(a)
  }
  out <- get_metadata(out)
  if (!is.null(threshold)) threshold <- rep_len(threshold, kNumGraphs)

  # Show a progress bar so you aren't left in the dark
  #---------------------------------------------------------
  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  loopfun <- make_brainGraph
  if (isTRUE(.progress)) {
    ncpu <- getOption('bg.ncpus')
    print(paste('Start time:', format(as.POSIXct(Sys.time()), '%Y-%m-%d %H:%M:%OS0')))
    progbar <- txtProgressBar(min=0, max=kNumGraphs, style=3)
    loopfun <- function(...) {
      curVal <- get('counter', envir=env) + ncpu
      assign('counter', curVal, envir=env)
      setTxtProgressBar(get('progbar', envir=env), curVal)
      flush.console()
      make_brainGraph(...)
    }
  }

  env <- environment()
  counter <- 0
  g <- foreach(i=seq_len(kNumGraphs)) %dopar% {
    res <- loopfun(x[, , i], atlas, type, level, set.attrs, modality, weighting, threshold[i],
                   name=gnames[i], Group=grpNames[i], subnet=subnet, mode=mode,
                   diag=diag, weighted=weighted, use.parallel=FALSE, ...)
  }
  if (isTRUE(.progress)) close(progbar)

  out$graphs <- g
  names(out$graphs) <- gnames
  class(out) <- c('brainGraphList', 'list')
  return(out)
}

#' @export
#' @rdname brainGraphList

make_brainGraphList.corr_mats <- function(x, atlas=x$atlas, type='observed',
                                          level='group',
                                          set.attrs=TRUE, modality=NULL,
                                          weighting=NULL, threshold=x$densities,
                                          gnames=names(x$r.thresh), grpNames=gnames,
                                          mode='undirected', weighted=NULL,
                                          diag=FALSE, .progress=getOption('bg.progress'), ...) {
  dims <- vapply(x$r.thresh, dim, numeric(3L))[3L, ]
  if (any(dims > 1L)) x <- x[1]
  A <- abind::abind(x$r.thresh)
  if (isTRUE(weighted)) A <- x$R * A
  out <- make_brainGraphList(A, atlas=atlas, type=type, level=level,
                             set.attrs=set.attrs, modality=modality, weighting=weighting,
                             threshold=threshold, mode=mode, weighted=weighted, diag=diag,
                             gnames=gnames, grpNames=grpNames, .progress=.progress, ...)
  return(out)
}

#' Create a graph list with GLM-specific attributes
#'
#' These methods create a \code{brainGraphList} with attributes specific to the
#' results of \code{\link{brainGraph_GLM}}, \code{\link{mtpc}}, or
#' \code{\link{NBS}}. The \code{graphs} element of the returned object will
#' contain one graph for each contrast.
#'
#' @note Only valid for \emph{vertex}-level and \emph{NBS} analyses.
#'
#' @param x A \code{bg_GLM}, \code{mtpc}, or \code{NBS} object
#' @param atlas Character string specifying the brain atlas to use
#' @inheritParams brainGraphList
#' @export
#'
#' @return A \code{brainGraphList} object, with a graph object for each contrast
#'   with additional attributes:
#'   \item{Graph}{\emph{name} (contrast name), \emph{outcome} (the outcome
#'     variable), \emph{alpha} (the significance level); for MTPC:
#'     \emph{tau.mtpc}, \emph{S.mtpc}, \emph{S.crit}, \emph{A.crit}}
#'   \item{Vertex}{\emph{size2} (t-statistic); \emph{size} (the t-stat
#'     transformed for visualization purposes); \emph{p} (equal to \eqn{1-p});
#'     \emph{p.fdr} (equal to \eqn{1-p_{FDR}}, the FDR-adjusted p-value);
#'     \emph{effect.size} (the contrast of parameter estimates for t-contrasts;
#'     the extra sum of squares for F-contrasts); \emph{se} (the
#'     standard error of \emph{gamma}); \emph{A.mtpc}, \emph{sig} (binary
#'     indicating whether \code{A.mtpc > A.crit}) (for MTPC)}
#' @name Creating_Graphs_GLM
#' @rdname glm_brainGraphList
#' @family Graph creation functions
#' @seealso \code{\link{brainGraph_GLM}, \link{mtpc}, \link{NBS}}

make_brainGraphList.bg_GLM <- function(x, atlas=x$atlas, type='observed',
                                       level='contrast', set.attrs=FALSE,
                                       modality=NULL, weighting=NULL,
                                       threshold=NULL, gnames=x$con.name, ...) {
  contrast <- p <- p.fdr <- p.perm <- se <- stat <- region <- ESS <- NULL
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
    if (x$con.type == 't') {
      V(g.diffs[[i]])$effect.size <- x$DT[contrast == i, gamma]
    } else {
      V(g.diffs[[i]])$effect.size <- x$DT[contrast == i, ESS]
    }
    V(g.diffs[[i]])$se <- x$DT[contrast == i, se]
    V(g.diffs[[i]])$size2 <- x$DT[contrast == i, stat]
    V(g.diffs[[i]])$size <- vec.transform(V(g.diffs[[i]])$size2, 0, 20)
    if (isTRUE(x$permute)) V(g.diffs[[i]])$p.perm <- 1 - x$DT[contrast == i, p.perm]
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
#' @rdname glm_brainGraphList

make_brainGraphList.mtpc <- function(x, atlas, type='observed', level='contrast',
                                     set.attrs=FALSE, modality=NULL, weighting=NULL,
                                     threshold=NULL, gnames=x$con.name, ...) {
  contrast <- region <- A.mtpc <- A.crit <- S.crit <- S.mtpc <- tau.mtpc <- NULL
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
#' @export
#' @return \code{make_brainGraphList.NBS} returns graphs with additional
#'   attributes:
#'   \item{Vertex}{\emph{comp} (integer vector indicating connected component
#'     membership), \emph{p.nbs} (P-value for each component)}
#'   \item{Edge}{\emph{stat} (the test statistic for each connection), \emph{p}
#'     (the P-value)}
#' @rdname glm_brainGraphList

make_brainGraphList.NBS <- function(x, atlas, type='observed', level='contrast',
                                    set.attrs=TRUE, modality=NULL, weighting=NULL,
                                    threshold=NULL, gnames=x$con.name,
                                    mode='undirected', weighted=TRUE, diag=FALSE,
                                    ...) {
  contrast <- p.perm <- csize <- NULL

  bgList <- make_brainGraphList(x$T.mat, atlas, type, level, set.attrs=set.attrs,
      modality=modality, weighting=weighting, threshold=threshold, mode=mode,
      weighted=weighted, diag=diag, gnames=gnames, ...)
  pList <- make_brainGraphList(x$p.mat, atlas, level=level, set.attrs=FALSE, weighted=TRUE, ...)

  for (i in seq_along(bgList[])) {
    class(bgList$graphs[[i]]) <- c('brainGraph_NBS', class(bgList[i]))
    bgList$graphs[[i]]$con.type <- x$con.type
    bgList$graphs[[i]]$alt <- x$alt

    if (ecount(bgList[i]) == 0) next
    E(bgList$graphs[[i]])$stat <- E(bgList[i])$weight
    E(bgList$graphs[[i]])$p <- 1 - E(pList[i])$weight
    if (any(E(bgList[i])$weight < 0)) bgList$graphs[[i]] <- delete_edge_attr(bgList[i], 'weight')

    clusts <- components(bgList[i])
    comps <- sort(unique(clusts$csize), decreasing=TRUE)
    memb <- clusts$membership
    x.tab <- table(memb)
    x.tab.st <- sort(x.tab, decreasing=TRUE)
    V(bgList$graphs[[i]])$comp <- match(memb, order(x.tab, decreasing=TRUE))
    V(bgList$graphs[[i]])$p.nbs <- 0
    xdt <- copy(x$components$observed)
    for (j in seq_along(comps)) {
      inds <- which(xdt[contrast == i, csize[j]] == x.tab.st)
      V(bgList$graphs[[i]])[V(bgList[i])$comp %in% inds]$p.nbs <- 1 - xdt[contrast == i, p.perm[j]]
    }
  }
  return(bgList)
}

################################################################################
# OTHER METHODS
################################################################################

#' Index subjects and/or groups in a brainGraphList
#'
#' The \code{[} method will let you subset/slice the graphs for individual
#' subjects and/or \emph{groups}.
#'
#' @section Subsetting/extracting:
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
#' # Return object with graphs from groups 1 and 3
#' my.bgl[g=c(1, 3), drop=FALSE]
#'
#' # Subset the first 10 subjects of group 2
#' my.bgl[1:10, 2]
#'}

`[.brainGraphList` <- function(x, i, g=NULL, drop=TRUE) {
  stop_msg <- 'Logical indexing vector must be of the same length as the number of groups.'
  warn_msg <- 'Indexing vector cannot have length greater than the number of groups.'
  # Subset groups; doesn't make sense for contrast lists
  if (!is.null(g) && x$level != 'contrast') {
    if (is.factor(g)) g <- as.character(g)
    group.vec <- groups(x)
    grpNames <- unique(group.vec)
    kNumGroups <- length(grpNames)
    if (is.logical(g) && length(g) != kNumGroups) stop(stop_msg)
    if (!is.logical(g) && length(g) > kNumGroups) {
      warning(warn_msg)
      g <- g[seq_len(kNumGroups)]
    }
    if (is.numeric(g) || is.logical(g)) g <- grpNames[g]
    if (is.character(g)) g <- which(group.vec %in% g)
    x$graphs <- x$graphs[g]
  }

  # Subset subjects
  kNumSubs <- length(x$graphs)
  if (missing(i)) {
    if (isTRUE(drop)) return(x$graphs)
    i <- seq_len(kNumSubs)
  }
  if (is.logical(i) && length(i) != kNumSubs) stop(sub('group', 'subject', stop_msg))
  if (!is.logical(i) && length(i) > kNumSubs) {
    warning(sub('group', 'subject', warn_msg))
    i <- i[seq_len(kNumSubs)]
  }

  if (is.character(i)) i <- which(names(x$graphs) %in% i)
  x$graphs <- if (length(i) == 1 && isTRUE(drop)) x$graphs[[i]] else x$graphs[i]
  out <- if (isTRUE(drop)) x$graphs else x
  return(out)
}

#' @export
print.brainGraphList <- function(x, ...) {
  kNumGraphs <- length(x$graphs)
  gnames <- names(x$graphs)
  print_title_summary(paste0('A "brainGraphList" object of *', x$type,
                             '* graphs containing ', kNumGraphs, ' ', x$level, 's.'))

  df <- print_bg_summary(x)
  df <- df[-c(6, 10, 11, 13, 14), ]
  print(df, right=FALSE, row.names=FALSE)
  cat('\n')

  # Print subject/group names
  msg <- switch(x$level, subject='IDs:', 'names:')
  message(paste(simpleCap(x$level), msg))
  if (kNumGraphs < 10) {
    print(gnames)
  } else {
    attrs.df <- print_text_vector(gnames, 6)
    print(attrs.df, row.names=FALSE)
    cat('\n')
  }
  if (x$level == 'subject' && 'Group' %in% graph_attr_names(x$graphs[[1]])) {
    message('Group membership:')
    print(table(groups(x)))
  }
  invisible(x)
}

#' Coerce list of graphs to a brainGraphList object
#'
#' \code{as_brainGraphList} coerces a list of graphs to a \code{brainGraphList}
#' object. It is assumed that certain metadata attributes -- threshold, package
#' version, atlas, imaging modality, edge weighting, and whether they are
#' random graphs -- are identical for all graphs in the list.
#'
#' To convert an object with 3 \dQuote{levels} (i.e., subject-level lists from
#' an older \code{brainGraph} version), see the code in the Examples below.
#'
#' @param g.list List of graph objects
#' @export
#'
#' @rdname brainGraphList
#' @examples
#' \dontrun{
#' ## Convert old version single-subject graph lists
#' ## g[[1]] is group 1, g[[1]][[1]] is threshold 1, g[[1]][[1]][[1]] is subj. 1
#' kNumThresholds <- length(g[[1]])
#' g.l <- vector('list', kNumThresholds)
#' for (i in seq_len(kNumThresholds)) {
#'   g.l[[i]] <- as_brainGraphList(do.call(Map, c(c, g))[[i]])
#' }
#' }

as_brainGraphList <- function(g.list, type=c('observed', 'random'),
                              level=c('subject', 'group', 'contrast')) {
  if (!inherits(g.list, 'list')) g.list <- list(g.list)

  type <- match.arg(type)
  level <- match.arg(level)
  if (type == 'observed' || level == 'group') {
    stopifnot(all(vapply(g.list, inherits, logical(1), 'brainGraph')))
    g1 <- g.list[[1]]
    ids <- vapply(g.list, graph_attr, character(1), 'name')
  } else {
    stopifnot(all(vapply(g.list, function(x) is_igraph(x[[1]]), logical(1))))
    g1 <- g.list[[1]][[1]]
    ids <- vapply(g.list, function(x) graph_attr(x[[1]], 'name'), character(1))
  }

  attrnames <- graph_attr_names(g1)
  attrs <- c('atlas', 'modality', 'weighting', 'threshold', 'version', 'sys', 'date')
  out <- setNames(vector('list', length=9), c('type', 'level', attrs))
  out$type <- type
  out$level <- level
  for (x in attrs) {
    if (x %in% attrnames) out[[x]] <- graph_attr(g1, x)
  }
  if (length(out$version < 3)) out$version <- list(r='', bg=out$version, ig='')
  for (x in c('sys', 'date')) if (is.null(out[[x]])) out[[x]] <- ''
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
#' @inheritParams plot.brainGraph
#' @export
#' @importFrom methods hasArg
#' @importFrom grDevices dev.off pdf
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot.brainGraphList <- function(x, plane, hemi, filename.base, diffs=FALSE, ...) {
  pdf(file=sprintf('%s.pdf', filename.base), onefile=TRUE)

  kNumGraphs <- length(x$graphs)
  plot(x[1], plane=plane, hemi=hemi, main=x[1]$name, ...)
  for (i in seq.int(2, kNumGraphs)) {
    plot(x[i], plane=plane, hemi=hemi, main=x[i]$name, ...)
    if (isTRUE(diffs)) {
      g.diff <- graph.difference(x$graphs[[i]], x$graphs[[i-1]])
      class(g.diff) <- c('brainGraph', class(g.diff))
      ecols <- rep('deeppink', ecount(g.diff))
      if (hasArg('edge.color')) ecols <- list(...)$edge.color
      plot(g.diff, add=TRUE, vertex.label=NA, vertex.shape='none',
           edge.width=5, edge.color=ecols, mni=FALSE, subtitle=NULL)
    }
  }
  dev.off()
}
