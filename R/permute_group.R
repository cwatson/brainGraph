#' Permutation test for group difference of graph measures
#'
#' \code{brainGraph_permute} draws permutations from linear model residuals to
#' determine the significance of between-group differences of a global or
#' vertex-wise graph measure. It is intended for structural covariance networks
#' (in which there is only one graph per group), but can be extended to other
#' types of data.
#'
#' If you would like to calculate differences in the area-under-the-curve (AUC)
#' across densities, then specify \code{auc=TRUE}.
#'
#' There are three possible \dQuote{levels}:
#' \enumerate{
#'   \item \emph{graph} Calculate modularity (Louvain algorithm), clustering
#'   coefficient, characteristic path length, degree assortativity, and global
#'   efficiency.
#'   \item \emph{vertex} Choose one of: centrality metrics (betweenness,
#'   closeness, communicability, eigenvector, leverage, pagerank, subgraph);
#'   k-core; degree; eccentricity; nodal or local efficiency; k-nearest neighbor
#'   degree; shortest path length; transitivity; or vulnerability.
#'   \item \emph{other} Supply your own function. This is useful if you want to
#'   calculate something that I haven't hard-coded. It must take as its own
#'   arguments: \code{g} (a list of lists of \code{igraph} graph objects); and
#'   \code{densities} (numeric vector).
#' }
#'
#' @param densities Numeric vector of graph densities
#' @param resids An object of class \code{brainGraph_resids} (the output from
#'   \code{\link{get.resid}})
#' @param N Integer; the number of permutations (default: 5e3)
#' @param perms Numeric matrix of permutations, if you would like to provide
#'   your own (default: \code{NULL})
#' @param auc Logical indicating whether or not to calculate differences in the
#'   area-under-the-curve of metrics (default: \code{FALSE})
#' @param level A character string for the attribute \dQuote{level} to calculate
#'   differences (default: \code{graph})
#' @param measure A character string specifying the vertex-level metric to
#'   calculate, only used if \code{level='vertex'} (default: \code{btwn.cent}).
#'   For the \code{summary} method, this is to focus on a single
#'   \emph{graph-level} measure (since multiple are calculated at once).
#' @param .function A custom function you can pass if \code{level='other'}
#' @export
#' @importFrom permute shuffleSet
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#'
#' @return An object of class \code{brainGraph_permute} with input arguments in
#'   addition to:
#'   \item{DT}{A data table with permutation statistics}
#'   \item{obs.diff}{A data table of the observed group differences}
#'   \item{Group}{Group names}
#'
#' @family Group analysis functions
#' @family Structural covariance network functions
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' myResids <- get.resid(lhrh, covars)
#' myPerms <- shuffleSet(n=nrow(myResids$resids.all), nset=1e3)
#' out <- brainGraph_permute(densities, m, perms=myPerms)
#' out <- brainGraph_permute(densities, m, perms=myPerms, level='vertex')
#' out <- brainGraph_permute(densities, m, perms=myPerms,
#'   level='other', .function=myFun)
#' }

brainGraph_permute <- function(densities, resids, N=5e3, perms=NULL, auc=FALSE,
                               level=c('graph', 'vertex', 'other'),
                               measure=c('btwn.cent', 'coreness', 'degree', 'eccentricity',
                                         'clo.cent', 'communicability', 'ev.cent', 'lev.cent',
                                         'pagerank', 'subg.cent', 'E.local', 'E.nodal',
                                         'knn', 'Lp', 'transitivity', 'vulnerability'),
                               .function=NULL) {
  stopifnot(inherits(resids, 'brainGraph_resids'))
  gID <- getOption('bg.group')
  measure <- match.arg(measure)
  level <- match.arg(level)
  if (level == 'other') {
    if (!is.function(.function)) stop('".function" must be a function!')
    measure <- 'other'
  } else if (level == 'graph') {
    measure <- NULL
  }

  if (is.null(perms)) perms <- shuffleSet(n=nobs(resids), nset=N)
  dims <- dim(perms)
  N <- dims[1L]
  perms <- rbind(perms, seq_len(dims[2L]))  # last row is observed metrics
  grps <- as.numeric(resids$resids.all[, get(gID)])

  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  # Should return a "data.table" (graph-level) or numeric matrix (vertex-level)
  # When "auc=FALSE", should have a "densities" column
  res.perm <- switch(level,
               vertex=permute_vertex_foreach(perms, densities, resids, grps, auc, measure),
               other=permute_other_foreach(perms, densities, resids, grps, .function),
               graph=permute_graph_foreach(perms, densities, resids, grps, auc))

  if (length(densities) == 1L) res.perm <- cbind(densities=densities, res.perm)
  if (level == 'vertex') {
    res.perm <- as.data.table(res.perm)
    start <- if (isTRUE(auc)) 1L else 2L
    setnames(res.perm, seq.int(start, dim(res.perm)[2L]), region.names(resids))
  }

  obs.ind <- if (isTRUE(auc)) N + 1L else (N + 1L) * seq_along(densities)
  obs.diff <- res.perm[obs.ind]
  res.perm <- res.perm[-obs.ind]
  out <- list(atlas=resids$atlas, auc=auc, N=N, level=level, measure=measure, densities=densities,
              resids=resids, DT=res.perm, obs.diff=obs.diff, Group=resids$Group)
  class(out) <- c('brainGraph_permute', class(out))
  return(out)
}

#==============================================================================
# Helper functions
#==============================================================================
make_graphs_perm <- function(densities, resids, inds, grps) {
  corrs <- lapply(unique(grps), function(x)
                  corr.matrix(resids[which(grps[inds] == x)],
                              densities=densities, rand=TRUE))
  sapply(corrs, lapply, function(x)
         apply(x$r.thresh, 3L, graph_from_adjacency_matrix, mode='undirected', diag=FALSE))
}

# Graph level
graph_attr_perm <- function(g, densities) {
  mod <- sapply(g, sapply, function(x) modularity(cluster_louvain(x)))
  Cp <- sapply(g, sapply, function(x) transitivity(x, type='localaverage'))
  Lp <- sapply(g, sapply, mean_distance)
  assort <- sapply(g, sapply, assortativity_degree)
  E.global <- sapply(g, sapply, efficiency, 'global')
  list(mod=mod, Cp=Cp, Lp=Lp, assort=assort, E.global=E.global)
}

permute_graph_foreach <- function(perms, densities, resids, grps, auc) {
  i <- NULL
  N <- dim(perms)[1L]
  if (isTRUE(auc)) {
    res.perm <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
      g <- make_graphs_perm(densities, resids, perms[i, ], grps)
      meas.list <- graph_attr_perm(g, densities)
      t(vapply(meas.list, function(y) auc_diff_perm(densities, y), numeric(1L)))
    }
    res.perm <- as.data.table(res.perm)
  } else {
    # Returns a list of (Ndensities * Nperms) X 2 matrices for all metrics
    res.perm <- foreach(i=seq_len(N), .combine=function(a, b) Map(rbind, a, b)) %dopar% {
      g <- make_graphs_perm(densities, resids, perms[i, ], grps)
      meas.list <- graph_attr_perm(g, densities)
    }
    res.perm <- data.table(densities=rep.int(densities, N),
                           sapply(res.perm, function(x) as.numeric(-diff(t(x)))))
    setkey(res.perm, densities)
  }
  return(res.perm)
}

# Vertex-level
vertex_attr_funs <- function(measure) {  # Move the "switch" outside the "foreach" loop
  switch(measure,
         coreness=coreness,
         degree=degree,
         eccentricity=eccentricity,
         clo.cent=function(x) centr_clo(x)$res,
         communicability=centr_betw_comm,
         ev.cent=function(x) centr_eigen(x)$vector,
         lev.cent=function(x) centr_lev(x)$res,
         pagerank=function(x) page_rank(x, weights=NA)$vector,
         subg.cent=subgraph_centrality,
         E.local=function(x) efficiency(x, type='local', weights=NA, use.parallel=FALSE),
         E.nodal=function(x) efficiency(x, type='nodal', weights=NA),
         knn=function(x) knn(x)$knn,
         Lp=function(x) mean_distance_wt(x, level='vertex', weights=NA),
         transitivity=function(x) transitivity(x, type='local', isolates='zero'),
         vulnerability=function(x) vulnerability(x, use.parallel=FALSE),
         function(x) centr_betw(x)$res)
}

permute_vertex_foreach <- function(perms, densities, resids, grps, auc, measure) {
  i <- NULL
  if (isTRUE(auc)) {
    diffFun <- function(densities, meas.list) {
      sapply(seq_len(dim(meas.list[[1L]])[2L]), function(x)
             auc_diff(densities, cbind(meas.list[[1L]][, x], meas.list[[2L]][, x])))
    }
  } else {
    diffFun <- function(densities, meas.list) cbind(densities, meas.list[[1L]] - meas.list[[2L]])
  }
  fun <- vertex_attr_funs(measure)
  res.perm <- foreach(i=seq_len(dim(perms)[1L]), .combine='rbind') %dopar% {
    g <- make_graphs_perm(densities, resids, perms[i, ], grps)
    meas.list <- lapply(g, function(x) t(sapply(x, fun)))
    diffFun(densities, meas.list)
  }
}

# Other-level
permute_other_foreach <- function(perms, densities, resids, grps, .function) {
  i <- NULL
  res.perm <- foreach(i=seq_len(dim(perms)[1L]), .combine='rbind') %dopar% {
    g <- make_graphs_perm(densities, resids, perms[i, ], grps)
    .function(g, densities)
  }
}

#==============================================================================
# Methods
#==============================================================================

#' Print a summary from a permutation analysis
#'
#' @param object,x A \code{brainGraph_permute} object (output by
#'   \code{\link{brainGraph_permute}}).
#' @param p.sig Character string specifying which p-value to use for displaying
#'   significant results (default: \code{p})
#' @param ... Unused
#' @inheritParams GLM
#' @export
#' @rdname brainGraph_permute

summary.brainGraph_permute <- function(object, measure=object$measure,
                                       alternative=c('two.sided', 'less', 'greater'),
                                       alpha=0.05, p.sig=c('p', 'p.fdr'), ...) {
  perm.diff <- p <- p.fdr <- region <- obs.diff <- ..measure <- NULL
  gID <- getOption('bg.group')
  if (object$level == 'other') {  # Hack to figure out which level it is when level="other"
    object$level <- if (dim(object$DT)[2L] > 6L) 'vertex' else 'graph'
  }

  if (object$level == 'graph' && is.null(measure)) measure <- 'mod'
  group_str <- paste0(measure, '.', object$Group)
  permDT <- copy(object$DT)
  g <- with(object, make_graphs_perm(densities, resids, seq_len(nobs(resids)),
                                     resids$resids.all[, as.numeric(get(gID))]))

  densities <- object$densities
  # VERTEX-LEVEL
  #-------------------------------------
  if (object$level == 'vertex') {
    obsDT <- copy(object$obs.diff)
    fun <- vertex_attr_funs(measure)
    obs <- lapply(g, function(x) t(sapply(x, fun)))

    if (isTRUE(object$auc)) {
      obs <- lapply(obs, apply, 2L, function(y) -auc_diff(densities, y))
      permDT[, densities := 1]
      obsDT[, densities := 1]
      densities <- 1
    }
    sum.dt <- data.table(densities=densities,
                         region=rep(V(g[[1L]][[1L]])$name, each=length(densities)),
                         g1=c(obs[[1L]]), g2=c(obs[[2L]]),
                         key=c('densities', 'region'))
    setnames(sum.dt, c('g1', 'g2'), group_str)
    obs.m <- melt(obsDT, id.vars='densities', variable.name='region', value.name='obs.diff')
    sum.dt <- merge(sum.dt, obs.m)
    permDT <- melt(permDT, id.vars='densities', variable.name='region', value.name=measure)
    setkeyv(permDT, key(sum.dt))
    permdiff <- permDT[, list(perm.diff=mean(get(measure))), by=key(sum.dt)]
    sum.dt <- merge(sum.dt, permdiff, by=key(sum.dt))

  # GRAPH-LEVEL
  #-------------------------------------
  } else if (object$level == 'graph') {
    stopifnot(hasName(permDT, measure))
    permDT[, region := 'graph']
    obs <- graph_attr_perm(g, densities)[[measure]]

    if (isTRUE(object$auc)) {
      obs <- t(apply(obs, 2L, function(y) -auc_diff(densities, y)))
      permDT[, densities := 1]
      densities <- 1
    }

    sum.dt <- data.table(densities=densities, region='graph', obs)
    setnames(sum.dt, c('V1', 'V2'), group_str)
    sum.dt[, obs.diff := object$obs.diff[[measure]]]
    sum.dt[, perm.diff := permDT[, mean(get(measure)), by=densities]$V1]
  }
  result.dt <- merge(permDT[, c('densities', 'region', ..measure)],
                     sum.dt[, c('densities', 'region', 'obs.diff')],
                     by=c('densities', 'region'))

  alt <- match.arg(alternative)
  compfun <- switch(alt,
                    two.sided=function(perm, obs) sum(abs(perm) >= abs(unique(obs))),
                    less=function(perm, obs) sum(perm <= unique(obs)),
                    greater=function(perm, obs) sum(perm >= unique(obs)))
  result.dt[, p := (compfun(get(measure), obs.diff) + 1L) / (.N + 1L), by=key(result.dt)]
  CI <- switch(alt, two.sided=c(alpha / 2, 1 - (alpha / 2)), less=c(alpha, 1), greater=c(1 / object$N, 1 - alpha))
  result.dt[, c('ci.low', 'ci.high') := as.list(sort(get(measure))[ceiling(.N * CI)]), by=key(result.dt)]
  result.dt <- result.dt[, .SD[1L], by=key(result.dt)]
  sum.dt <- merge(sum.dt, result.dt[, !c(measure, 'obs.diff'), with=FALSE], by=key(result.dt))
  setcolorder(sum.dt, c('densities', 'region', group_str, 'obs.diff', 'perm.diff', 'ci.low', 'ci.high', 'p'))
  if (isFALSE(object$auc)) {
    groupby <- if (object$level == 'graph') NULL else 'densities'
    sum.dt[, p.fdr := p.adjust(p, 'fdr'), by=groupby]
  }

  meas.full <- print_full_metric(measure)
  p.sig <- match.arg(p.sig)
  perm.sum <- c(object, list(DT.sum=sum.dt, meas.full=meas.full, alt=alt, alpha=alpha, p.sig=p.sig))
  perm.sum$measure <- measure
  class(perm.sum) <- c('summary.brainGraph_permute', class(perm.sum))
  perm.sum
}

#' @aliases summary.brainGraph_permute
#' @method print summary.brainGraph_permute
#' @export

print.summary.brainGraph_permute <- function(x, ...) {
  print_title_summary('Permutation analysis')
  cat('# of permutations:', prettyNum(x$N, ','), '\n')
  cat('Level: ', x$level, '\n')
  cat('Graph metric: ', x$meas.full, '\n')
  if (isTRUE(x$auc)) cat('Area-under-the-curve (AUC) calculated across', length(x$densities), 'densities:\n', x$densities, '\n')

  symb <- switch(x$alt, two.sided='!=', greater='>', less='<')
  alt <- sprintf('%s - %s %s 0', x$Group[1L], x$Group[2L], symb)
  cat('Alternative hypothesis: ', '\t', alt, '\n')
  cat('Alpha: ', x$alpha, '\n\n')
  if (with(x, dim(DT.sum[get(p.sig) < alpha])[1L]) == 0L) {
    cat('No significant results!\n')
  } else {
    clp <- 100 * (1 - x$alpha)
    setnames(x$DT.sum, c('ci.low', 'ci.high'), paste0(clp, '% CI ', c('low', 'high')))
    with(x, print(DT.sum[get(p.sig) < alpha]))
  }
  invisible(x)
}

#' Plot results from permutation testing
#'
#' @param ptitle Character string specifying a title for the plot (default:
#'   \code{NULL})
#' @export
#' @return The \code{plot} method returns a \emph{list} of \code{ggplot} objects
#' @rdname brainGraph_permute

plot.brainGraph_permute <- function(x, measure=x$measure,
                                    alternative=c('two.sided', 'less', 'greater'),
                                    alpha=0.05, p.sig=c('p', 'p.fdr'), ptitle=NULL, ...) {
  densities <- Group <- sig <- trend <- yloc <- obs <- mylty <- ci.low <- ci.high <-
    variable <- value <- reg.num <- region <- perm.diff <- obs.diff <- plot2 <- NULL
  p.sig <- match.arg(p.sig)
  perm.sum <- summary(x, measure=measure, alternative=alternative, alpha=alpha)
  measure <- perm.sum$measure
  sum.dt <- perm.sum$DT.sum
  if (is.null(ptitle)) ptitle <- perm.sum$meas.full
  ylabel2 <- sprintf('Observed and permutation difference (%s - %s)', x$Group[1L], x$Group[2L])

  # GRAPH LEVEL
  #-------------------------------------
  if (x$level == 'graph') {
    plot.dt <- melt(sum.dt, id.vars=setdiff(names(sum.dt), paste0(measure, '.', x$Group)),
                    variable.name='Group', value.name='obs')
    plot.dt[, Group := factor(Group, labels=x$Group)]
    idvars <- c('densities', 'region', 'p', 'Group', 'obs')
    if (isFALSE(x$auc)) idvars <- c(idvars, 'p.fdr')
    plot.dt <- melt(plot.dt, id.vars=idvars)
    plot.dt[, c('sig', 'trend') := '']
    plot.dt[get(p.sig) < alpha, sig := '*']
    plot.dt[get(p.sig) >= alpha & get(p.sig) < 2 * alpha, trend := '*']
    plot.dt[, yloc := extendrange(obs, f=0.07)[1L]]
    plot.dt[, mylty := 0]
    plot.dt[variable == 'obs.diff', mylty := 1]
    plot.dt[variable %in% c('ci.low', 'ci.high'), mylty := 2]
    plot.dt[variable == 'perm.diff', mylty := 3]
    # Line plot of observed group values
    p <- ggplot(plot.dt, aes(x=densities))
    plot1 <- p +
      geom_line(aes(y=obs, col=Group)) +
      geom_text(aes(y=yloc, label=sig), col='red', size=8) +
      geom_text(aes(y=yloc, label=trend), col='blue', size=8) +
      theme(legend.position='bottom')

    # Line plot of observed difference and CI's
    plot2 <- p +
      geom_line(data=plot.dt[variable =='obs.diff'], aes(y=value, lty=factor(mylty)), col='red') +
      geom_point(data=plot.dt[variable == 'obs.diff'], aes(y=value), col='red', size=3) +
      geom_line(data=plot.dt[variable == 'perm.diff'], aes(y=value, lty=factor(mylty))) +
      geom_line(data=plot.dt[variable == 'ci.low'], aes(y=value, lty=factor(mylty))) +
      geom_line(data=plot.dt[variable == 'ci.high'], aes(y=value, lty=factor(mylty))) +
      scale_linetype_manual(name='',
                            labels=c('Observed diff.', '95% CI', 'Mean permutation diff.'),
                            values=1:3) +
      theme(legend.position='bottom', legend.spacing.x=unit(9, 'pt'),
            legend.key=element_rect(fill='white'), legend.background=element_rect(fill='gray79'))

    ylabel1 <- 'Observed values'
    plot1 <- plot1 +
      labs(title=ptitle, y=ylabel1, x='Density') +
      theme(plot.title=element_text(hjust=0.5, face='bold'))
    plot2 <- plot2 +
      labs(title=ptitle, y=ylabel2, x='Density') +
      theme(plot.title=element_text(hjust=0.5, face='bold'))

  # VERTEX LEVEL
  #-------------------------------------
  } else {
    plot.dt <- droplevels(sum.dt[get(p.sig) < alpha])
    if (dim(plot.dt)[1L] == 0L) stop('No significant results!')
    plot.dt[, reg.num := seq_len(.N), by=densities]
    plot1 <- ggplot(plot.dt, aes(x=region))
      geom_point(aes(y=perm.diff)) +
      geom_errorbar(aes(ymin=ci.low, ymax=ci.high)) +
      geom_segment(aes(x=reg.num - .25, xend=reg.num + .25, y=obs.diff, yend=obs.diff),
                   col='red', size=1.25) +
      facet_wrap(~ densities, scales='free') +
      labs(title=ptitle, y=ylabel2, x='Region') +
      theme(plot.title=element_text(hjust=0.5, face='bold'))
  }
  return(list(plot1, plot2))
}
