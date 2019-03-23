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
#' There are three possible "levels":
#' \enumerate{
#'   \item \emph{graph} Calculate modularity (Louvain algorithm), clustering
#'   coefficient, characteristic path length, degree assortativity, global
#'   efficiency, lobe assortativity, and edge asymmetry.
#'   \item \emph{vertex} Choose one of: betweenness centrality, degree, nodal
#'   efficiency, k-nearest neighbor degree, transitivity, or vulnerability.
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
#' @param level A character string for the attribute "level" to calculate
#'   differences (default: \code{graph})
#' @param atlas Character string of the atlas name; required if
#'   \code{level='graph'}
#' @param measure A character string specifying the vertex-level metric to
#'   calculate, only used if \code{level='vertex'} (default: \code{btwn.cent}).
#'   For the \code{summary} method, this is to focus on a single
#'   \emph{graph-level} measure (since multiple are calculated at once).
#' @param .function A custom function you can pass if \code{level='other'}
#' @export
#'
#' @return An object of class \code{brainGraph_permute} with input arguments in
#'   addition to:
#'   \item{DT}{A data table with permutation statistics}
#'   \item{obs.diff}{A data table of the observed group differences}
#'   \item{groups}{Group names}
#'
#' @family Group analysis functions
#' @family Structural covariance network functions
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' myResids <- get.resid(lhrh, covars)
#' myPerms <- shuffleSet(n=nrow(myResids$resids.all), nset=1e3)
#' out <- brainGraph_permute(densities, m, perms=myPerms, atlas='dk')
#' out <- brainGraph_permute(densities, m, perms=myPerms, level='vertex')
#' out <- brainGraph_permute(densities, m, perms=myPerms,
#'   level='other', .function=myFun)
#' }

brainGraph_permute <- function(densities, resids, N=5e3, perms=NULL, auc=FALSE,
                               level=c('graph', 'vertex', 'other'),
                               measure=c('btwn.cent', 'degree', 'E.nodal', 'ev.cent',
                                         'knn', 'transitivity', 'vulnerability'),
                               atlas=resids$atlas, .function=NULL) {
  Group <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'))
  measure <- match.arg(measure)
  level <- match.arg(level)
  if (level == 'other') {
    if (is.null(.function) | !is.function(.function)) {
      stop(paste('Argument ".function" must be a function!'))
    }
    measure <- 'other'
  } else if (level == 'graph') {
    stopifnot(!is.null(atlas))
    measure <- NULL
  } else if (level == 'vertex') {
    if (isTRUE(auc)) {
      diffFun <- function(densities, meas.list) {
        sapply(seq_len(ncol(meas.list[[1]])), function(x)
               auc_diff(densities, cbind(meas.list[[1]][, x], meas.list[[2]][, x])))
      }
    } else {
      diffFun <- function(densities, meas.list) cbind(densities, meas.list[[1]] - meas.list[[2]])
    }
  }

  if (is.null(perms)) perms <- shuffleSet(n=nrow(resids$resids.all), nset=N)
  N <- nrow(perms)
  perms <- rbind(perms, 1:ncol(perms))  # last row is observed metrics
  groups <- as.numeric(resids$resids.all$Group)

  # Loop through the permutation matrix
  res.perm <- switch(level,
               vertex=permute_vertex_foreach(perms, densities, resids, groups, measure, diffFun),
               other=permute_other_foreach(perms, densities, resids, groups, .function),
               graph=permute_graph_foreach(perms, densities, resids, groups, atlas, auc))

  if (length(densities) == 1) res.perm <- cbind(densities=densities, res.perm)

  if (level == 'vertex') {
    res.perm <- as.data.table(res.perm)
    regions <- names(resids$resids.all[, !c('Study.ID', 'Group')])
    if (isTRUE(auc)) {
      setnames(res.perm, 1:ncol(res.perm), regions)
    } else {
      setnames(res.perm, 2:ncol(res.perm), regions)
    }
  }

  if (!isTRUE(auc)) {
    setkey(res.perm, densities)
    obs.ind <- (N + 1) * 1:length(densities)
    obs.diff <- res.perm[obs.ind]
    res.perm <- res.perm[-obs.ind]
  } else {
    obs.diff <- res.perm[.N]
    res.perm <- res.perm[-.N]
  }
  out <- list(atlas=atlas, auc=auc, N=N, level=level, measure=measure, densities=densities,
              resids=resids, DT=res.perm, obs.diff=obs.diff, groups=resids$groups)
  class(out) <- c('brainGraph_permute', class(out))
  return(out)
}

#==============================================================================
# Helper functions
#==============================================================================
make_graphs_perm <- function(densities, resids, inds, groups) {
  corrs <- lapply(unique(groups), function(x)
                  corr.matrix(resids[which(groups[inds] == x)],
                              densities=densities, rand=TRUE))
  sapply(corrs, lapply, function(x)
         apply(x$r.thresh, 3, graph_from_adjacency_matrix, mode='undirected', diag=F))
}

# Graph level
diffFun_graph_noAUC <- function(densities, meas.list) {
  tmp <- data.table(densities=densities, meas=meas.list[, 1] - meas.list[, 2], key='densities')
  setnames(tmp, 'meas', deparse(substitute(meas.list)))
  return(tmp)
}
graph_attr_perm <- function(g, densities, atlas) {
  g <- lapply(g, lapply, make_brainGraph, atlas, type='random')

  mod <- sapply(g, sapply, function(x) modularity(cluster_louvain(x)))
  Cp <- sapply(g, sapply, function(x) transitivity(x, type='localaverage'))
  Lp <- sapply(g, sapply, mean_distance)
  assort <- sapply(g, sapply, assortativity_degree)
  E.global <- sapply(g, sapply, efficiency, 'global')
  assort.lobe <- sapply(g, sapply, function(x)
                        assortativity_nominal(x, as.integer(factor(V(x)$lobe))))
  asymm <- sapply(g, sapply, function(x) edge_asymmetry(x)$asymm)
  list(mod=mod, Cp=Cp, Lp=Lp, assort=assort, E.global=E.global, assort.lobe=assort.lobe, asymm=asymm)
}
graph_attr_perm_diffs <- function(densities, meas.list, auc) {
  if (isTRUE(auc)) {
    tmp <- data.table(t(sapply(meas.list, function(y) auc_diff(densities, y))))
  } else {
    meas.dt <- lapply(meas.list, function(y) diffFun_graph_noAUC(densities, y))
    for (i in 1:7) setnames(meas.dt[[i]], 'y', names(meas.dt)[i])
    tmp <- Reduce(merge, meas.dt)
  }
  return(tmp)
}
permute_graph_foreach <- function(perms, densities, resids, groups, atlas, auc) {
  i <- NULL
  res.perm <- foreach(i=seq_len(nrow(perms)), .combine='rbind') %dopar% {
    g <- make_graphs_perm(densities, resids, perms[i, ], groups)
    meas.list <- graph_attr_perm(g, densities, atlas)
    graph_attr_perm_diffs(densities, meas.list, auc)
  }
}

# Vertex-level
vertex_attr_perm <- function(measure, g, densities) {
  switch(measure,
    vulnerability=lapply(g, function(x) t(sapply(x, vulnerability))),
    degree=lapply(g, function(x) t(sapply(x, degree))),
    E.nodal=lapply(g, function(x) t(sapply(x, efficiency, 'nodal'))),
    ev.cent=lapply(g, function(x) t(sapply(x, function(y) centr_eigen(y)$vector))),
    knn=lapply(g, function(x) t(sapply(x, function(y) graph.knn(y)$knn))),
    transitivity=lapply(g, function(x) t(sapply(x, transitivity,type='local', isolates='zero'))),
    lapply(g, function(x) t(sapply(x, function(y) centr_betw(y)$res))))
}
permute_vertex_foreach <- function(perms, densities, resids, groups, measure, diffFun) {
  i <- NULL
  res.perm <- foreach(i=seq_len(nrow(perms)), .combine='rbind') %dopar% {
    g <- make_graphs_perm(densities, resids, perms[i, ], groups)
    meas.list <- vertex_attr_perm(measure, g, densities)
    diffFun(densities, meas.list)
  }
}

# Other-level
permute_other_foreach <- function(perms, densities, resids, groups, .function) {
  i <- NULL
  res.perm <- foreach(i=seq_len(nrow(perms)), .combine='rbind') %dopar% {
    g <- make_graphs_perm(densities, resids, perms[i, ], groups)
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
#' @method summary brainGraph_permute
#' @rdname brainGraph_permute

summary.brainGraph_permute <- function(object, measure=NULL,
                                       alternative=c('two.sided', 'less', 'greater'),
                                       alpha=0.05, p.sig=c('p', 'p.fdr'), ...) {
  perm.diff <- p <- N <- p.fdr <- region <- obs.diff <- NULL

  permDT <- copy(object$DT)
  g <- with(object, make_graphs_perm(densities, resids, 1:nrow(resids$resids.all),
                                     resids$resids.all[, as.numeric(Group)]))
  # OTHER
  #-------------------------------------
  if (object$level == 'other') {  # Hack to figure out which level it is when level="other"
    if (ncol(permDT) > 8) {
      object$level <- 'vertex'
    } else {
      object$level <- 'graph'
    }
  }

  # VERTEX-LEVEL
  #-------------------------------------
  if (object$level == 'vertex') {
    densities <- object$densities
    measure <- object$measure
    obsDT <- copy(object$obs.diff)
    meas.list <- with(object, vertex_attr_perm(measure, g, densities))

    if (isTRUE(object$auc)) {
      obs <- lapply(meas.list, apply, 2, function(y) sum(diff(object$densities) * (head(y, -1) + tail(y, -1))) / 2)
      permDT[, densities := 1]
      obsDT[, densities := 1]
      densities <- 1
    } else {
      obs <- meas.list
    }
    sum.dt <- data.table(densities=densities,
                         region=rep(V(g[[1]][[1]])$name, each=length(densities)),
                         g1=c(obs[[1]]),
                         g2=c(obs[[2]]),
                         key='densities')
    setnames(sum.dt, c('g1', 'g2'), paste0(measure, '.', object$groups))
    obs <- melt(obsDT, id.vars='densities', variable.name='region', value.name='obs.diff')
    sum.dt <- merge(sum.dt, obs, by=c('densities', 'region'))
    permDT <- melt(permDT, id.vars='densities', variable.name='region', value.name=measure)
    setkeyv(permDT, key(sum.dt))
    permdiff <- permDT[, list(perm.diff=mean(get(measure))), by=list(densities, region)]
    sum.dt <- merge(sum.dt, permdiff, by=key(sum.dt))

  # GRAPH-LEVEL
  #-------------------------------------
  } else if (object$level == 'graph') {
    if (is.null(measure)) measure <- 'mod'
    stopifnot(measure %in% names(permDT))
    permDT[, region := 'graph']
    if (measure %in% c('asymm', 'assortativity.lobe')) {
      g <- lapply(g, lapply, make_brainGraph, object$atlas, type='random')
    }
    meas.list <- with(object, graph_attr_perm(g, densities, atlas))

    obs <- meas.list[[measure]]
    if (isTRUE(object$auc)) {
      obs <- apply(obs, 2, function(y) sum(diff(object$densities) * (head(y, -1) + tail(y, -1))) / 2)
      sum.dt <- data.table(densities=1, region='graph')
      sum.dt[, (paste0(measure, '.', object$groups)) := as.list(obs)]
      sum.dt[, obs.diff := object$obs.diff[[measure]]]
      sum.dt[, perm.diff := permDT[, mean(get(measure))]]

      permDT[, densities := 1]

    } else {
      sum.dt <- data.table(densities=object$densities, region='graph')
      for (i in seq_along(object$groups)) {
        sum.dt[, paste0(measure, '.', object$groups[i]) := obs[, i]]
      }
      sum.dt[, obs.diff := object$obs.diff[[measure]]]
      sum.dt[, perm.diff := permDT[, mean(get(measure)), by=list(densities, region)]$V1]
    }
  }
  result.dt <- merge(permDT[, c('densities', 'region', measure), with=F],
                     sum.dt[, c('densities', 'region', 'obs.diff'), with=F],
                     by=c('densities', 'region'))

  alt <- match.arg(alternative)
  if (alt == 'two.sided') {
    result.dt[, p := (sum(abs(get(measure)) >= abs(unique(obs.diff))) + 1) / (.N + 1), by=key(result.dt)]
    CI <- c(alpha / 2, 1 - (alpha / 2))
  } else if (alt == 'less') {
    result.dt[, p := (sum(get(measure) <= unique(obs.diff)) + 1) / (.N + 1), by=key(result.dt)]
    CI <- c(alpha, 1)
  } else if (alt == 'greater') {
    result.dt[, p := (sum(get(measure) >= unique(obs.diff)) + 1) / (.N + 1), by=key(result.dt)]
    CI <- c(1 / N, 1 - alpha)
  }
  result.dt[, c('ci.low', 'ci.high') := as.list(sort(get(measure))[ceiling(.N * CI)]), by=key(result.dt)]
  result.dt <- result.dt[, .SD[1], by=key(result.dt)]
  sum.dt <- merge(sum.dt, result.dt[, !c(measure, 'obs.diff'), with=F], by=key(result.dt))
  setcolorder(sum.dt,
              c('densities', 'region', paste0(measure, '.', object$groups), 'obs.diff',
                'ci.low', 'ci.high', 'perm.diff', 'p'))
  if (!isTRUE(object$auc)) {
    if (object$level == 'graph') {
      sum.dt[, p.fdr := p.adjust(p, 'fdr')]
    } else {
      sum.dt[, p.fdr := p.adjust(p, 'fdr'), by=densities]
    }
  }

  meas.full <- switch(measure,
                      mod='Modularity',
                      E.global='Global efficiency',
                      E.global.wt='Global efficiency (weighted)',
                      Cp='Clustering coefficient',
                      Lp='Characteristic path length',
                      assort='Degree assortativity',
                      assort.lobe='Lobe assortativity',
                      asymm='Edge asymmetry',
                      btwn.cent='Betweenness centrality',
                      vulnerability='Vulnerability',
                      degree='Degree',
                      E.nodal='Nodal efficiency',
                      ev.cent='Eigenvector centrality',
                      knn='K-nearest neighbor degree',
                      transitivity='Local transitivity')

  p.sig <- match.arg(p.sig)
  perm.sum <- with(object, list(auc=auc, N=N, level=level, densities=densities,
                                DT.sum=sum.dt, meas.full=meas.full, groups=groups,
                                alt=alt, alpha=alpha, p.sig=p.sig))
  class(perm.sum) <- c('summary.brainGraph_permute', class(perm.sum))
  perm.sum
}

#' @aliases summary.brainGraph_permute
#' @method print summary.brainGraph_permute

print.summary.brainGraph_permute <- function(x, ...) {
  print_title_summary('Permutation analysis')
  cat('# of permutations:', prettyNum(x$N, ','), '\n')
  cat('Level: ', x$level, '\n')
  cat('Graph metric: ', x$meas.full, '\n')
  if (isTRUE(x$auc)) cat('Area-under-the-curve (AUC) calculated across', length(x$densities), 'densities:\n', x$densities, '\n')

  alt <- switch(x$alt,
                two.sided=with(x, sprintf('%s - %s != 0', groups[1], groups[2])),
                greater=with(x, sprintf('%s - %s > 0', groups[1], groups[2])),
                less=with(x, sprintf('%s - %s < 0', groups[1], groups[2])))
  cat('Alternative hypothesis: ', alt, '\n')
  cat('Alpha: ', x$alpha, '\n\n')
  if (with(x, nrow(DT.sum[get(p.sig) < alpha])) == 0) {
    cat ('No significant results!\n')
  } else {
    with(x, print(DT.sum[get(p.sig) < alpha]))
  }
  invisible(x)
}

#' Plot results from permutation testing
#'
#' @param ptitle Character string specifying a title for the plot (default:
#'   \code{NULL})
#' @export
#' @method plot brainGraph_permute
#' @return The \code{plot} method returns a \emph{list} of \code{ggplot} objects
#' @rdname brainGraph_permute

plot.brainGraph_permute <- function(x, measure=NULL,
                                    alternative=c('two.sided', 'less', 'greater'),
                                    alpha=0.05, p.sig=c('p', 'p.fdr'), ptitle=NULL, ...) {
  densities <- Group <- sig <- trend <- yloc <- obs <- mylty <- ci.low <- ci.high <-
    variable <- value <- reg.num <- region <- perm.diff <- obs.diff <- NULL
  p.sig <- match.arg(p.sig)
  if (x$level == 'graph') {
    if (is.null(measure)) measure <- 'mod'
    stopifnot(measure %in% names(x$DT))
  } else {
    measure <- x$measure
  }
  perm.sum <- summary(x, measure=measure, alternative=alternative, alpha=alpha)
  sum.dt <- perm.sum$DT.sum
  if (is.null(ptitle)) ptitle <- perm.sum$meas.full
  ylabel2 <- sprintf('Observed and permutation difference (%s - %s)', x$groups[1], x$groups[2])

  # GRAPH LEVEL
  #-------------------------------------
  if (x$level == 'graph') {
    plot.dt <- melt(sum.dt, id.vars=setdiff(names(sum.dt), paste0(measure, '.', x$groups)),
                    variable.name='Group', value.name='obs')
    plot.dt[, Group := factor(Group, labels=x$groups)]
    idvars <- c('densities', 'region', 'p', 'Group', 'obs')
    if (!isTRUE(x$auc)) idvars <- c(idvars, 'p.fdr')
    plot.dt <- melt(plot.dt, id.vars=idvars)
    plot.dt[, c('sig', 'trend') := '']
    plot.dt[get(p.sig) < alpha, sig := '*']
    plot.dt[get(p.sig) >= alpha & get(p.sig) < 2 * alpha, trend := '*']
    plot.dt[, yloc := round(min(obs) - 0.05 * diff(range(obs)), 3)]
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
      geom_line(data=plot.dt[variable=='obs.diff'], aes(y=value, lty=factor(mylty)), col='red') +
      geom_point(data=plot.dt[variable == 'obs.diff'], aes(y=value), col='red', size=3) +
      geom_line(data=plot.dt[variable == 'perm.diff'], aes(y=value, lty=factor(mylty))) +
      geom_line(data=plot.dt[variable == 'ci.low'], aes(y=value, lty=factor(mylty))) +
      geom_line(data=plot.dt[variable == 'ci.high'], aes(y=value, lty=factor(mylty))) +
      scale_linetype_manual(name='',
                            labels=c('  Observed diff.     ',
                                     '  95% CI     ',
                                     '  Mean permutation diff.     '),
                            values=1:3) +
      theme(legend.position='bottom',
            legend.key=element_rect(fill='white'),
            legend.background=element_rect(fill='gray79'))

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
    if (nrow(sum.dt[get(p.sig) < alpha]) == 0) stop('No significant results!')
    plot.dt <- droplevels(sum.dt[get(p.sig) < alpha])
    plot.dt[, reg.num := seq_len(.N), by=densities]
    p <- ggplot(plot.dt[get(p.sig) < alpha], aes(x=region))
    plot1 <- p +
      geom_point(aes(y=perm.diff)) +
      geom_errorbar(aes(ymin=ci.low, ymax=ci.high)) +
      geom_segment(aes(x=reg.num - .25, xend=reg.num + .25, y=obs.diff, yend=obs.diff),
                   col='red', size=1.25) +
      facet_wrap(~ densities, scales='free') +
      labs(title=ptitle, y=ylabel2, x='Region') +
      theme(plot.title=element_text(hjust=0.5, face='bold'))
    plot2 <- NULL
  }
  return(list(plot1, plot2))
}
