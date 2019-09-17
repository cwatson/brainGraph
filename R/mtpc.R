#' Multi-threshold permutation correction
#'
#' Applies the \emph{multi-threshold permutation correction (MTPC)} method to
#' perform inference in graph theory analyses of brain MRI data.
#'
#' This is a multi-step procedure: (steps 3-4 are the time-consuming steps)
#' \enumerate{
#'   \item Apply thresholds \eqn{\tau} to the networks, and compute network
#'     metrics for all networks and thresholds. (already done beforehand)
#'   \item Compute test statistics \eqn{S_{obs}} for each threshold. (done by
#'     \code{\link{brainGraph_GLM}})
#'   \item Permute group assignments and compute test statistics for each
#'     permutation and threshold. (done by \code{\link{brainGraph_GLM}})
#'   \item Build a null distribution of the maximum statistic across thresholds
#'     (and across brain regions) for each permutation. (done by
#'     \code{\link{brainGraph_GLM}})
#'   \item Determine the critical value, \eqn{S_{crit}} from the null
#'     distribution of maximum statistics.
#'   \item Identify clusters where \eqn{S_{obs} > S_{crit}} and compute the AUC
#'     for these clusters (denoted \eqn{A_{MTPC}}).
#'   \item Compute a critical AUC (\eqn{A_{crit}}) from the mean of the
#'     supra-critical AUC's for the permuted tests.
#'   \item Reject \eqn{H_0} if \eqn{A_{MTPC} > A_{crit}}.
#' }
#'
#' @param g.list A list of \code{brainGraphList} objects for all thresholds
#' @param thresholds Numeric vector of the thresholds applied to the raw
#'   connectivity matrices.
#' @param clust.size Integer indicating the size of \dQuote{clusters} (i.e.,
#'   consecutive thresholds for which the observed statistic exceeds the null)
#'   (default: \code{3L})
#' @param res.glm A list of \code{bg_GLM} objects, as output by a previous run
#'   of \code{mtpc}. Useful if you want to change the cluster size without
#'   re-running all of the GLM's and permutations (default: \code{NULL})
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM}} and/or
#'   \code{\link{brainGraph_GLM_design}}
#' @inheritParams GLM
#' @export
#'
#' @return An object of class \code{mtpc} with some input arguments plus the
#'   following elements:
#'   \item{res.glm}{List with length equal to the number of thresholds; each
#'     list element is the output from \code{\link{brainGraph_GLM}}}
#'   \item{DT}{A \code{data.table} for all thresholds, combined from the outputs
#'     of \code{\link{brainGraph_GLM}}}
#'   \item{stats}{A data.table containing \code{S.mtpc} (the max. observed
#'     statistic), \code{tau.mtpc} (the threshold of the max. observed
#'     statistic), \code{S.crit} (the critical statistic value), and
#'     \code{A.crit} (the critical AUC)}
#'   \item{null.dist}{Numeric matrix with \code{N} rows and number of columns
#'     equal to the number of thresholds. Each element is the maximum statistic
#'     for that permutation and threshold.}
#'   \item{perm.order}{Numeric matrix; the permutation set applied for all
#'     thresholds (each row is a separate permutation)}
#'
#' @family Group analysis functions
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Drakesmith, M. and Caeyenberghs, K. and Dutt, A. and Lewis, G. and
#'   David, A.S. and Jones, D.K. (2015) Overcoming the effects of false
#'   positives and threshold bias in graph theoretical analyses of neuroimaging
#'   data. \emph{NeuroImage}, \bold{118}, 313--333.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2015.05.011}
#' @examples
#' \dontrun{
#' diffs.mtpc <- mtpc(g.list=g.norm, thresholds=thresholds, N=N,
#'      covars=covars.dti, measure='E.nodal.wt', coding='effects',
#'      contrasts=c(0, 0, 0, 0, -2), alt='greater',
#'      binarize=c('Sex', 'Scanner'), con.name='Group 1 > Group 2')
#' sig.regions <- diffs.mtpc$DT[A.mtpc > A.crit]
#' }

mtpc <- function(g.list, thresholds, covars, measure, contrasts, con.type=c('t', 'f'),
                 outcome=NULL, con.name=NULL, level=c('vertex', 'graph'),
                 clust.size=3L, perm.method=c('freedmanLane', 'terBraak', 'smith'),
                 part.method=c('beckmann', 'guttman', 'ridgway'), N=500L, perms=NULL,
                 alpha=0.05, res.glm=NULL, long=TRUE, ...) {
  A.crit <- A.mtpc <- contrast <- V1 <- S.crit <- DT <- region <- stat <- threshold <- values <- null.out <- NULL

  # Check if components are 'brainGraphList' objects
  matches <- vapply(g.list, inherits, logical(1), 'brainGraphList')
  if (any(!matches)) stop("Input must be a list of 'brainGraphList' objects.")
  stopifnot(length(g.list) == length(thresholds))

  if (is.null(perms)) perms <- shuffleSet(n=length(g.list[[1]]$graphs), nset=N)

  if (is.null(res.glm)) {
    res.glm <- lapply(g.list, function(z)
      brainGraph_GLM(z, covars, measure, contrasts, con.type, outcome, con.name=con.name, N=N,
                     level=level, permute=TRUE, perm.method=perm.method, part.method=part.method,
                     perms=perms, alpha=alpha, long=TRUE, ...))
  }
  alt <- res.glm[[1]]$alt
  for (i in seq_along(thresholds)) res.glm[[i]]$DT[, threshold := thresholds[i]]
  mtpc.all <- rbindlist(lapply(res.glm, with, DT))
  setkey(mtpc.all, contrast, region)
  mtpc.all[, c('S.crit', 'A.mtpc', 'A.crit') := 0]

  con.type <- match.arg(con.type)
  kNumContrasts <- if (con.type == 't') nrow(res.glm[[1]]$contrasts) else length(res.glm[[1]]$contrasts)
  null.dist.all <- null.dist.max <- vector('list', length=kNumContrasts)
  Scrit <- Acrit <- rep(0, kNumContrasts)
  myMax <- maxfun(alt)
  mySort <- sortfun(alt)
  for (i in seq_len(kNumContrasts)) {
    null.dist.all[[i]] <- vapply(res.glm, function(x) x$perm$null.dist[contrast == i, V1], numeric(N))
    null.dist.max[[i]] <- apply(null.dist.all[[i]], 1, myMax)
    Scrit[i] <- mySort(null.dist.max[[i]])[floor(N * (1 - alpha)) + 1]
    mtpc.all[contrast == i, S.crit := Scrit[i]]
  }

  rlefun <- switch(alt,
                   two.sided=function(x, y) rle(abs(x) > abs(y)),
                   less=function(x, y) rle(x < y),
                   greater=function(x, y) rle(x > y))
  crit.regions <- mtpc.all[, rlefun(stat, S.crit), by=key(mtpc.all)][values == TRUE & lengths >= clust.size, list(region=unique(region)), by=contrast]
  setkey(crit.regions, contrast, region)
  mtpc.all[crit.regions,
           A.mtpc := auc_rle(stat, S.crit, thresholds, alt, clust.size),
           by=list(contrast, region)]

  maxobsfun <- switch(alt,
                      two.sided=function(x) max(is.finite(abs(x)) * abs(x), na.rm=TRUE),
                      less=function(x) min(is.finite(x) * x, na.rm=TRUE),
                      greater=function(x) max(is.finite(x) * x, na.rm=TRUE))
  S.mtpc <- mtpc.all[, maxobsfun(stat), by=contrast]
  S.mtpc <- S.mtpc[, unique(V1), by=contrast]
  tau.mtpc <- mtpc.all[stat %in% S.mtpc$V1, threshold, by=contrast]
  tau.mtpc <- tau.mtpc[, .SD[1], by=contrast]

  for (i in seq_len(kNumContrasts)) {
    null.crit <- which(apply(null.dist.all[[i]], 1, function(x)
                             any(with(rlefun(x, Scrit[i]), values == TRUE & lengths >= clust.size))))
    if (length(null.crit) > 0) {
      null.crits <- null.dist.all[[i]][null.crit, , drop=FALSE]
      Acrit[i] <- mean(apply(null.crits, 1, auc_rle, Scrit[i], thresholds, alt, clust.size))
      mtpc.all[contrast == i, A.crit := Acrit[i]]
    }
  }

  glm.attr <- res.glm[[1]]
  glm.attr[c('y', 'DT', 'permute', 'perm')] <- NULL
  for (i in seq_along(thresholds)) res.glm[[i]]$perm$null.dist <- NULL
  mtpc.stats <- data.table(contrast=seq_len(kNumContrasts), tau.mtpc=tau.mtpc$threshold,
                           S.mtpc=S.mtpc$V1, S.crit=Scrit, A.crit=Acrit)

  if (isTRUE(long)) null.out <- null.dist.all
  out <- c(glm.attr, list(res.glm=res.glm, clust.size=clust.size, DT=mtpc.all, stats=mtpc.stats, null.dist=null.out, perm.order=perms,
                          perm.method=glm.attr$perm.method, part.method=glm.attr$part.method))
  class(out) <- c('mtpc', class(out))
  return(out)
}

auc_rle <- function(t.stat, S.crit, thresholds, alt, clust.size) {
  inds <- get_rle_inds(clust.size, alt, t.stat, S.crit, thresholds)
  x <- Map(function(a, b) thresholds[a:b], inds$starts, inds$ends)
  y <- Map(function(a, b) t.stat[a:b], inds$starts, inds$ends)
  auc <- mapply(auc_diff, x, y)
  return(sum(abs(auc) * is.finite(auc), na.rm=TRUE))
}

get_rle_inds <- function(clust.size, alt, t.stat, S.crit, thresholds) {
  compfun <- switch(alt,
                    two.sided=function(x, y) rle(abs(x) > abs(y)),
                    less=function(x, y) rle(x < y),
                    greater=function(x, y) rle(x > y))
  runs <- compfun(t.stat, S.crit)
  myruns <- which(runs$values == TRUE & runs$lengths >= clust.size)
  runs.len.cumsum <- cumsum(runs$lengths)
  ends <- runs.len.cumsum[myruns]
  newindex <- ifelse(myruns > 1, myruns - 1, 0)
  starts <- runs.len.cumsum[newindex] + 1
  if (0 %in% newindex) starts <- c(1, starts)
  return(list(starts=starts, ends=ends))
}

################################################################################
# S3 METHODS
################################################################################

#' Print a summary of MTPC results
#'
#' @param object,x A \code{mtpc} object
#' @inheritParams summary.bg_GLM
#' @export
#' @method summary mtpc
#' @rdname mtpc

summary.mtpc <- function(object, contrast=NULL, digits=max(3L, getOption('digits') - 2L), print.head=TRUE, ...) {
  A.mtpc <- A.crit <- stat <- region <- S.mtpc <- NULL
  object$printCon <- contrast
  object$digits <- digits
  object$print.head <- print.head

  # Summary table
  whichmaxfun <- switch(object$alt, two.sided=function(x) which.max(abs(x)), less=which.min, greater=which.max)
  DT.sum <- object$DT[A.mtpc > A.crit, .SD[whichmaxfun(is.finite(stat) * stat)], by=list(contrast, region)]
  DT.sum <- DT.sum[, c('region', 'contrast', 'stat', 'Outcome', 'Contrast', 'threshold', 'S.crit', 'A.mtpc', 'A.crit'), with=FALSE]
  setcolorder(DT.sum, c('Outcome', 'Contrast', 'region', 'threshold', 'stat', 'S.crit', 'A.mtpc', 'A.crit', 'contrast'))
  setnames(DT.sum, c('region', 'stat', 'threshold'), c('Region', 'S.mtpc', 'tau.mtpc'))
  setorder(DT.sum, contrast, -S.mtpc)
  object$DT.sum <- DT.sum

  class(object) <- c('summary.mtpc', class(object))
  return(object)
}

#' @aliases summary.mtpc
#' @method print summary.mtpc
#' @keywords internal

print.summary.mtpc <- function(x, ...) {
  A.mtpc <- A.crit <- contrast <- region <- NULL
  print_title_summary('MTPC results')
  cat('Level: ', x$level, '\n')

  print_measure_summary(x)
  print_permutation_summary(x)

  cat('# of thresholds: ', length(x$res.glm), '\n')
  cat('Cluster size (across thresholds): ', x$clust.size, '\n\n')

  print_contrast_type_summary(x)
  print_subs_summary(x)

  message('\n', 'Statistics', '\n', rep('-', 20))
  clp <- 100 * (1 - x$alpha)
  cat('tau.mtpc: threshold of the maximum observed statistic\n')
  cat('S.mtpc: maximum observed statistic\n')
  cat('S.crit: the critical', paste0('(', clp, 'th percentile)'), 'value of the null max. statistic\n')
  cat('A.crit: critical AUC from the supra-critical null AUCs\n\n')
  print(x$stats)
  cat('\n')

  # Print results for each contrast
  print_contrast_stats_summary(x)

  invisible(x)
}

#' Plot statistics from an MTPC analysis
#'
#' Plot the statistics from an MTPC analysis, along with the maximum permuted
#' statistics. The output is similar to Figure 11 in Drakesmith et al. (2015).
#'
#' @param only.sig.regions Logical indicating whether to plot only significant
#'   regions (default: \code{TRUE})
#' @param show.null Logical indicating whether to plot points of the maximum
#'   null statistics (per permutation)
#' @param caption.stats Logical indicating whether to print the MTPC statistics
#'   in the caption of the plot (default: \code{FALSE})
#' @inheritParams plot.bg_GLM
#' @export
#' @method plot mtpc
#' @rdname mtpc
#'
#' @return The \code{plot} method returns a \emph{list} of
#'   \code{\link[ggplot2]{ggplot}} objects
#' @examples
#' \dontrun{
#' mtpcPlots <- plot(mtpc.diffs)
#'
#' ## Arrange plots into 3x3 grids
#' ml <- marrangeGrob(mtpcPlots, nrow=3, ncol=3)
#' ggsave('mtpc.pdf', ml)
#' }

plot.mtpc <- function(x, contrast=1L, region=NULL, only.sig.regions=TRUE,
                      show.null=TRUE, caption.stats=FALSE, ...) {
  stat_ribbon <- stat <- threshold <- S.crit <- A.mtpc <- A.crit <- NULL

  stopifnot(inherits(x, 'mtpc'))
  mycontrast <- contrast
  DT <- x$DT[contrast == mycontrast, !'p.fdr']
  DT[, stat_ribbon := stat]  # For filling in supra-threshold areas
  thresholds <- DT[, unique(threshold)]
  DT$nullthresh <- x$stats[contrast == mycontrast, unique(S.crit)]
  whichmaxfun <- switch(x$alt, two.sided=function(y) which.max(abs(y)), less=which.min, greater=which.max)
  myMax <- maxfun(x$alt)
  thr <- apply(x$null.dist[[mycontrast]], 1, whichmaxfun)
  thr.y <- apply(x$null.dist[[mycontrast]], 1, myMax)
  nullcoords <- data.table(threshold=thresholds[thr], y=thr.y)

  # Local function to plot for a single region
  plot_single <- function(x, DT, nullcoords, show.null) {
    stat <- S.crit <- threshold <- stat_ribbon <- nullthresh <- y <- A.mtpc <- A.crit <- NULL

    inds <- DT[, get_rle_inds(x$clust.size, x$alt, stat, S.crit, threshold)]
    n <- dim(inds)[1L]
    if (n == 0) {
      DT[, stat_ribbon := NA]
    } else {
      if (n == 1) {
        all.inds <- c(apply(inds, 1, function(z) z[1]:z[2]))
      } else if (n > 1) {
        all.inds <- c(unlist(apply(inds, 1, function(z) z[1]:z[2])))
      } else {
        all.inds <- 1:dim(DT)[1L]
      }
      DT[-all.inds, stat_ribbon := NA]
    }

    lineplot <- ggplot(data=DT, mapping=aes(x=threshold)) +
      geom_line(aes(y=stat), col='red4', size=1.25, na.rm=TRUE) +
      geom_hline(aes(yintercept=nullthresh), lty=2)
    if (n > 0) {
      if (x$alt == 'less') {
        lineplot <- lineplot +
          geom_ribbon(aes(ymax=stat_ribbon, ymin=nullthresh), fill='red4', alpha=0.3)
      } else {
        lineplot <- lineplot +
          geom_ribbon(aes(ymin=stat_ribbon, ymax=nullthresh), fill='red4', alpha=0.3)
      }
      lineplot <- lineplot + geom_point(aes(y=stat_ribbon), col='red4', size=2, na.rm=TRUE)
    }

    if (isTRUE(all(diff(thresholds) < 0))) lineplot <- lineplot + scale_x_reverse()
    p.region <- if (x$level == 'graph') 'graph-level' else DT[, as.character(unique(region))]
    p.title <- paste0('Region: ', p.region)
    p.subtitle <- paste0('Outcome: ', x$outcome)

    if (isTRUE(show.null)) {
      lineplot <- lineplot + geom_point(data=nullcoords, aes(y=y), col='darkgreen', alpha=0.4, na.rm=TRUE)
    }
    if (isTRUE(caption.stats)) {
      Smtpc <- bquote('S'['mtpc']*' = '~ .(DT[, format(myMax(stat))]))
      Scrit <- bquote('S'['crit']*' = '~ .(DT[, format(unique(S.crit))]))
      Amtpc <- bquote('A'['mtpc']*' = '~ .(DT[, format(max(A.mtpc))]))
      Acrit <- bquote('A'['crit']*' = '~ .(DT[, format(unique(A.crit))]))
      statslabel <- bquote(atop(.(Scrit)~ ';'~ .(paste('\t'))~ .(Acrit),
                                .(Smtpc)~ ';'~ .(paste('\t'))~ .(Amtpc)))

      if (DT[, unique(A.mtpc) > unique(A.crit)]) statslabel <- bquote(.(statslabel)~ .(paste0('\t(p < ', x$alpha, ')')))
      lineplot <- lineplot +
        labs(caption=statslabel) +
        theme(plot.caption=element_text(hjust=0))
    }
    lineplot <- lineplot +
      labs(title=p.title, subtitle=p.subtitle, x='Threshold', y=paste0(toupper(x$con.type), '-statistic')) +
      theme(plot.title=element_text(hjust=0.5, face='bold'),
            plot.subtitle=element_text(hjust=0.5))

    return(lineplot)
  }

  if (x$level == 'graph') {
    lineplot <- plot_single(x, DT, nullcoords, show.null)
    return(lineplot)
  } else if (x$level == 'vertex') {
    if (is.null(region)) {
      if (!isTRUE(only.sig.regions)) {
        region <- DT[, levels(region)]
      } else {
        region <- droplevels(DT[A.mtpc > A.crit])[, levels(region)]
      }
    }

    lineplots <- setNames(vector('list', length(region)), region)
    for (z in region) {
      lineplots[[z]] <- plot_single(x, DT[region == z], nullcoords, show.null)
    }
    return(lineplots)
  }
}
