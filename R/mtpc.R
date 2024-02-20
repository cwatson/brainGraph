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
#' @importFrom permute shuffleSet
#'
#' @return An object of class \code{mtpc} with some input arguments plus the
#'   following elements:
#'   \item{X,qr,cov.unscaled}{Design matrix, QR decomposition, and unscaled
#'     covariance matrix, if the design is the same across thresholds}
#'   \item{contrasts}{The contrast matrix or list of matrices}
#'   \item{con.name}{Contrast names}
#'   \item{removed.subs}{Named integer vector of subjects with incomplete data}
#'   \item{atlas}{The atlas of the input graphs}
#'   \item{rank,df.residual}{The model rank and residual degrees of freedom}
#'   \item{res.glm}{List with length equal to the number of thresholds; each
#'     list element is the output from \code{\link{brainGraph_GLM}}}
#'   \item{DT}{A \code{data.table} for all thresholds, combined from the outputs
#'     of \code{\link{brainGraph_GLM}}}
#'   \item{stats}{A data.table containing \code{S.mtpc} (the max. observed
#'     statistic), \code{tau.mtpc} (the threshold of the max. observed
#'     statistic), \code{S.crit} (the critical statistic value), and
#'     \code{A.crit} (the critical AUC)}
#'   \item{null.dist}{Numeric array with \code{N} columns and number of rows
#'     equal to the number of thresholds. The 3rd dimension is for each
#'     contrast. Each element of the array is the maximum statistic
#'     for that permutation, threshold, and contrast combination.}
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
#'   \doi{10.1016/j.neuroimage.2015.05.011}
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
                 clust.size=3L, perm.method=c('freedmanLane', 'terBraak', 'smith',
                                              'draperStoneman', 'manly', 'stillWhite'),
                 part.method=c('beckmann', 'guttman', 'ridgway'), N=500L, perms=NULL,
                 alpha=0.05, res.glm=NULL, long=TRUE, ...) {
  S.mtpc <- tau.mtpc <- A.crit <- A.mtpc <- Contrast <- S.crit <- region <- stat <-
    threshold <- values <- null.out <- NULL

  # Check if components are 'brainGraphList' objects
  matches <- vapply(g.list, is.brainGraphList, logical(1L))
  if (any(!matches)) stop("Input must be a list of 'brainGraphList' objects.")
  stopifnot(length(g.list) == length(thresholds))

  if (is.null(res.glm)) {
    # If 'outcome != measure' make sure the same perm order is used throughout
    if (is.null(perms)) perms <- shuffleSet(n=length(g.list[[1L]]$graphs), nset=N)
    if (!is.null(outcome) && (outcome != measure)) {
      n <- sum(complete.cases(covars))
      if (dim(perms)[1L] != n) perms <- shuffleSet(n=n, nset=N)
    }
    res.glm <- lapply(g.list, function(z)
      brainGraph_GLM(z, covars, measure, contrasts, con.type, outcome, con.name=con.name, N=N,
                     level=level, permute=TRUE, perm.method=perm.method, part.method=part.method,
                     perms=perms, alpha=alpha, long=TRUE, ...))
  }
  alt <- res.glm[[1L]]$alt
  for (i in seq_along(thresholds)) res.glm[[i]]$DT[, threshold := thresholds[i]]
  mtpc.all <- rbindlist(lapply(res.glm, function(x) x$DT))
  setkey(mtpc.all, Contrast, region)
  mtpc.all[, c('S.crit', 'A.mtpc', 'A.crit') := 0]

  conNames <- res.glm[[1L]]$con.name
  myMax <- switch(alt, two.sided=colMaxAbs, less=colMin, greater=colMax)
  mySort <- sortfun(alt)
  index <- floor((1 - alpha) * N) + 1L
  # Dimensions will be "thresholds X Nperms X contrasts"
  null.max.all <- aperm(sapply(res.glm, function(x) x$perm$null.max, simplify='array'),
                        c(3L, 1L, 2L))
  dimnames(null.max.all) <- list(NULL, NULL, conNames)
  null.max.perms <- apply(null.max.all, 3L, myMax)
  Scrit <- structure(apply(null.max.perms, 2L, mySort, index), names=conNames)
  for (x in conNames) mtpc.all[Contrast == x, S.crit := Scrit[x]]

  rlefun <- rle_comp(alt)
  rle_obs <- mtpc.all[, rlefun(stat, S.crit), by=key(mtpc.all)]
  crit.regions <- rle_obs[(values)][lengths >= clust.size, list(region=unique(region)), by=Contrast]
  setkeyv(crit.regions, key(mtpc.all))
  mtpc.all[crit.regions,
           A.mtpc := auc_rle(stat, S.crit, thresholds, alt, clust.size),
           by=key(mtpc.all)]

  null.crit <- setNames(vector('list', length(conNames)), conNames)
  Acrit <- setNames(rep.int(0, length(conNames)), conNames)
  for (x in conNames) {
    null.crit[[x]] <- which(apply(null.max.all[, , x], 2L, function(y)
                                  any(with(rlefun(y, Scrit[x]), values == TRUE & lengths >= clust.size))))
    if (length(null.crit) > 0L) {
      null.crits <- abind::adrop(null.max.all[, null.crit[[x]], x, drop=FALSE], drop=3L)
      Acrit[x] <- mean(apply(null.crits, 2L, auc_rle, Scrit[x], thresholds, alt, clust.size))
      mtpc.all[Contrast == x, A.crit := Acrit[x]]
    }
  }

  maxobsfun <- switch(alt,
                      two.sided=function(x) max(is.finite(abs(x)) * abs(x), na.rm=TRUE),
                      less=function(x) min(is.finite(x) * x, na.rm=TRUE),
                      greater=function(x) max(is.finite(x) * x, na.rm=TRUE))
  mtpc.stats <- mtpc.all[, list(S.mtpc=unique(maxobsfun(stat))), by=Contrast]
  mtpc.stats[, tau.mtpc := mtpc.all[mtpc.stats, .SD[S.mtpc == stat, threshold]]]
  mtpc.stats[, S.crit := Scrit]
  mtpc.stats[, A.crit := Acrit]
  setcolorder(mtpc.stats, c('Contrast', 'tau.mtpc', 'S.mtpc', 'S.crit', 'A.crit'))

  glm.attr <- res.glm[[1L]]
  #TODO: may have to remove also 'var.covar'
  #TODO: is it possible for 'removed.subs' and 'df.residual' to be different?
  glm.attr[c('y', 'DT.Xy', 'DT', 'permute', 'runX', 'runY', 'coefficients', 'residuals',
             'sigma', 'fitted.values', 'se', 'perm', 'perm.order')] <- NULL
  if (glm.attr$outcome != glm.attr$measure) glm.attr[c('X', 'qr', 'cov.unscaled')] <- NULL
  for (i in seq_along(thresholds)) {
    res.glm[[i]]$perm[c('null.dist', 'null.max')] <- NULL
    res.glm[[i]]$perm.order <- NULL
  }

  if (isTRUE(long)) null.out <- null.max.all
  out <- c(glm.attr, list(res.glm=res.glm, clust.size=clust.size, DT=mtpc.all, stats=mtpc.stats,
                          null.dist=null.out, perm.order=perms))
  class(out) <- c('mtpc', class(out))
  return(out)
}

#' Run length encoding of the comparison of 2 vectors
#'
#' Returns a function that will compare 2 vectors and return the run length
#' encoding (RLE) of the logical vector from that comparison.
#' @noRd

rle_comp <- function(alt) {
  switch(alt,
         two.sided=function(x, y) rle(abs(x) > abs(y)),
         less=function(x, y) rle(x < y),
         greater=function(x, y) rle(x > y))
}

#' Calculate the AUC for the observed T- or F-statistic
#'
#' @noRd

auc_rle <- function(t.stat, S.crit, thresholds, alt, clust.size) {
  inds <- get_rle_inds(clust.size, alt, t.stat, S.crit, thresholds)
  x <- Map(function(a, b) thresholds[a:b], inds$starts, inds$ends)
  y <- Map(function(a, b) t.stat[a:b], inds$starts, inds$ends)
  auc <- mapply(auc_diff, x, y)
  return(sum(abs(auc) * is.finite(auc), na.rm=TRUE))
}

#' Get the starting and ending points for significant runs
#'
#' Given a vector of observed statistics and a critical value, calculate the
#' \dQuote{runs} (using \code{\link{rle}}) based on comparing these values. Then
#' return the starting and ending points of the runs.
#' @noRd

get_rle_inds <- function(clust.size, alt, t.stat, S.crit, thresholds) {
  compfun <- rle_comp(alt)
  runs <- compfun(t.stat, S.crit)
  sigruns <- which(runs$values == TRUE & runs$lengths >= clust.size)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
#  ends <- runs.len.cumsum[myruns]
#  newindex <- ifelse(sigruns > 1L, sigruns - 1L, 0L)
#  if (0L %in% newindex) starts <- c(1L, starts)
  return(list(starts=starts[sigruns], ends=ends[sigruns]))
}

################################################################################
# S3 METHODS
################################################################################

#' Print a summary of MTPC results
#'
#' @param object,x A \code{mtpc} object
#' @inheritParams summary.bg_GLM
#' @export
#' @rdname mtpc

summary.mtpc <- function(object, contrast=NULL, digits=max(3L, getOption('digits') - 2L), print.head=TRUE, ...) {
  Contrast <- A.mtpc <- A.crit <- stat <- region <- S.mtpc <- NULL
  object$printCon <- contrast
  object$digits <- digits
  object$print.head <- print.head

  # Summary table
  whichmaxfun <- switch(object$alt, two.sided=function(x) which.max(abs(x)), less=which.min, greater=which.max)
  DT.sum <- object$DT[A.mtpc > A.crit, .SD[whichmaxfun(is.finite(stat) * stat)], by=list(Contrast, region)]
  DT.sum <- DT.sum[, c('region', 'stat', 'Contrast', 'threshold', 'S.crit', 'A.mtpc', 'A.crit'), with=FALSE]
  setcolorder(DT.sum, c('Contrast', 'region', 'threshold', 'stat', 'S.crit', 'A.mtpc', 'A.crit'))
  setnames(DT.sum, c('region', 'stat', 'threshold'), c('Region', 'S.mtpc', 'tau.mtpc'))
  setorder(DT.sum, Contrast, -S.mtpc)
  object$DT.sum <- DT.sum

  class(object) <- c('summary.mtpc', class(object))
  return(object)
}

#' @aliases summary.mtpc
#' @method print summary.mtpc
#' @export

print.summary.mtpc <- function(x, ...) {
  print_title_summary('MTPC results (', x$level, '-level)')

  print_model_summary(x)
  print_contrast_type_summary(x)
  print_subs_summary(x)

  cat('# of thresholds: ', length(x$res.glm), '\n')
  cat('Cluster size (across thresholds): ', x$clust.size, '\n')

  print_permutation_summary(x)
  message('\n', 'Statistics', '\n', rep.int('-', getOption('width') / 4))
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
#'   in the caption of the plot. Default: \code{FALSE}
#' @inheritParams plot.bg_GLM
#' @export
#' @rdname mtpc
#'
#' @return The \code{plot} method returns a \code{trellis} object or a list of
#'   \code{ggplot} objects
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
  Contrast <- stat_ribbon <- stat <- threshold <- S.crit <- A.mtpc <- A.crit <- panel.num <-
    xstart <- starts <- xend <- NULL

  stopifnot(inherits(x, 'mtpc'))
  conNames <- x$con.name
  mycontrast <- conNames[contrast]
  DT <- x$DT[Contrast == mycontrast, !'p.fdr']
  DT[, stat_ribbon := stat]  # For filling in supra-threshold areas
  thresholds <- DT[, unique(threshold)]
  DT$nullthresh <- x$stats[Contrast == mycontrast, unique(S.crit)]
  whichmaxfun <- switch(x$alt, two.sided=function(y) which.max(abs(y)), less=which.min, greater=which.max)
  myMax <- switch(x$alt, two.sided=colMaxAbs, less=colMin, greater=colMax)
  thr <- apply(x$null.dist[, , mycontrast], 2L, whichmaxfun)
  thr.y <- myMax(x$null.dist[, , mycontrast])
  nullcoords <- data.table(threshold=thresholds[thr], y=thr.y)

  # 'ggplot2'; Local function to plot for a single region
  plot_single <- function(x, DT, nullcoords, show.null) {
    stat <- S.crit <- threshold <- stat_ribbon <- nullthresh <- y <- A.mtpc <- A.crit <- NULL

    myMax <- switch(x$alt, two.sided=function(y) max(abs(y)), less=min, greater=max)
    inds <- DT[, get_rle_inds(x$clust.size, x$alt, stat, S.crit, threshold)]
    n <- dim(inds)[1L]
    if (n == 0L) {
      DT[, stat_ribbon := NA]
    } else {
      if (n == 1L) {
        all.inds <- c(apply(inds, 1L, function(z) z[1L]:z[2L]))
      } else if (n > 1L) {
        all.inds <- c(unlist(apply(inds, 1L, function(z) z[1L]:z[2L])))
      } else {
        all.inds <- 1L:dim(DT)[1L]
      }
      DT[-all.inds, stat_ribbon := NA]
    }

    lineplot <- ggplot2::ggplot(data=DT, mapping=ggplot2::aes(x=threshold)) +
      ggplot2::geom_line(ggplot2::aes(y=stat), col='red4', size=1.25, na.rm=TRUE) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=nullthresh), lty=2)
    if (n > 0L) {
      if (x$alt == 'less') {
        lineplot <- lineplot +
          ggplot2::geom_ribbon(ggplot2::aes(ymax=stat_ribbon, ymin=nullthresh), fill='red4', alpha=0.3)
      } else {
        lineplot <- lineplot +
          ggplot2::geom_ribbon(ggplot2::aes(ymin=stat_ribbon, ymax=nullthresh), fill='red4', alpha=0.3)
      }
      lineplot <- lineplot + ggplot2::geom_point(ggplot2::aes(y=stat_ribbon), col='red4', size=2, na.rm=TRUE)
    }

    if (isTRUE(all(diff(thresholds) < 0))) lineplot <- lineplot + ggplot2::scale_x_reverse()
    p.region <- if (x$level == 'graph') 'graph-level' else DT[, as.character(unique(region))]
    p.title <- paste0('Region: ', p.region)
    p.subtitle <- paste0('Outcome: ', x$outcome)

    if (isTRUE(show.null)) {
      lineplot <- lineplot + ggplot2::geom_point(data=nullcoords, ggplot2::aes(y=y), col='darkgreen', alpha=0.4, na.rm=TRUE)
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
        ggplot2::labs(caption=statslabel) +
        ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0))
    }
    lineplot <- lineplot +
      ggplot2::labs(title=p.title, subtitle=p.subtitle, x='Threshold', y=paste0(toupper(x$con.type), '-statistic')) +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, face='bold'),
            plot.subtitle=ggplot2::element_text(hjust=0.5))

    return(lineplot)
  }

  # 'base' plotting
  #-------------------------------------------------------
  if (!requireNamespace('ggplot2', quietly=TRUE)) {
    if (x$level == 'vertex') {
      if (is.null(region)) {
        if (isFALSE(only.sig.regions)) {
          myregions <- DT[, levels(region)]
        } else {
          myregions <- droplevels(DT[A.mtpc > A.crit])[, levels(region)]
        }
      } else {
        myregions <- region
      }
    }

    DT <- droplevels(DT[region %in% myregions])
    DT[, panel.num := as.numeric(region)]
    inds <- droplevels(DT[, get_rle_inds(x$clust.size, x$alt, stat, S.crit, threshold), by=region])
    inds[, panel.num := as.numeric(region)]
    inds[, xstart := thresholds[starts]]
    inds[, xend := thresholds[ends]]
    if (isTRUE(show.null)) {
      panelfun <- function(x, y, ...) {
        panel.num <- starts <- nullthresh <- NULL
        panel.points(-nullcoords$threshold, nullcoords$y, pch=19, fill='darkgreen')
        for (i in seq_len(inds[panel.num == panel.number(), .N])) {
          xcoords <- -thresholds[inds[panel.num == panel.number()][i, seq(starts, ends)]]
          ycoords <- DT[panel.num == panel.number() & threshold %in% -xcoords, c(stat, nullthresh)]
          panel.polygon(x=c(xcoords, rev(xcoords)), y=ycoords, col='lightgrey', border=NULL)
        }
        panel.abline(h=DT[, nullthresh], lty=2)
        panel.xyplot(x, y, ..., col='red', lwd=2)
      }
    } else {
      panelfun <- function(x, y, ...) {
        panel.num <- starts <- nullthresh <- NULL
        for (i in seq_len(inds[panel.num == panel.number(), .N])) {
          xcoords <- thresholds[inds[panel.num == panel.number()][i, seq(starts, ends)]]
          ycoords <- DT[panel.num == panel.number() & threshold %in% xcoords, c(stat, nullthresh)]
          panel.polygon(x=c(xcoords, rev(xcoords)), y=ycoords, col='lightgrey', border=NULL)
        }
        panel.abline(h=DT[, nullthresh], lty=2)
        panel.xyplot(x, y, ..., col='red', lwd=2)
      }
    }

    lineplot <- xyplot(stat ~ -threshold | region, data=DT, type=c('l', 'p'), pch=19,
                       xlab='Threshold', ylab=paste0(toupper(x$con.type), '-statistic'),
                       main=paste0('Outcome: ', x$outcome),
                       scales=list(y=list(relation='free')),
                       panel=panelfun)
    return(lineplot)

  # 'ggplot2' plotting
  #-------------------------------------------------------
  } else {
    if (x$level == 'graph') {
      lineplot <- plot_single(x, DT, nullcoords, show.null)
      return(lineplot)
    } else if (x$level == 'vertex') {
      if (is.null(region)) {
        if (isFALSE(only.sig.regions)) {
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
}

#' @export
#' @rdname mtpc
nobs.mtpc <- function(object, ...) nobs(object$res.glm[[1L]])

#' @export
#' @rdname mtpc
terms.mtpc <- function(x, ...) terms(x$res.glm[[1]])

#' @export
#' @rdname mtpc
formula.mtpc <- formula.bg_GLM

#' @export
#' @rdname mtpc
labels.mtpc <- labels.bg_GLM

#' @method case.names mtpc
#' @export
#' @rdname mtpc
case.names.mtpc <- function(object, ...) case.names(object$res.glm[[1]])

#' @method variable.names mtpc
#' @export
#' @rdname mtpc
variable.names.mtpc <- function(object, ...) variable.names(object$res.glm[[1]])

#' @method df.residual mtpc
#' @export
#' @rdname mtpc
df.residual.mtpc <- df.residual.bg_GLM
