#' Approaches to estimate individual network contribution
#'
#' \code{loo} calculates the individual contribution to group network data for
#' each subject in each group using a \dQuote{leave-one-out} approach. The
#' residuals of a single subject are excluded, and a correlation matrix is
#' created. This is compared to the original correlation matrix using the Mantel
#' test.
#'
#' @param resids An object of class \code{brainGraph_resids} (the output from
#'   \code{\link{get.resid}})
#' @param corrs List of lists of correlation matrices (as output by
#'   \code{\link{corr.matrix}}).
#' @param level Character string; the level at which you want to calculate
#'   contributions (either \code{global} or \code{regional})
#' @export
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#'
#' @return A \code{data.table} with columns for
#'   \item{Study.ID}{Subject identifier}
#'   \item{Group}{Group membership}
#'   \item{region}{If \code{level='regional'}}
#'   \item{IC,RC}{The value of the individual/regional contributions}
#'
#' @family Structural covariance network functions
#' @name IndividualContributions
#' @rdname individ_contrib
#' @examples
#' \dontrun{
#' IC <- loo(resids.all, corrs)
#' RC <- loo(resids.all, corrs, level='regional')
#' }
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Saggar, M. and Hosseini, S.M.H. and Buno, J.L. and Quintin, E.
#'   and Raman, M.M. and Kesler, S.R. and Reiss, A.L. (2015) Estimating
#'   individual contributions from group-based structural correlations networks.
#'   \emph{NeuroImage}, \bold{120}, 274--284.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2015.07.006}

loo <- function(resids, corrs, level=c('global', 'regional')) {
  ..sID <- ..gID <- NULL
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'), inherits(corrs, 'corr_mats'))
  level <- match.arg(level)
  group.vec <- groups(resids)
  n <- nobs(resids)
  if (level == 'global') {
    if (!requireNamespace('ade4', quietly=TRUE)) stop('Must install the "ade4" package.')
    combFun <- c
    diffFun <- function(a, b) 1 - ade4::mantel.rtest(as.dist(a), as.dist(b), nrepet=1e3)$obs
  } else if (level == 'regional') {
    combFun <- rbind
    diffFun <- function(a, b) colSums(abs(a - b))
  }
  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  IC <- foreach(i=seq_len(n), .combine=combFun) %dopar% {
    resids.excl <- resids[-i]
    new.corrs <- corr.matrix(resids.excl[group.vec[i]], densities=0.1)
    diffFun(corrs$R[, , group.vec[i]], new.corrs$R[, , 1L])
  }

  DT <- cbind(resids$resids.all[, c(..sID, ..gID)], IC)
  if (level == 'regional') {
    DT <- melt(DT, id.vars=c(sID, gID), variable.name='region', value.name='RC')
  }
  out <- list(method='Leave one out', level=level, DT=DT)
  class(out) <- c('IC', class(out))
  return(out)
}

#' Add-one-patient approach to estimate individual network contribution
#'
#' \code{aop} calculates the individual contribution using an
#' \dQuote{add-one-patient} approach. The residuals of a single patient are
#' added to those of a control group, and a correlation matrix is created. This
#' is repeated for all individual patients and each patient group.
#'
#' @note For \code{aop}, it is assumed by default that the control group is the
#'   first group.
#'
#' @param control.value Integer or character string specifying the control group
#'   (default: 1L)
#' @export
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#'
#' @rdname individ_contrib
#' @examples
#' \dontrun{
#' IC <- aop(resids.all, corrs)
#' RC <- aop(resids.all, corrs, level='regional')
#' }

aop <- function(resids, corrs, level=c('global', 'regional'), control.value=1L) {
  ..sID <- ..gID <- NULL
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'), inherits(corrs, 'corr_mats'))

  corr.mat <- corrs[, control.value]$R[, , 1L]
  grps <- groups(resids)
  kNumSubj <- table(grps)
  grps <- unique(grps)
  if (is.numeric(control.value)) control.value <- grps[control.value]
  patient.str <- setdiff(grps, control.value)

  control.inds <- resids$resids.all[get(gID) == control.value, which=TRUE]
  level <- match.arg(level)
  if (level == 'global') {
    if (!requireNamespace('ade4', quietly=TRUE)) stop('Must install the "ade4" package.')
    combFun <- c
    diffFun <- function(a, b) 1 - ade4::mantel.rtest(as.dist(a), as.dist(b), nrepet=1e3)$obs
  } else if (level == 'regional') {
    combFun <- rbind
    diffFun <- function(a, b) data.table(t(colSums(abs(a - b))))
  }

  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  IC <- setNames(vector('list', length(patient.str)), patient.str)
  for (j in patient.str) {
    pat.inds <- resids$resids.all[get(gID) == j, which=TRUE]
    IC[[j]] <- foreach(i=seq_len(kNumSubj[j]), .combine=combFun) %dopar% {
      resids.aop <- resids[c(control.inds, pat.inds[i])]
      resids.aop$resids.all[, eval(gID) := control.value]
      resids.aop$resids.all <- droplevels(resids.aop$resids.all)
      resids.aop$Group <- resids.aop$resids.all[, levels(get(gID))]
      setkeyv(resids.aop$resids.all, gID)
      new.corr <- corr.matrix(resids.aop, densities=0.1)$R[, , 1L]
      diffFun(corr.mat, new.corr)
    }
    IC[[j]] <- cbind(resids$resids.all[j, c(..sID, ..gID)], IC[[j]])
  }
  DT <- rbindlist(IC)
  if (level == 'global') {
    setnames(DT, 'V2', 'IC')
  } else if (level == 'regional') {
    DT <- melt(DT, id.vars=c(sID, gID), variable.name='region', value.name='RC')
  }
  out <- list(method='Add one patient', level=level, DT=DT)
  class(out) <- c('IC', class(out))
  return(out)
}

#' Print a summary of individual contribution estimates
#'
#' The \code{summary} method prints the group/region-wise means and standard
#' deviations.
#'
#' @param object,x A \code{IC} object
#' @param region Character vector specifying which regions' IC's to print. Only
#'   relevant if \code{method='Leave one out'}
#' @param ... Unused
#' @inheritParams summary.bg_GLM
#' @export
#' @rdname individ_contrib

summary.IC <- function(object, region=NULL, digits=max(3L, getOption('digits') - 2L), ...) {
  avg <- diff_from_mean <- IC <- Max <- med <- Min <- RC <- se <- stdev <- NULL
  gID <- getOption('bg.group')
  object$digits <- digits
  DT.sum <- copy(object$DT)
  if (object$level == 'regional') {
    regions <- if (is.null(region)) region.names(DT.sum) else region
    object$regions <- regions
    DT.sum <- droplevels(DT.sum[region %in% regions])
    DT.sum[, Min := min(RC), by=list(get(gID), region)]
    DT.sum[, med := median(RC), by=list(get(gID), region)]
    DT.sum[, avg := mean(RC), by=list(get(gID), region)]
    DT.sum[, Max := max(RC), by=list(get(gID), region)]
    DT.sum[, stdev := sd(RC), by=list(get(gID), region)]
    DT.sum[, se := stdev / sqrt(.N), by=list(get(gID), region)]
    outliers <- DT.sum[, .SD[RC > avg + 2 * stdev], by=list(get(gID), region)]
    setnames(outliers, 'get', gID)
    outliers.reg <- outliers[, .N, by=region]
    outliers.reg.vec <- with(outliers.reg, structure(N, names=as.character(region)))
    object$outliers <- list(DT=outliers, region=outliers.reg.vec)
  } else if (object$level == 'global') {
    DT.sum[, avg := mean(IC), by=gID]
    outliers <- DT.sum[, .SD[IC > avg + 2 * sd(IC)], by=gID]
    outliers[, diff_from_mean := IC - avg]
    outliers[, avg := NULL]
    DT.sum[, avg := NULL]
    object$outliers$DT <- outliers
  }
  object$DT.sum <- DT.sum

  # Calculate some group descriptive statistics
  if (object$level == 'global') {
    grps <- object$DT.sum[, levels(factor(get(gID)))]
    sums <- setNames(vector('list', length(grps)), grps)
    for (g in grps) {
      sums[[g]] <- object$DT.sum[get(gID) == g, quantile(IC, c(0, .1, .25, .5, .75, .9, 1))]
      sums[[g]] <- append(sums[[g]], object$DT.sum[get(gID) == g, mean(IC)], after=4L)
      names(sums[[g]]) <- c('Min.', '10%', '1st Qu.', 'Median', 'Mean', '3rd Qu.', '90%', 'Max.')
    }
    sums <- t(abind::abind(sums, along=2L))
    object$sums <- sums
  }
  class(object) <- c('summary.IC', class(object))
  return(object)
}

#' @aliases summary.IC
#' @method print summary.IC
#' @export

print.summary.IC <- function(x, ...) {
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  region <- IC <- NULL
  print_title_summary('Individual contributions')
  cat('Method: ', x$method, '\n')
  cat('Level: ', x$level, '\n\n')

  if (x$level == 'regional') {
    print(x$DT.sum[region == levels(region)[1L], table(get(gID))])
  } else {
    print(x$DT.sum[, table(get(gID))])
  }
  cat('\n')

  width <- getOption('width')
  dashes <- rep.int('-', width / 4)
  if (x$level == 'global') {
    message('Group summaries\n', dashes)
    print(x$sums, digits=x$digits)
    cat('\n')
  }

  if (x$level == 'regional') {
    message('# of outliers per region: (sorted in descending order)\n', dashes)
    print(sort(x$outliers$region, decreasing=TRUE))
    cat('\n')
    DT <- x$DT.sum[, .SD[1L, !c(sID, 'RC'), with=FALSE], by=list(get(gID), region)]
    setnames(DT, c('get', 'med', 'avg', 'stdev', 'se'), c(gID, 'Median', 'Mean', 'Std. Dev', 'Std. Err'))
    message('Region summaries\n', dashes)
    print(DT, digits=x$digits)
  } else {
    message('Outliers\n', dashes)
    print(x$outliers$DT[order(get(gID), -IC)], digits=x$digits)
  }
  invisible(x)
}


#' Plot regional contributions estiamtes
#'
#' The \code{plot} method is only valid for \emph{regional} contribution
#' estimates, and plots the average regional contribution for each
#' vertex/region.
#'
#' @param plot.type Character string indicating the type of plot; the default is
#'   to plot the mean (along with standard errors)
#' @param ids Logical indicating whether to plot Study ID's for outliers.
#'   Otherwise plots the integer index
#' @export
#' @rdname individ_contrib

plot.IC <- function(x, plot.type=c('mean', 'smooth', 'boxplot'), region=NULL, ids=TRUE, ...) {
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  RC <- avg <- se <- ind <- mark <- IC <- NULL
  DT <- summary(x)$DT.sum
  kNumGroups <- DT[, nlevels(factor(get(gID)))]
  grps <- DT[, levels(factor(get(gID)))]
  leg.pos <- if (kNumGroups == 1L) 'none' else 'bottom'

  if (x$level == 'regional') {
    xlabel <- 'Region'
    ylabel <- 'Regional contribution'
    regions <- if (is.null(region)) region.names(DT) else region
    txtsize <- if (length(regions) > 50) 6 else 9

    # 'base' plotting
    if (!requireNamespace('ggplot2', quietly=TRUE)) {
      DT2 <- DT[region %in% regions, .SD[1L], by=c(gID, 'region')]
      ymin <- DT2[, min(avg - se)]
      ymax <- DT2[, max(avg + se)]
      if (kNumGroups > 1L) par(mar=c(8.6, 4.1, 2.1, 0.4), xpd=TRUE)
      plot(0, type='n', xlim=c(0, length(regions)), ylim=extendrange(c(ymin, ymax)),
           xlab=xlabel, ylab=ylabel, xaxt='n')
      for (i in seq_len(kNumGroups)) {
        DT2[get(gID) == grps[i], lines(avg, col=plot.cols[i])]
        DT2[get(gID) == grps[i],
            polygon(x=c(region, rev(region)), y=c(avg + se, rev(avg - se)),
                    col=adjustcolor(plot.cols[i], alpha.f=0.4), border=plot.cols[i])]
      }
      text(x=seq_along(regions), par('usr')[3L], labels=DT2$region,
           srt=45, pos=1, offset=0.8, cex=0.5, xpd=TRUE)
      if (kNumGroups > 1L) {
        legend('bottom', title=gID, grps, fill=plot.cols[1L:kNumGroups],
               inset=c(0, -0.35), horiz=TRUE)
      }

    # 'ggplot2' plotting
    } else {
      plot.type <- match.arg(plot.type)
      if (plot.type == 'boxplot') {
        p <- ggplot2::ggplot(DT[region %in% regions], ggplot2::aes(x=region, y=RC)) +
          ggplot2::geom_boxplot(ggplot2::aes(fill=get(gID), group=interaction(get(gID), region)))

      } else {
        p <- ggplot2::ggplot(DT[region %in% regions], ggplot2::aes(x=region, col=get(gID), group=get(gID)))
        if (plot.type == 'smooth') {
          p <- p + ggplot2::stat_smooth(method='loess', ggplot2::aes(y=RC))

        } else if (plot.type == 'mean') {
          DT[, avg := mean(RC), by=list(get(gID), region)]
          DT[, se := sd(RC) / sqrt(.N), by=list(get(gID), region)]
          p <- p + ggplot2::geom_line(ggplot2::aes(y=avg)) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin=avg-se, ymax=avg+se, fill=get(gID)), alpha=0.5)
        }
      }
    }

  } else {
    xlabel <- 'Subject'
    ylabel <- 'Individual contribution'
    n <- dim(DT)[1L]
    txtsize <- if (n > 50) 9 else 12
    if (isTRUE(ids)) {
      DT[, ind := get(sID)]
    } else {
      spec <- paste0('%0', floor(log10(n) + 1), 'i')
      DT[, ind := sprintf(spec, .I)]
    }
    DT[, mark := ifelse(IC > mean(IC) + 2*sd(IC), 1, 0)]
    DT[mark == 0, ind := '']
    DT[, mark := as.factor(mark)]

    # 'base' plotting
    if (!requireNamespace('ggplot2', quietly=TRUE)) {
      if (kNumGroups > 1L) par(mar=c(8.6, 4.1, 2.1, 0.4), xpd=TRUE, xaxt='n')
      plot(0, type='n', xlim=c(0, n), ylim=extendrange(DT$IC),
           xlab=xlabel, ylab=ylabel, xaxt='n')
      DT[, plot(IC, type='p', pch=19, col=plot.cols[get(gID)])]
      text(x=seq_len(n), par('usr')[3L], labels=DT[, get(sID)],
           srt=45, pos=1, offset=0.8, cex=0.5, xpd=TRUE)
      if (kNumGroups > 1L) {
        legend('bottom', title=gID, grps, fill=plot.cols[1L:kNumGroups],
               inset=c(0, -0.35), horiz=TRUE)
      }

    # 'ggplot2' plotting
    } else {
      textfun <- if (!requireNamespace('ggrepel', quietly=TRUE)) ggplot2::geom_text
        else ggrepel::geom_text_repel
      p <- ggplot2::ggplot(DT, ggplot2::aes(x=get(sID), y=IC, col=get(gID))) +
        textfun(ggplot2::aes(label=ind), size=3) +
        ggplot2::geom_point(ggplot2::aes(shape=mark, size=mark)) +
        ggplot2::scale_color_manual(name=gID, labels=grps, values=plot.cols[1L:kNumGroups]) +
        ggplot2::scale_shape_manual(name=gID, labels=grps, values=c(20, 17)) +
        ggplot2::scale_size_manual(name=gID, labels=grps, values=c(2, 3))
    }
  }

  ptitle <- paste0(ylabel, 's, ', tolower(x$method), ' method')
  if (!requireNamespace('ggplot2', quietly=TRUE)) {
    title(ptitle)
  } else {
    p <- p +
      ggplot2::labs(x=xlabel, y=ylabel, title=ptitle, fill=gID, col=gID) +
      ggplot2::theme(legend.position=leg.pos,
                     panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(size=txtsize, angle=45, vjust=0.5),
                     plot.title=ggplot2::element_text(hjust=0.5, face='bold'))
    return(p)
  }
}
