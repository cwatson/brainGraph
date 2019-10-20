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
#' @importFrom ade4 mantel.rtest
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
  Group <- Study.ID <- i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'), inherits(corrs, 'corr_mats'))
  level <- match.arg(level)
  group.vec <- groups(resids)
  n <- nobs(resids)
  if (level == 'global') {
    combFun <- c
    diffFun <- function(a, b) 1 - mantel.rtest(as.dist(a), as.dist(b), nrepet=1e3)$obs
  } else if (level == 'regional') {
    combFun <- rbind
    diffFun <- function(a, b) colSums(abs(a - b))
  }
  IC <- foreach(i=seq_len(n), .combine=combFun) %dopar% {
    resids.excl <- resids[-i]
    new.corrs <- corr.matrix(resids.excl[group.vec[i]], densities=0.1)
    diffFun(corrs$R[, , group.vec[i]], new.corrs$R[, , 1L])
  }

  DT <- cbind(resids$resids.all[, list(Study.ID, Group)], IC)
  if (level == 'regional') {
    DT <- melt(DT, id.vars=c('Study.ID', 'Group'), variable.name='region', value.name='RC')
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
#'   (default: 1)
#' @export
#' @importFrom ade4 mantel.rtest
#'
#' @rdname individ_contrib
#' @examples
#' \dontrun{
#' IC <- aop(resids.all, corrs)
#' RC <- aop(resids.all, corrs, level='regional')
#' }

aop <- function(resids, corrs, level=c('global', 'regional'), control.value=1L) {
  Group <- i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'), inherits(corrs, 'corr_mats'))

  corr.mat <- corrs[, control.value]$R[, , 1L]
  grps <- groups(resids)
  kNumSubj <- table(grps)
  grps <- unique(grps)
  if (is.numeric(control.value)) control.value <- grps[control.value]
  patient.str <- setdiff(grps, control.value)

  control.inds <- resids$resids.all[Group == control.value, which=TRUE]
  level <- match.arg(level)
  if (level == 'global') {
    combFun <- c
    diffFun <- function(a, b) 1 - mantel.rtest(as.dist(a), as.dist(b), nrepet=1e3)$obs
  } else if (level == 'regional') {
    combFun <- rbind
    diffFun <- function(a, b) data.table(t(colSums(abs(a - b))))
  }
  IC <- setNames(vector('list', length(patient.str)), patient.str)
  for (j in patient.str) {
    pat.inds <- resids$resids.all[Group == j, which=TRUE]
    IC[[j]] <- foreach(i=seq_len(kNumSubj[j]), .combine=combFun) %dopar% {
      resids.aop <- resids[c(control.inds, pat.inds[i])]
      resids.aop$resids.all[, Group := control.value]
      resids.aop$resids.all <- droplevels(resids.aop$resids.all)
      resids.aop$Group <- resids.aop$resids.all[, factor(levels(Group))]
      setkey(resids.aop$resids.all, Group)
      new.corr <- corr.matrix(resids.aop, densities=0.1)$R[, , 1L]
      diffFun(corr.mat, new.corr)
    }
    IC[[j]] <- cbind(resids$resids.all[j, c('Study.ID', 'Group')], IC[[j]])
  }
  DT <- rbindlist(IC)
  if (level == 'global') {
    setnames(DT, 'V2', 'IC')
  } else if (level == 'regional') {
    DT <- melt(DT, id.vars=c('Study.ID', 'Group'), variable.name='region', value.name='RC')
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
  avg <- diff_mean <- Group <- IC <- Max <- med <- Min <- RC <- se <- stdev <- NULL
  object$digits <- digits
  DT.sum <- copy(object$DT)
  if (object$level == 'regional') {
    regions <- if (is.null(region)) region.names(DT.sum) else region
    object$regions <- regions
    DT.sum <- droplevels(DT.sum[region %in% regions])
    DT.sum[, Min := min(RC), by=list(Group, region)]
    DT.sum[, med := median(RC), by=list(Group, region)]
    DT.sum[, avg := mean(RC), by=list(Group, region)]
    DT.sum[, Max := max(RC), by=list(Group, region)]
    DT.sum[, stdev := sd(RC), by=list(Group, region)]
    DT.sum[, se := stdev / sqrt(.N), by=list(Group, region)]
    outliers <- DT.sum[, .SD[RC > avg + 2 * stdev], by=list(Group, region)]
    outliers.reg <- outliers[, .N, by=region]
    outliers.reg.vec <- with(outliers.reg, structure(N, names=as.character(region)))
    object$outliers <- list(DT=outliers, region=outliers.reg.vec)
  } else if (object$level == 'global') {
    DT.sum[, avg := mean(IC), by=Group]
    outliers <- DT.sum[, .SD[IC > avg + 2 * sd(IC)], by=Group]
    outliers[, diff_mean := IC - avg]
    outliers[, avg := NULL]
    DT.sum[, avg := NULL]
    object$outliers$DT <- outliers
  }
  object$DT.sum <- DT.sum

  # Calculate some group descriptive statistics
  if (object$level == 'global') {
    grps <- object$DT.sum[, levels(Group)]
    sums <- setNames(vector('list', length(grps)), grps)
    for (g in grps) {
      sums[[g]] <- object$DT.sum[Group == g, quantile(IC, c(0, .1, .25, .5, .75, .9, 1))]
      sums[[g]] <- append(sums[[g]], object$DT.sum[Group == g, mean(IC)], after=4)
      names(sums[[g]]) <- c('Min.', '10%', '1st Qu.', 'Median', 'Mean', '3rd Qu.', '90%', 'Max.')
    }
    sums <- t(abind::abind(sums, along=2))
    object$sums <- sums
  }
  class(object) <- c('summary.IC', class(object))
  return(object)
}

#' @aliases summary.IC
#' @method print summary.IC
#' @keywords internal

print.summary.IC <- function(x, ...) {
  Group <- region <- IC <- NULL
  print_title_summary('Individual contributions')
  cat('Method: ', x$method, '\n')
  cat('Level: ', x$level, '\n\n')

  if (x$level == 'regional') {
    print(x$DT.sum[region == levels(region)[1], table(Group)])
  } else {
    print(x$DT.sum[, table(Group)])
  }
  cat('\n')

  width <- getOption('width')
  dashes <- rep('-', width / 4)
  if (x$level == 'global') {
    message('Group summaries\n', dashes)
    print(x$sums, digits=x$digits)
    cat('\n')
  }

  if (x$level == 'regional') {
    message('# of outliers per region: (sorted in descending order)\n', dashes)
    print(sort(x$outliers$region, decreasing=TRUE))
    cat('\n')
    DT <- x$DT.sum[, .SD[1, !c('Study.ID', 'RC')], by=list(Group, region)]
    setnames(DT, c('med', 'avg', 'stdev', 'se'), c('Median', 'Mean', 'Std. Dev', 'Std. Err'))
    message('Region summaries\n', dashes)
    print(DT, digits=x$digits)
  } else {
    message('Outliers\n', dashes)
    print(x$outliers$DT[order(Group, -IC)], digits=x$digits)
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
#' @importFrom ggrepel geom_text_repel
#' @rdname individ_contrib

plot.IC <- function(x, plot.type=c('mean', 'smooth', 'boxplot'), region=NULL, ids=TRUE, ...) {
  RC <- Group <- avg <- se <- ind <- mark <- IC <- Study.ID <- NULL
  DT <- summary(x)$DT.sum
  kNumGroups <- DT[, nlevels(Group)]
  grps <- DT[, levels(Group)]
  leg.pos <- if (kNumGroups == 1) 'none' else 'bottom'
  if (x$level == 'regional') {
    xlabel <- 'Region'
    ylabel <- 'Regional contribution'
  } else {
    xlabel <- 'Subject'
    ylabel <- 'Individual contribution'
  }
  ptitle <- paste0(ylabel, 's, ', tolower(x$method), ' method')

  if (x$level == 'regional') {
    regions <- if (is.null(region)) region.names(DT) else region
    txtsize <- if (length(regions) > 50) 6 else 9

    plot.type <- match.arg(plot.type)
    if (plot.type == 'boxplot') {
      p <- ggplot(DT[region %in% regions], aes(x=region, y=RC)) +
        geom_boxplot(aes(fill=Group, group=interaction(Group, region)))

    } else {
      p <- ggplot(DT[region %in% regions], aes(x=region, col=Group, group=Group))
      if (plot.type == 'smooth') {
        p <- p + stat_smooth(method='loess', aes(y=RC))

      } else if (plot.type == 'mean') {
        DT[, avg := mean(RC), by=list(Group, region)]
        DT[, se := sd(RC) / sqrt(.N), by=list(Group, region)]
        p <- p + geom_line(aes(y=avg)) +
          geom_ribbon(aes(ymin=avg-se, ymax=avg+se, fill=Group), alpha=0.5)
      }
    }

  } else {
    n <- dim(DT)[1L]
    txtsize <- if (n > 50) 9 else 12
    if (isTRUE(ids)) {
      DT[, ind := Study.ID]
    } else {
      spec <- paste0('%0', floor(log10(n) + 1), 'i')
      DT[, ind := sprintf(spec, .I)]
    }
    DT[, mark := ifelse(IC > mean(IC) + 2*sd(IC), 1, 0)]
    DT[mark == 0, ind := '']
    DT[, mark := as.factor(mark)]
    # From "ggsci"; the "npg" palette
    cols <- c('#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
              '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85')
    p <- ggplot(DT, aes(x=Study.ID, y=IC, col=Group)) +
      geom_text_repel(aes(label=ind), size=3) +
      geom_point(aes(shape=mark, size=mark)) +
      scale_color_manual(name='Group', labels=grps, values=cols[1:kNumGroups]) +
      scale_shape_manual(name='Group', labels=grps, values=c(20, 17)) +
      scale_size_manual(name='Group', labels=grps, values=c(2, 3))
  }
  p <- p +
    labs(x=xlabel, y=ylabel, title=ptitle) +
    theme(legend.position=leg.pos,
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=txtsize, angle=45, vjust=0.5),
          plot.title=element_text(hjust=0.5, face='bold'))
  return(p)
}
