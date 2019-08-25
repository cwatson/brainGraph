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
#' @aliases loo
#' @rdname individ_contrib
#' @examples
#' \dontrun{
#' IC <- loo(resids.all, corrs)
#' RC <- loo(resids.all, corrs, level='regional')
#' }
#' @family Group analysis functions
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
  group.vec <- resids$resids.all$Group
  group.num <- as.integer(group.vec)
  group.vec <- as.character(group.vec)
  if (level == 'global') {
    IC <- foreach(i=seq_len(nrow(resids$resids.all)), .combine='c') %dopar% {
      resids.excl <- resids[-i]
      new.corrs <- corr.matrix(resids.excl[group.vec[i]], densities=0.1)

      1 - mantel.rtest(as.dist(corrs$R[, , group.num[i]]),
                       as.dist(new.corrs$R[, , 1]),
                       nrepet=1e3)$obs
    }

    DT <- data.table(resids$resids.all[, list(Study.ID, Group)], IC=IC)
  } else if (level == 'regional') {
    RC <- foreach(i=seq_len(nrow(resids$resids.all)), .combine='rbind') %dopar% {
      resids.excl <- resids[-i]
      new.corrs <- corr.matrix(resids.excl[group.vec[i]], densities=0.1)
      colSums(abs(corrs$R[, , group.num[i]] - new.corrs$R[, , 1]))
    }
    RC.dt <- cbind(resids$resids.all[, list(Study.ID, Group)], RC)
    DT <- melt(RC.dt, id.vars=c('Study.ID', 'Group'),
                 variable.name='region', value.name='RC')
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
#' @aliases aop
#' @rdname individ_contrib
#' @examples
#' \dontrun{
#' IC <- aop(resids.all, corrs)
#' RC <- aop(resids.all, corrs, level='regional')
#' }

aop <- function(resids, corrs, level=c('global', 'regional'), control.value=1L) {
  Group <- Study.ID <- i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'), inherits(corrs, 'corr_mats'))

  corr.mat <- corrs[, control.value]$R[, , 1]
  groups <- resids$groups
  kNumSubj <- resids$resids.all[, tabulate(Group)]
  if (is.numeric(control.value)) {
    control.int <- control.value
    control.str <- groups[control.int]
  } else {
    control.int <- which(groups %in% control.value)
    control.str <- control.value
  }
  patient.str <- groups[-control.int]
  patient.int <- seq_along(groups)[-control.int]

  control.inds <- resids$resids.all[, which(Group == control.str)]
  level <- match.arg(level)
  if (level == 'global') {
    IC <- sapply(groups[-control.int], function(x) NULL)
    for (j in patient.int) {
      pat.inds <- resids$resids.all[, which(Group == patient.str)]
      IC[[groups[j]]] <- foreach(i=seq_len(kNumSubj[j]), .combine='c') %dopar% {
        resids.aop <- resids[c(control.inds, pat.inds[i])]
        resids.aop$resids.all[, Group := control.str]
        resids.aop$resids.all <- droplevels(resids.aop$resids.all)
        resids.aop$groups <- resids.aop$resids.all[, factor(levels(Group))]
        setkey(resids.aop$resids.all, Group)
        new.corr <- corr.matrix(resids.aop, densities=0.1)$R[, , 1]
        1 - mantel.rtest(as.dist(corr.mat),
                         as.dist(new.corr),
                         nrepet=1e3)$obs
      }
      IC[[groups[j]]] <- cbind(resids$resids.all[groups[j], c('Study.ID', 'Group')], IC[[groups[j]]])
    }
    DT <- rbindlist(IC)
    setnames(DT, 'V2', 'IC')

  } else if (level == 'regional') {
    RC <- sapply(groups[-control.int], function(x) NULL)
    for (j in patient.int) {
      pat.inds <- resids$resids.all[, which(Group == patient.str)]
      RC[[groups[j]]] <- foreach(i=seq_len(kNumSubj[j]), .combine='rbind') %dopar% {
        resids.aop <- resids[c(control.inds, pat.inds[i])]
        resids.aop$resids.all[, Group := control.str]
        resids.aop$resids.all <- droplevels(resids.aop$resids.all)
        resids.aop$groups <- resids.aop$resids.all[, factor(levels(Group))]
        setkey(resids.aop$resids.all, Group)
        new.corr <- corr.matrix(resids.aop, densities=0.1)$R[, , 1]
        data.table(t(colSums(abs(corr.mat - new.corr))))
      }
      RC[[groups[j]]] <- cbind(resids$resids.all[groups[j], c('Study.ID', 'Group')], RC[[groups[j]]])
    }
    RC.dt <- rbindlist(RC)
    DT <- melt(RC.dt, id.vars=c('Study.ID', 'Group'),
                variable.name='region', value.name='RC')
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
#' @export
#' @method summary IC
#' @rdname individ_contrib

summary.IC <- function(object, region=NULL, digits=max(3L, getOption('digits') - 2L), ...) {
  avg <- RC <- Group <- stdev <- se <- IC <- NULL
  object$digits <- digits
  DT.sum <- copy(object$DT)
  if (object$level == 'regional') {
    regions <- if (is.null(region)) DT.sum[, levels(region)] else region
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
  } else {
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
    groups <- object$DT.sum[, levels(Group)]
    sums <- setNames(vector('list', length(groups)), groups)
    for (g in groups) {
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
#' @method plot IC
#' @rdname individ_contrib

plot.IC <- function(x, plot.type=c('mean', 'smooth', 'boxplot'), region=NULL, ids=TRUE, ...) {
  RC <- Group <- avg <- se <- ind <- mark <- IC <- Study.ID <- NULL
  DT <- summary(x)$DT.sum
  kNumGroups <- DT[, nlevels(Group)]
  groups <- DT[, levels(Group)]
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
    regions <- if (is.null(region)) DT[, levels(region)] else region
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
    n <- nrow(DT)
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
    cols <- c('#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F3B7F',
              '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85')
    p <- ggplot(DT, aes(x=Study.ID, y=IC, col=Group)) +
      geom_text_repel(aes(label=ind), size=3) +
      geom_point(aes(shape=mark, size=mark)) +
      scale_color_manual(name='Group', labels=groups, values=cols[1:kNumGroups]) +
      scale_shape_manual(name='Group', labels=groups, values=c(20, 17)) +
      scale_size_manual(name='Group', labels=groups, values=c(2, 3))
  }
  p <- p +
    labs(x=xlabel, y=ylabel, title=ptitle) +
    theme(legend.position=leg.pos,
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=txtsize, angle=45, vjust=0.5),
          plot.title=element_text(hjust=0.5, face='bold'))
  return(p)
}
