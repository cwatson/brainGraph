#' Approaches to estimate individual network contribution
#'
#' \code{loo} calculates the individual contribution to group network data for
#' each subject in each group using a "leave-one-out" approach. The residuals of
#' a single subject are excluded, and a correlation matrix is created. This is
#' compared to the original correlation matrix using the Mantel test.
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
#' @references Saggar M., Hosseini S.M.H., Buno J.L., Quintin E., Raman M.M.,
#'   Kesler S.R., Reiss A.L. (2015) \emph{Estimating individual contributions
#'   from group-based structural correlations networks}. NeuroImage,
#'   120:274-284. \url{https://dx.doi.org/10.1016/j.neuroimage.2015.07.006}

loo <- function(resids, corrs, level=c('global', 'regional')) {
  Group <- Study.ID <- i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'), inherits(corrs, 'corr_mats'))
  level <- match.arg(level)
  group.vec <- resids$resids.all$Group
  group.num <- as.integer(group.vec)
  group.vec <- as.character(group.vec)
  if (level == 'global') {
    IC <- foreach (i=seq_len(nrow(resids$resids.all)), .combine='c') %dopar% {
      resids.excl <- resids[-i]
      new.corrs <- corr.matrix(resids.excl[group.vec[i]], densities=0.1)

      1 - mantel.rtest(as.dist(corrs$R[, , group.num[i]]),
                       as.dist(new.corrs$R[, , 1]),
                       nrepet=1e3)$obs
    }

    DT <- data.table(resids$resids.all[, list(Study.ID, Group)], IC=IC)
  } else if (level == 'regional') {
    RC <- foreach (i=seq_len(nrow(resids$resids.all)), .combine='rbind') %dopar% {
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

#' "Add-one-patient" approach to estimate individual network contribution
#'
#' \code{aop} calculates the individual contribution using an "add-one-patient"
#' approach. The residuals of a single patient are added to those of a control
#' group, and a correlation matrix is created. This is repeated for all
#' individual patients and each patient group.
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
#' @param object A \code{IC} object
#' @param region Character vector of regions to plot; default is to plot for all
#'   regions
#' @param ... Unused
#' @export
#' @method summary IC
#' @rdname individ_contrib

summary.IC <- function(object, region=NULL, ...) {
  avg <- RC <- Group <- stdev <- se <- IC <- NULL
  DT.sum <- copy(object$DT)
  if (object$level == 'regional') {
    if (is.null(region)) {
      regions <- DT.sum[, levels(region)]
    } else {
      regions <- region
    }
    DT.sum[, avg := mean(RC), by=list(Group, region)]
    DT.sum[, stdev := sd(RC), by=list(Group, region)]
    DT.sum[, se := stdev / sqrt(.N), by=list(Group, region)]
    outliers <- DT.sum[, .SD[RC > mean(RC) + 2 * stdev], by=list(Group, region)]
    outliers.reg <- outliers[, .N, by=region]
    outliers.reg.vec <- structure(outliers.reg$N, names=as.character(outliers.reg$region))
    object$outliers <- list(DT=outliers, region=outliers.reg.vec)
  } else {
    object$outliers$DT <- DT.sum[, .SD[IC > mean(IC) + 2 * sd(IC)], by=Group]
  }
  object$DT.sum <- DT.sum
  class(object) <- c('summary.IC', class(object))
  return(object)
}

#' @aliases summary.IC
#' @method print summary.IC
#' @keywords internal

print.summary.IC <- function(x, ...) {
  Group <- region <- IC <- NULL
  title <- 'Individual contributions'
  message('\n', title, '\n', rep('-', getOption('width') / 2))
  cat('Method: ', x$method, '\n')
  cat('Level: ', x$level, '\n\n')

  if (x$level == 'regional') {
    cat('Number of outliers per region: (sorted in descending order)\n')
    print(sort(x$outliers$region, decreasing=TRUE))
    cat('\n')
    print(x$DT.sum[, .SD[1, !c('Study.ID', 'RC')], by=list(Group, region)])
  } else {
    cat('"Outliers"\n')
    print(x$outliers$DT[order(Group, -IC)])
  }
  invisible(x)
}


#' Plot regional contributions estiamtes
#'
#' The \code{plot} method is only valid for \emph{regional} contribution
#' estimates, and plots the average regional contribution for each
#' vertex/region.
#'
#' @param x A \code{IC} object
#' @param plot.type Character string indicating the type of plot; the default is
#'   to plot the mean (along with standard errors)
#' @export
#' @importFrom ggrepel geom_text_repel
#' @method plot IC
#' @rdname individ_contrib

plot.IC <- function(x, plot.type=c('mean', 'smooth', 'boxplot'), region=NULL, ...) {
  RC <- Group <- avg <- se <- ind <- mark <- IC <- Study.ID <- NULL
  DT <- copy(x$DT)
  if (x$level == 'regional') {
    if (is.null(region)) {
      regions <- DT[, levels(region)]
    } else {
      regions <- region
    }
    txtsize <- ifelse(length(regions) > 50, 6, 9)

    plot.type <- match.arg(plot.type)
    if (plot.type == 'boxplot') {
      p <- ggplot(DT[region %in% regions], aes(x=region, y=RC)) +
        geom_boxplot(aes(fill=Group, group=interaction(Group, region)))

    } else {
      if (plot.type == 'smooth') {
        p <- ggplot(DT[region %in% regions], aes(x=region, col=Group, group=Group)) +
          stat_smooth(method='loess', aes(y=RC))

      } else if (plot.type == 'mean') {
        DT[, avg := mean(RC), by=list(Group, region)]
        DT[, se := sd(RC) / sqrt(.N), by=list(Group, region)]
        p <- ggplot(DT[region %in% regions], aes(x=region, col=Group, group=Group)) +
          geom_line(aes(y=avg)) +
          geom_ribbon(aes(ymin=avg-se, ymax=avg+se, fill=Group), alpha=0.5)
      }
    }
    p <- p +
      theme(legend.position='bottom',
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            axis.text.x=element_text(size=txtsize, angle=45, vjust=0.5),
            plot.title=element_text(hjust=0.5, face='bold')) +
      labs(x='Region', y='Regional contribution',
           title=paste0('Regional contributions, ', tolower(x$method), ' method'))

  } else {
    txtsize <- ifelse(nrow(DT) > 50, 9, 12)
    DT[, ind := as.character(.SD[, .I])]
    DT[, mark := ifelse(IC > mean(IC) + 2*sd(IC), 1, 0)]
    DT[mark == 0, ind := '']
    DT[, mark := as.factor(mark)]
    p <- ggplot(DT, aes(x=Study.ID, y=IC, col=Group)) +
      geom_text_repel(aes(label=ind), size=3) +
      geom_point(aes(shape=mark, size=mark)) +
      scale_shape_manual(values=c(20, 17)) +
      scale_size_manual(values=c(2, 3)) +
      labs(x='Subject', y='Individual contribution',
           title=paste0('Individual contributions, ', tolower(x$method), ' method')) +
      theme(legend.position='none',
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            axis.text.x=element_text(size=txtsize, angle=45, vjust=0.5),
            plot.title=element_text(hjust=0.5, face='bold'))
  }
  return(p)
}
