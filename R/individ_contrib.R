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
#'   from group-based structural correlations networks}. NeuroImage, 120:274-284.
#'   doi:10.1016/j.neuroimage.2015.07.006

loo <- function(resids, corrs, level=c('global', 'regional')) {
  Group <- Study.ID <- i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'))
  level <- match.arg(level)
  group.vec <- resids$resids.all$Group
  group.num <- as.integer(group.vec)
  group.vec <- as.character(group.vec)
  if (level == 'global') {
    IC <- foreach (i=seq_len(nrow(resids$resids.all)), .combine='c') %dopar% {
      resids.excl <- resids[-i]
      new.corrs <- corr.matrix(resids.excl[group.vec[i]], densities=0.1)

      1 - mantel.rtest(as.dist(corrs[[group.num[i]]]$R),
                       as.dist(new.corrs[[1]]$R),
                       nrepet=1e3)$obs
    }

    DT <- data.table(resids$resids.all[, list(Study.ID, Group)], IC=IC)
  } else if (level == 'regional') {
    RC <- foreach (i=seq_len(nrow(resids$resids.all)), .combine='rbind') %dopar% {
      resids.excl <- resids[-i]
      new.corrs <- corr.matrix(resids.excl[group.vec[i]], densities=0.1)
      colSums(abs(corrs[[group.num[i]]]$R - new.corrs[[1]]$R))
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
#' @param corr.mat Numeric; correlation matrix of the \emph{control} group
#' @param control.value Integer or character string specifying the control group
#'   (default: 1)
#' @export
#' @importFrom ade4 mantel.rtest
#'
#' @aliases aop
#' @rdname individ_contrib
#' @examples
#' \dontrun{
#' IC <- aop(resids.all, corrs[[1]]$R)
#' RC <- aop(resids.all, corrs[[1]]$R, level='regional')
#' }

aop <- function(resids, corr.mat, level=c('global', 'regional'), control.value=1) {
  Group <- Study.ID <- i <- NULL
  stopifnot(inherits(resids, 'brainGraph_resids'))

  groups <- resids$resids.all[, levels(Group)]
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
        setkey(resids.aop$resids.all, Group)
        new.corr <- corr.matrix(resids.aop, densities=0.1)[[1]]$R
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
        setkey(resids.aop$resids.all, Group)
        new.corr <- corr.matrix(resids.aop, densities=0.1)[[1]]$R
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
  DT.sum <- copy(object$DT)
  if (is.null(region)) {
    regions <- DT.sum[, levels(region)]
  } else {
    regions <- region
  }
  DT.sum[, avg := mean(RC), by=list(Group, region)]
  DT.sum[, stdev := sd(RC), by=list(Group, region)]
  DT.sum[, se := stdev / .N, by=list(Group, region)]

  outliers <- DT.sum[, .SD[RC > mean(RC) + 2 * stdev], by=list(Group, region)]
  outliers.reg <- outliers[, .N, by=region]
  outliers.reg.vec <- structure(outliers.reg$N, names=as.character(outliers.reg$region))
  object$DT.sum <- DT.sum
  object$outliers <- list(DT=outliers, region=outliers.reg.vec)
  class(object) <- c('summary.IC', class(object))
  return(object)
}

#' @aliases summary.IC
#' @method print summary.IC
#' @keywords internal

print.summary.IC <- function(x, ...) {
  title <- 'Individual contributions'
  message('\n', title, '\n', rep('-', getOption('width') / 2))
  cat('Method: ', x$method, '\n')
  cat('Level: ', x$level, '\n\n')

  cat('Number of outliers per region: (sorted in descending order)\n')
  print(sort(x$outliers$region, decreasing=TRUE))
  cat('\n')
  print(x$DT.sum[, .SD[1, !c('Study.ID', 'RC')], by=list(Group, region)])
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
#' @method plot IC
#' @rdname individ_contrib

plot.IC <- function(x, plot.type=c('mean', 'smooth', 'boxplot'), region=NULL, ...) {
  stopifnot(x$level == 'regional')
  DT <- copy(x$DT)
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
      DT[, se := sd(RC) / .N, by=list(Group, region)]
      p <- ggplot(DT[region %in% regions], aes(x=region, col=Group, group=Group)) +
        geom_line(aes(y=avg)) +
        geom_ribbon(aes(ymin=avg-se, ymax=avg+se, fill=Group), alpha=0.5)
    }
  }
  p <- p +
    theme(legend.position='bottom',
          axis.text.x=element_text(size=txtsize, angle=45, vjust=0.5),
          plot.title=element_text(hjust=0.5, face='bold')) +
    labs(x='Region', y='Regional contribution',
         title=paste0('Regional contributions, ', tolower(x$method), ' method'))
  return(p)
}
