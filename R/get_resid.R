#' Linear model residuals in structural covariance networks
#'
#' \code{get.resid} runs linear models across brain regions listed in a
#' \code{data.table} (e.g. cortical thickness), adjusting for variables in
#' \code{covars} (e.g. age, sex, etc.), and calculates the
#' \emph{externally Studentized} (or \emph{leave-one-out}) residuals.
#'
#' You can choose to run models for each of your subject groups separately or
#' combined (the default) via the \code{method} argument. You may also choose
#' whether or not to include the mean, per-hemisphere structural measure in the
#' models. Finally, you can specify variables that are present in \code{covars}
#' but you would like to exclude from the models.
#'
#' @param dt.vol A \code{data.table} containing all the volumetric measure of
#'   interest (i.e., the object \code{lhrh} as ouptut by
#'   \code{\link{brainGraph_init}})
#' @param covars A \code{data.table} of the covariates of interest
#' @param method Character string indicating whether to test models for subject
#'   groups separately or combined (default: \code{comb.groups})
#' @param use.mean Logical should we control for the mean hemispheric brain
#'   value (e.g. mean LH/RH cortical thickness) (default: \code{FALSE})
#' @param exclude Character vector of covariates to exclude (default:
#'   \code{NULL})
#' @param ... Arguments passed to \code{\link{brainGraph_GLM_design}} (optional)
#' @export
#'
#' @return \code{get.resid} - an object of class \code{brainGraph_resids} with
#'   elements:
#'   \item{X}{The \emph{design matrix}}
#'   \item{method}{The input argument \code{method}}
#'   \item{use.mean}{The input argument \code{use.mean}}
#'   \item{all.dat.tidy}{The tidied \code{data.table} of volumetric data (e.g.,
#'     mean regional cortical thickness) and covariates, along with
#'     \emph{resids} column added}
#'   \item{resids.all}{The "wide" \code{data.table} of residuals}
#'   \item{groups}{Group names}
#' @seealso \code{\link{rstudent}}
#' @family Structural covariance network functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

get.resid <- function(dt.vol, covars, method=c('comb.groups', 'sep.groups'),
                      use.mean=FALSE, exclude=NULL, ...) {
  region <- resids <- Group <- Study.ID <- value <- NULL

  stopifnot('Group' %in% names(covars))
  if (!'Study.ID' %in% names(covars)) covars$Study.ID <- as.character(seq_len(nrow(covars)))
  method <- match.arg(method)
  groups <- covars[, levels(factor(Group))]
  exclude <- c('Study.ID', exclude)
  DT.cov <- merge(covars, dt.vol, by='Study.ID')
  DT.m <- melt(DT.cov, id.vars=names(covars), variable.name='region')
  setkey(DT.m, region, Study.ID)

  if (isTRUE(use.mean)) {
    if (length(grep('^l.*', names(dt.vol))) == 0) {
      lh <- '.*\\.L$'
      rh <- '.*\\.R$'
    } else {
      lh <- '^l.*'
      rh <- '^r.*'
    }
    mean.lh <- dt.vol[, rowMeans(.SD), .SDcols=names(dt.vol)[grep(lh, names(dt.vol))], by=Study.ID]
    mean.rh <- dt.vol[, rowMeans(.SD), .SDcols=names(dt.vol)[grep(rh, names(dt.vol))], by=Study.ID]
    covars.lh <- cbind(covars, mean.lh)
    covars.rh <- cbind(covars, mean.rh)

    if (method == 'comb.groups') {
      lhvars <- get_lm_vars(covars.lh, exclude, ...)
      rhvars <- get_lm_vars(covars.rh, exclude, ...)

      DT.m[grep(lh, region), resids := rstudent_mat(lhvars, value), by=region]
      DT.m[grep(rh, region), resids := rstudent_mat(rhvars, value), by=region]
      X <- list(lh=lhvars$X, rh=rhvars$X)
    } else {
      covars.lh <- split(covars.lh, by='Group')
      covars.rh <- split(covars.rh, by='Group')
      DT.m <- split(DT.m, by='Group')
      X <- sapply(groups, function(x) NULL)
      for (g in groups) {
        lhvars <- get_lm_vars(covars.lh[[g]], exclude, ...)
        rhvars <- get_lm_vars(covars.rh[[g]], exclude, ...)

        DT.m[[g]][grep(lh, region), resids := rstudent_mat(lhvars, value), by=region]
        DT.m[[g]][grep(rh, region), resids := rstudent_mat(rhvars, value), by=region]
        X[[g]] <- list(lh=lhvars$X, rh=rhvars$X)
      }
      DT.m <- rbindlist(DT.m)
    }

  } else {
    if (method == 'comb.groups') {
      lmvars <- get_lm_vars(covars, exclude, ...)
      DT.m[, resids := rstudent_mat(lmvars, value), by=region]
      X <- lmvars$X
    } else {
      covars <- split(covars, by='Group')
      DT.m <- split(DT.m, by='Group')
      X <- sapply(groups, function(x) NULL)
      for (g in groups) {
        lmvars <- get_lm_vars(covars[[g]], exclude, ...)
        DT.m[[g]][, resids := rstudent_mat(lmvars, value), by=region]
        X[[g]] <- lmvars$X
      }
      DT.m <- rbindlist(DT.m)
      covars <- rbindlist(covars)
    }
  }

  # Return data to "wide" format with just the residuals
  resids.all <- dcast(DT.m, 'Study.ID + Group ~ region', value.var='resids')
  setkey(resids.all, Group, Study.ID)

  out <- list(X=X, method=method, use.mean=use.mean, all.dat.tidy=DT.m, resids.all=resids.all, groups=groups)
  class(out) <- c('brainGraph_resids', class(out))
  return(out)
}

#' Get some variables for LM
#'
#' @inheritParams get.resid
#' @keywords internal
#' @return A list containing:
#'   \item{X}{The design matrix}
#'   \item{lev}{The leverage}
#'   \item{n}{The number of observations}
#'   \item{p}{The number of parameters}

get_lm_vars <- function(covars, exclude, ...) {
  X <- brainGraph_GLM_design(covars[, !exclude, with=F], ...)
  H <- X %*% solve(crossprod(X)) %*% t(X)
  lev <- diag(H)
  n <- nrow(X)
  p <- ncol(X)
  return(list(X=X, lev=lev, n=n, p=p))
}

#' Calculate studentized residuals with matrix input
#'
#' @param lmvars List containing: \code{X} (design matrix); \code{lev}
#'   (leverage); \code{n} (num. observations); \code{p} (num. parameters)
#' @param y Numeric vector; the outcome variable
#' @keywords internal

rstudent_mat <- function(lmvars, y) {
  est <- fastLmPure(lmvars$X, y, method=2)
  var.hat <- rep(0, lmvars$n)
  for (i in seq_len(lmvars$n)) var.hat[i] <- (1 / (lmvars$n - lmvars$p - 1)) * crossprod(est$residuals[-i])
  resids <- est$residuals / (sqrt(var.hat * (1 - lmvars$lev)))
}

#' Indexing for structural covariance residuals
#'
#' The \code{[} method will let you reorder or subset residuals based on a given
#' numeric vector. However, this is used in bootstrap and permutation analysis
#' and should generally not be called directly by the user.
#'
#' @param x A \code{brainGraph_resids} object
#' @param i Numeric vector of the indices
#' @param g Character string indicating the group (default: \code{NULL})
#' @export
#' @method [ brainGraph_resids
#'
#' @name Extract.brainGraph_resids
#' @rdname get.resid

`[.brainGraph_resids` <- function(x, i, g=NULL) {
  Group <- Study.ID <- NULL
  if (!is.null(g)) x$resids.all <- x$resids.all[g]
  if (missing(i)) i <- seq_len(nrow(x$resids.all))
  x$resids.all <- x$resids.all[i]
  setkey(x$resids.all, Group, Study.ID)
  x$all.dat.tidy <- x$all.dat.tidy[Study.ID %in% x$resids.all$Study.ID]
  return(x)
}

#' Print a summary of residuals for structural covariance data
#'
#' The \code{summary} method prints the number of outliers per region, and the
#' number of times a given subject was an outlier (i.e., across regions).
#'
#' @param object A \code{brainGraph_resids} object
#' @param regions Character vector of region(s) to focus on; default behavior is
#'   to show summary for all regions
#' @export
#' @method summary brainGraph_resids
#' @return \code{\link{summary.brainGraph_resids}} returns a list with two
#'   data tables, one of the residuals, and one of only the outlier regions
#' @rdname get.resid

summary.brainGraph_resids <- function(object, regions=NULL, ...) {
  region <- resids <- Study.ID <- NULL
  myregions <- regions
  if (is.null(myregions)) myregions <- object$all.dat.tidy[, levels(region)]
  DT <- droplevels(object$all.dat.tidy[region %in% myregions,
                   c('Study.ID', 'Group', 'region', 'resids')])
  setkey(DT, region, resids)
  outliers <- DT[, .SD[abs(resids) > abs(mean(resids) + 2*sd(resids))], by=region]
  outliers.reg <- outliers[, .N, by=region]
  outliers.reg.vec <- structure(outliers.reg$N, names=as.character(outliers.reg$region))

  outliers.sub <- outliers[, .N, by=Study.ID]
  outliers.sub.vec <- structure(outliers.sub$N, names=as.character(outliers.sub$Study.ID))

  out <- list(DT.sum=DT,
              outliers=list(DT=outliers, region=outliers.reg.vec, subject=outliers.sub.vec))
  class(out) <- c('summary.brainGraph_resids', class(out))
  return(out)
}

#' @aliases summary.brainGraph_resids
#' @method print summary.brainGraph_resids
#' @keywords internal

print.summary.brainGraph_resids <- function(x, ...) {
  title <- 'Structural covariance residuals'
  message('\n', title, '\n', rep('-', getOption('width') / 2))
  cat('Number of outliers per region: (sorted in descending order)\n')
  print(sort(x$outliers$region, decreasing=TRUE))

  cat('\n\nNumber of times each subject was an outlier: (sorted in descending order)\n')
  print(sort(x$outliers$subject, decreasing=TRUE))
  cat('\n')
  invisible(x)
}

#' Plot model residuals for each brain region
#'
#' Check the model residuals for each brain region in a structural covariance
#' analysis. It shows a \emph{qqplot} of the studentized residuals, as output
#' from \code{\link{get.resid}}.
#'
#' @param x A \code{brainGraph_resids} object (from \code{\link{get.resid}})
#' @param cols Logical indicating whether to color by group (default: FALSE)
#' @inheritParams summary.brainGraph_resids
#' @export
#' @method plot brainGraph_resids
#' @importFrom ggrepel geom_text_repel
#'
#' @return A list of \code{\link[ggplot2]{ggplot}} objects
#'
#' @family Structural covariance network functions
#' @seealso \code{\link[stats]{qqnorm}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' myresids <- get.resids(lhrh, covars)
#' residPlots <- plot(myresids, cols=TRUE)
#'
#' ## Save as a multi-page PDF
#' ml <- marrangeGrob(residPlots, nrow=3, ncol=3)
#' ggsave('residuals.pdf', ml)
#' }

plot.brainGraph_resids <- function(x, regions=NULL, cols=FALSE, ...) {
  region <- ind <- mark <- Study.ID <- resids <- NULL
  DT <- summary(x, regions=regions)$DT.sum
  regions <- DT[, levels(region)]

  setkey(DT, region, Study.ID)
  DT[, ind := as.character(.SD[, .I]), by=region]
  setkey(DT, region, resids)
  DT[, x := qnorm(ppoints(resids)), by=region]
  DT[, mark := ifelse(abs(resids) > abs(mean(resids) + 2*sd(resids)), 1, 0), by=region]
  DT[mark == 0, ind := '']
  DT[, mark := as.factor(mark)]

  # Local function to plot for a single region
  plot_single <- function(DT.resids, cols) {
    Group <- resids <- NULL
    p <- ggplot(DT.resids, aes(x=x, y=resids, col=Group)) +
      geom_text_repel(aes(label=ind), size=3) +
      geom_point(aes(shape=mark, size=mark)) +
      geom_line(aes(x=x, y=x), col='gray50') +
      scale_shape_manual(values=c(20, 17)) +
      scale_size_manual(values=c(2, 3)) +
      theme(plot.title=element_text(hjust=0.5, size=10, face='bold'),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title=paste0('Normal Q-Q: ', DT.resids[, unique(region)]),
           x='Theoretical Quantiles', y='Sample Quantiles')
    if (!isTRUE(cols)) {
      p <- p + scale_color_manual(values=rep('black', DT.resids[, length(unique(Group))]))
    }
    return(p)
  }

  p.all <- sapply(regions, function(x) NULL)
  for (z in regions) {
    p.all[[z]] <- plot_single(DT[region == z], cols)
  }
  return(p.all)
}
