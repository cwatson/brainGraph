#' Linear model residuals across brain regions
#'
#' Runs linear models across brain regions listed in a \code{data.table} (e.g.
#' cortical thickness), adjusting for variables in \code{covars} (e.g. age, sex,
#' etc.). It also calculates the \emph{externally Studentized} (or
#' \emph{leave-one-out}) residuals.
#'
#' @param dt.vol A \code{data.table} containing all the volumetric measure of
#'   interest (i.e., the object \code{lhrh} as ouptut by
#'   \code{\link{brainGraph_init}})
#' @param covars A \code{data.table} of the covariates of interest
#' @param use.mean Logical should we control for the mean hemispheric brain
#'   value (e.g. mean LH/RH cortical thickness) (default: \code{FALSE})
#' @param exclude Character vector of covariates to exclude (default:
#'   \code{NULL})
#' @param ... Arguments passed to \code{\link{brainGraph_GLM_design}} (optional)
#' @export
#'
#' @return An object of class \code{brainGraph_resids} with elements:
#'   \item{X}{The \emph{design matrix}}
#'   \item{all.dat.tidy}{The tidied \code{data.table} of volumetric data (e.g.,
#'     mean regional cortical thickness) and covariates, along with
#'     \emph{resids} column added}
#'   \item{resids.all}{The "wide" \code{data.table} of residuals}
#' @seealso \code{\link{rstudent}}
#' @family Structural covariance network functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

get.resid <- function(dt.vol, covars, use.mean=FALSE, exclude=NULL, ...) {
  region <- resids <- Group <- Study.ID <- value <- NULL

  exclude <- c('Study.ID', exclude)
  DT.cov <- merge(covars, dt.vol, by='Study.ID')
  DT.m <- melt(DT.cov, id.vars=names(covars), variable.name='region')
  setkey(DT.m, region, Study.ID)

  # Adjust by mean hemispheric values
  if (isTRUE(use.mean)) {
    if (length(grep('^l.*', names(dt.vol))) == 0) {
      lh.string <- '.*\\.L$'
      rh.string <- '.*\\.R$'
    } else {
      lh.string <- '^l.*'
      rh.string <- '^r.*'
    }
    mean.lh <- dt.vol[, rowMeans(.SD), .SDcols=names(dt.vol)[grep(lh.string, names(dt.vol))], by=Study.ID]
    mean.rh <- dt.vol[, rowMeans(.SD), .SDcols=names(dt.vol)[grep(rh.string, names(dt.vol))], by=Study.ID]
    covars.lh <- cbind(covars, mean.lh)
    covars.rh <- cbind(covars, mean.rh)
    X.lh <- brainGraph_GLM_design(covars.lh[, !exclude, with=F], ...)
    H.lh <- X.lh %*% solve(crossprod(X.lh)) %*% t(X.lh)
    lev.lh <- diag(H.lh)
    n.lh <- nrow(X.lh)
    p.lh <- ncol(X.lh)

    X.rh <- brainGraph_GLM_design(covars.rh[, !exclude, with=F], ...)
    H.rh <- X.rh %*% solve(crossprod(X.rh)) %*% t(X.rh)
    lev.rh <- diag(H.rh)
    n.rh <- nrow(X.rh)
    p.rh <- ncol(X.rh)

    DT.m[grep(lh.string, region),
              resids := rstudent_mat(X.lh, value, lev.lh, n.lh, p.lh), by=region]
    DT.m[grep(rh.string, region),
              resids := rstudent_mat(X.rh, value, lev.rh, n.rh, p.rh), by=region]
    X <- list(lh=X.lh, rh=X.rh)

  } else {
    if ('mean.lh' %in% names(covars)) exclude <- c(exclude, 'mean.lh')
    if ('mean.rh' %in% names(covars)) exclude <- c(exclude, 'mean.rh')
    X <- brainGraph_GLM_design(covars[, !exclude, with=F], ...)
    H <- X %*% solve(crossprod(X)) %*% t(X)
    lev <- diag(H)
    n <- nrow(X)
    p <- ncol(X)
    DT.m[, resids := rstudent_mat(X, value, lev, n, p), by=region]
  }

  # Return data to "wide" format with just the residuals
  all.dat.wide <- dcast(DT.m,
                        paste(paste(names(covars), collapse='+'), '~ region'),
                        value.var='resids')
  resids.all <- all.dat.wide[, !setdiff(names(covars), c('Study.ID', 'Group')), with=F]
  setkey(resids.all, Group, Study.ID)

  out <- list(X=X, all.dat.tidy=DT.m, resids.all=resids.all)
  class(out) <- c('brainGraph_resids', class(out))
  return(out)
}

#' Calculate studentized residuals with matrix input
#'
#' @param X Numeric matrix; the design matrix
#' @param y Numeric vector; the outcome variable
#' @param lev Numeric vector; the leverage of the "hat" matrix
#' @param n Integer; number of rows in \code{X}
#' @param p Integer; number of columns in \code{X}
#' @keywords internal

rstudent_mat <- function(X, y, lev, n, p) {
  est <- fastLmPure(X, y, method=2)
  var.hat <- rep(0, n)
  for (i in seq_len(n)) var.hat[i] <- (1 / (n - p - 1)) * crossprod(est$residuals[-i])
  resids <- est$residuals / (sqrt(var.hat * (1 - lev)))
}

#' Print a summary of residuals for structural covariance data
#'
#' @param object A \code{brainGraph_resids} object
#' @param regions Character vector of region(s) to focus on; default behavior is
#'   to show summary for all regions
#' @param ... Unused
#' @export
#' @method summary brainGraph_resids
#'
#' @return A list with the full \code{data.table} of residuals, in addition to a
#'   \code{data.table} of only the outliers

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

  out <- list(DT=DT,
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
#' ## First get the resid's and store the plots as an object
#' myresids <- get.resids(lhrh, covars)
#' p.resids <- plot(myresids, cols=TRUE)
#'
#' ## Open a new plot device for each set of 9
#' lapply(p.resids, function(x) {dev.new(); print(x)})
#'
#' ## Save as a multi-page PDF; requires \code{gridExtra}
#' ml <- marrangeGrob(p.resids, nrow=1, ncol=1)
#' ggsave('residuals.pdf', ml)
#' }

plot.brainGraph_resids <- function(x, regions=NULL, cols=FALSE, ...) {
  region <- ind <- mark <- Study.ID <- resids <- NULL
  res.sum <- summary(x, regions=regions)
  DT <- res.sum$DT
  kNumRegions <- DT[, nlevels(region)]

  a <- ifelse(kNumRegions < 9, 1, kNumRegions %/% 9)
  b <- kNumRegions %% 9

  setkey(DT, region, Study.ID)
  DT[, ind := as.character(.SD[, .I]), by=region]
  setkey(DT, region, resids)
  DT[, x := qnorm(ppoints(resids)), by=region]
  DT[, mark := ifelse(abs(resids) > abs(mean(resids) + 2*sd(resids)), 1, 0), by=region]
  DT[mark == 0, ind := '']
  DT[, mark := as.factor(mark)]

  ggQQ <- function(DT.resids, cols) {
    Group <- resids <- NULL
    p <- ggplot(DT.resids, aes(x=x, y=resids)) +
      labs(x='Theoretical Quantiles', y='Sample Quantiles')
    if (isTRUE(cols)) {
      p <- p +
        geom_text_repel(aes(x=x, y=resids, label=ind, col=Group), size=3) +
        geom_line(aes(x=x, y=x), col='gray50') +
        geom_point(aes(col=Group, shape=mark, size=mark)) +
        scale_shape_manual(values=c(20, 8)) +
        scale_size_manual(values=c(1.5, 3)) +
        theme(legend.position='none')
    } else {
      p <- p +
        geom_point() +
        geom_line(aes(x=x, y=x))
    }
    p <- p + facet_wrap( ~ region, nrow=3, ncol=3, scales='free')
    return(p)
  }

  all.p <- vector('list', length=(a + (!b == 0)))
  for (j in seq_len(a)) {
    N1 <- 9 * (j - 1) + 1
    N2 <- min(N1 + 8, kNumRegions)
    all.p[[j]] <- ggQQ(DT[levels(region)[N1:N2]], cols)
  }

  if (kNumRegions > 9 & ! b == 0) {
    all.p[[j+1]] <- ggQQ(DT[levels(region)[(N2+1):kNumRegions]], cols)
  }
  return(all.p)
}
