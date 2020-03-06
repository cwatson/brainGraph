#' Linear model residuals in structural covariance networks
#'
#' \code{get.resid} runs linear models across brain regions listed in a
#' \code{data.table} (e.g., cortical thickness), adjusting for variables in
#' \code{covars} (e.g. age, sex, etc.), and calculates the
#' \emph{externally Studentized} (or \emph{leave-one-out}) residuals.
#'
#' You can choose to run models for each of your subject groups separately or
#' combined (the default) via the \code{method} argument. You may also choose
#' whether or not to include the mean, per-hemisphere structural measure in the
#' models. Finally, you can specify variables that are present in \code{covars}
#' but you would like to exclude from the models.
#'
#' If you do not explicitly specify the atlas name, then it will be guessed from
#' the size of your data. This could cause problems if you are using a custom
#' atlas, with or without the same number of regions as a dataset in the
#' package.
#'
#' @param dt.vol A \code{data.table} containing all the volumetric measure of
#'   interest (i.e., the object \code{lhrh} as output by
#'   \code{\link{import_scn}})
#' @param covars A \code{data.table} of the covariates of interest
#' @param method Character string indicating whether to test models for subject
#'   groups separately or combined. Default: \code{comb.groups}
#' @param use.mean Logical should we control for the mean hemispheric brain
#'   value (e.g. mean LH/RH cortical thickness). Default: \code{FALSE}
#' @param exclude.cov Character vector of covariates to exclude. Default:
#'   \code{NULL}
#' @param atlas Character string indicating the brain atlas
#' @param ... Arguments passed to \code{\link{brainGraph_GLM_design}} (optional)
#' @export
#'
#' @return \code{get.resid} - an object of class \code{brainGraph_resids} with
#'   elements:
#'   \item{X}{The \emph{design matrix}, if using default arguments. If
#'   \code{use.mean=TRUE} then it will be a \emph{named list} with a separate
#'   matrix for the left and right hemispheres. If \code{method='sep.groups'}, a
#'   nested named list for each group and hemisphere.}
#'   \item{method}{The input argument \code{method}}
#'   \item{use.mean}{The input argument \code{use.mean}}
#'   \item{all.dat.long}{A \code{data.table} in \dQuote{long} format of the
#'     input structural data merged with the covariates, and a \emph{resids}
#'     column added}
#'   \item{resids.all}{The \dQuote{wide} \code{data.table} of residuals}
#'   \item{Group}{Group names}
#'   \item{atlas}{The atlas name}
#' @name Residuals
#' @rdname residuals
#' @seealso \code{\link[stats]{influence.measures}}, \code{\link[stats]{qqnorm}}
#' @family Structural covariance network functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

get.resid <- function(dt.vol, covars, method=c('comb.groups', 'sep.groups'),
                      use.mean=FALSE, exclude.cov=NULL, atlas=NULL, ...) {
  Region <- resids <- value <- NULL

  gID <- getOption('bg.group')
  stopifnot(hasName(covars, gID))
  grps <- covars[, levels(factor(get(gID)))]
  sID <- getOption('bg.subject_id')
  if (!hasName(covars, sID)) covars[, eval(sID) := seq_len(dim(covars)[1L])]
  covars[, eval(sID) := check_sID(get(sID))]
  DT.cov <- merge(covars, dt.vol, by=sID)
  DT.m <- melt(DT.cov, id.vars=names(covars), variable.name='Region')
  setkeyv(DT.m, c('Region', sID))

  method <- match.arg(method)
  if (isTRUE(use.mean)) {
    if (length(grep('^l.*', names(dt.vol))) == 0) {
      lh <- '.*\\.L$'
      rh <- '.*\\.R$'
    } else {
      lh <- '^l.*'
      rh <- '^r.*'
    }
    mean.lh <- dt.vol[, rowMeans(.SD), .SDcols=patterns(lh)]
    mean.rh <- dt.vol[, rowMeans(.SD), .SDcols=patterns(rh)]
    covars.lh <- cbind(covars, mean.lh)
    covars.rh <- cbind(covars, mean.rh)

    if (method == 'comb.groups') {
      lhvars <- get_lm_vars(covars.lh, exclude.cov, ...)
      rhvars <- get_lm_vars(covars.rh, exclude.cov, ...)

      DT.m[grep(lh, Region), resids := rstudent_mat(lhvars, value), by=Region]
      DT.m[grep(lh, Region), resids := resids / sqrt(1 / lhvars$df)]
      DT.m[grep(rh, Region), resids := rstudent_mat(rhvars, value), by=Region]
      DT.m[grep(rh, Region), resids := resids / sqrt(1 / rhvars$df)]
      X <- list(lh=lhvars$X, rh=rhvars$X)
    } else {
      covars.lh <- split(covars.lh, by=gID)
      covars.rh <- split(covars.rh, by=gID)
      DT.m <- split(DT.m, by=gID)
      X <- setNames(vector('list', length(grps)), grps)
      for (g in grps) {
        lhvars <- get_lm_vars(covars.lh[[g]], exclude.cov, ...)
        rhvars <- get_lm_vars(covars.rh[[g]], exclude.cov, ...)

        DT.m[[g]][grep(lh, Region), resids := rstudent_mat(lhvars, value), by=Region]
        DT.m[[g]][grep(lh, Region), resids := resids / sqrt(1 / lhvars$df)]
        DT.m[[g]][grep(rh, Region), resids := rstudent_mat(rhvars, value), by=Region]
        DT.m[[g]][grep(rh, Region), resids := resids / sqrt(1 / rhvars$df)]
        X[[g]] <- list(lh=lhvars$X, rh=rhvars$X)
      }
      DT.m <- rbindlist(DT.m)
    }

  } else {
    if (method == 'comb.groups') {
      lmvars <- get_lm_vars(covars, exclude.cov, ...)
      DT.m[, resids := rstudent_mat(lmvars, value), by=Region]
      DT.m[, resids := resids / sqrt(1 / lmvars$df)]
      X <- lmvars$X
    } else {
      covars <- split(covars, by=gID)
      DT.m <- split(DT.m, by=gID)
      X <- setNames(vector('list', length(grps)), grps)
      for (g in grps) {
        lmvars <- get_lm_vars(covars[[g]], exclude.cov, ...)
        DT.m[[g]][, resids := rstudent_mat(lmvars, value), by=Region]
        DT.m[[g]][, resids := resids / sqrt(1 / lmvars$df)]
        X[[g]] <- lmvars$X
      }
      DT.m <- rbindlist(DT.m)
      covars <- rbindlist(covars)
    }
  }

  # Return data to "wide" format with just the residuals
  resids.all <- dcast(DT.m, paste(sID, '+', gID, '~ Region'), value.var='resids')
  setkeyv(resids.all, c(gID, sID))

  out <- list(X=X, method=method, use.mean=use.mean, all.dat.long=DT.m,
              resids.all=resids.all, Group=grps, atlas=NULL)
  out$atlas <- if (is.null(atlas)) guess_atlas(dt.vol) else atlas

  class(out) <- c('brainGraph_resids', class(out))
  return(out)
}

#' Get some variables for LM
#'
#' @inheritParams Residuals
#' @keywords internal
#' @return A list containing:
#'   \item{X}{The design matrix}
#'   \item{lev}{The leverage}
#'   \item{n}{The number of observations}
#'   \item{p}{The number of parameters}

get_lm_vars <- function(covars, exclude.cov, ...) {
  if (!is.null(exclude.cov)) covars <- covars[, !exclude.cov, with=FALSE]
  X <- brainGraph_GLM_design(covars, ...)
  XtX <- crossprod(X)
  U <- chol.default(XtX)
  Z <- forwardsolve(t(U), t(X))
  lev <- colSums(Z^2)
  dims <- dim(X)
  return(list(X=X, XtX=XtX, lev=lev, df=dims[1L]-dims[2L]-1L))
}

#' Calculate studentized residuals with matrix input
#'
#' @param lmvars List containing: \code{X} (design matrix); \code{XtX}
#'   (cross-product of the design matrix); \code{lev} (leverage); \code{df}
#'   (model degrees of freedom minus 1)
#' @param y Numeric vector; the outcome variable
#' @keywords internal

rstudent_mat <- function(lmvars, y) {
  XtY <- crossprod(lmvars$X, y)
  beta <- solve.default(lmvars$XtX, XtY)
  res <- c(y - lmvars$X %*% beta)

  var.hat <- c(crossprod(res)) - res^2
  resids <- res / (sqrt(var.hat * (1 - lmvars$lev)))
}

#' Indexing for structural covariance residuals
#'
#' The \code{[} method reorders or subsets residuals based on a given
#' numeric vector. However, this is used in bootstrap and permutation analysis
#' and should generally not be called directly by the user.
#'
#' @param x,object A \code{brainGraph_resids} object
#' @param i Numeric vector of the indices
#' @param g Character string indicating the group. Default: \code{NULL}
#' @export
#'
#' @name Extract.brainGraph_resids
#' @rdname residuals

`[.brainGraph_resids` <- function(x, i, g=NULL) {
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  if (!is.null(g)) x$resids.all <- droplevels(x$resids.all[g])
  if (missing(i)) i <- seq_len(dim(x$resids.all)[1L])
  x$resids.all <- x$resids.all[i]
  x$Group <- x$resids.all[, levels(factor(get(gID)))]
  setkeyv(x$resids.all, c(gID, sID))
  x$all.dat.long <- x$all.dat.long[get(sID) %in% case.names(x)]
  return(x)
}

#' Print a summary of residuals for structural covariance data
#'
#' The \code{summary} method prints the number of outliers per region, and the
#' number of times a given subject was an outlier (i.e., across regions).
#'
#' @param region Character vector of region(s) to focus on; default behavior is
#'   to show summary for all regions
#' @export
#' @return \code{\link{summary.brainGraph_resids}} returns a list with two
#'   data tables, one of the residuals, and one of only the outlier regions
#' @rdname residuals

summary.brainGraph_resids <- function(object, region=NULL, ...) {
  Region <- resids <- NULL
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  if (is.null(region)) region <- region.names(object)
  DT <- droplevels(object$all.dat.long[Region %in% region,
                                       c(sID, gID, 'Region', 'resids'), with=FALSE])
  setkey(DT, Region, resids)
  outliers <- DT[, .SD[abs(resids) > abs(mean(resids) + 2*sd(resids))], by=Region]
  outliers.reg <- outliers[, .N, by=Region]
  outliers.reg.vec <- structure(outliers.reg$N, names=as.character(outliers.reg$Region))

  outliers.sub <- outliers[, .N, by=sID]
  outliers.sub.vec <- structure(outliers.sub$N, names=as.character(outliers.sub[, get(sID)]))

  out <- list(Region=region, DT.sum=DT,
              outliers=list(DT=outliers, Region=outliers.reg.vec, subject=outliers.sub.vec))
  class(out) <- c('summary.brainGraph_resids', class(out))
  return(out)
}

#' @aliases summary.brainGraph_resids
#' @method print summary.brainGraph_resids
#' @export

print.summary.brainGraph_resids <- function(x, ...) {
  resids <- NULL
  print_title_summary('Structural covariance residuals')
  cat('\n')
  width <- getOption('width')
  dashes <- rep.int('-', width / 4)
  message('# of outliers per region: (sorted in descending order)\n', dashes)
  print(sort(x$outliers$Region, decreasing=TRUE))
  cat('\n')

  message('Number of times each subject was an outlier: (sorted in descending order)\n', dashes)
  print(sort(x$outliers$subject, decreasing=TRUE))
  cat('\n')

  message('Outliers\n', dashes)
  print(x$outliers$DT[order(resids)])
  invisible(x)
}

#' Plot model residuals for each brain region
#'
#' The \code{plot} method lets you check the model residuals for each brain
#' region in a structural covariance analysis. It shows a \emph{qqplot} of the
#' studentized residuals, as output from \code{\link{get.resid}}.
#'
#' @param cols Logical indicating whether to color by group. Default:
#'   \code{FALSE}
#' @inheritParams plot.bg_GLM
#' @export
#' @importFrom ggrepel geom_text_repel
#'
#' @return The \code{plot} method returns a \emph{list} of
#'   \code{\link[ggplot2]{ggplot}} objects
#'
#' @rdname residuals
#' @examples
#' \dontrun{
#' myresids <- get.resids(lhrh, covars)
#' residPlots <- plot(myresids, cols=TRUE)
#'
#' ## Save as a multi-page PDF
#' ml <- marrangeGrob(residPlots, nrow=3, ncol=3)
#' ggsave('residuals.pdf', ml)
#' }

plot.brainGraph_resids <- function(x, region=NULL, cols=FALSE, ids=TRUE, ...) {
  Region <- ind <- mark <- resids <- NULL
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  DT <- summary(x, region=region)$DT.sum
  regions <- DT[, levels(Region)]
  setkeyv(DT, c(gID, sID))

  rnames <- if (isTRUE(ids)) case.names(x) else as.character(seq_len(nobs(x)))
  DT[, ind := rnames, by=Region]
  setkey(DT, Region, resids)
  DT[, x := qnorm(ppoints(resids)), by=Region]
  DT[, mark := ifelse(abs(resids) > abs(mean(resids) + 2*sd(resids)), 1, 0), by=Region]
  DT[mark == 0, ind := '']
  DT[, mark := as.factor(mark)]

  # Local function to plot for a single region
  plot_single <- function(DT.resids, cols) {
    Region <- resids <- NULL
    p <- ggplot(DT.resids, aes(x=x, y=resids, col=get(gID))) +
      geom_text_repel(aes(label=ind), size=3) +
      geom_point(aes(shape=mark, size=mark)) +
      geom_line(aes(x=x, y=x), col='gray50') +
      scale_shape_manual(values=c(20, 17)) +
      scale_size_manual(values=c(2, 3)) +
      theme(plot.title=element_text(hjust=0.5, size=10, face='bold'),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title=paste0('Normal Q-Q: ', DT.resids[, unique(Region)]),
           x='Theoretical Quantiles', y='Sample Quantiles')
    if (isFALSE(cols)) {
      p <- p + scale_color_manual(values=rep.int('black', DT.resids[, length(unique(get(gID)))]))
    }
    return(p)
  }

  p.all <- setNames(vector('list', length(regions)), regions)
  for (z in regions) {
    p.all[[z]] <- plot_single(DT[Region == z], cols)
  }
  return(p.all)
}

#' @export
#' @rdname residuals
nobs.brainGraph_resids <- function(object, ...) dim(object$resids.all)[1L]

#' @export
#' @method case.names brainGraph_resids
#' @rdname residuals
case.names.brainGraph_resids <- function(object, ...) {
  sID <- getOption('bg.subject_id')
  object$resids.all[, as.character(get(sID))]
}
