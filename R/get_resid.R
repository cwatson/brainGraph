#' Linear model residuals in structural covariance networks
#'
#' \code{get.resid} runs linear models across brain regions listed in a
#' \code{data.table} (e.g., cortical thickness), adjusting for variables in
#' \code{covars} (e.g. age, sex, etc.), and calculates the
#' \emph{externally Studentized} (or \emph{leave-one-out}) residuals.
#'
#' You can choose to run models for each of your subject groups separately or
#' combined (the default) via the \code{method} argument. You may also choose
#' whether to include the mean, per-hemisphere structural measure in the
#' models. Finally, you can specify variables that are present in \code{covars}
#' which you would like to exclude from the models. Optional arguments can be
#' provided that get passed to \code{\link{brainGraph_GLM_design}}.
#'
#' If you do not explicitly specify the atlas name, then it will be guessed from
#' the size of your data. This could cause problems if you are using a custom
#' atlas, with or without the same number of regions as a dataset in the
#' package.
#'
#' @note It is assumed that \code{dt.vol} was created using
#' \code{\link{import_scn}}. In older versions, there were issues when the Study
#' ID was specified as an integer and was not \dQuote{zero-padded}. This is done
#' automatically by \code{\link{import_scn}}, so if you are using an external
#' program, please be sure that the Study ID column is matched in both
#' \code{dt.vol} and \code{covars}.
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
#'   \item{data}{A data.table with the input volume/thickness/etc. data as well
#'     as the covariates used in creating the design matrix.}
#'   \item{X}{The \emph{design matrix}, if using default arguments. If
#'   \code{use.mean=TRUE} then it will be a \emph{named list} with a separate
#'   matrix for the left and right hemispheres. If \code{method='sep.groups'}, a
#'   nested named list for each group and hemisphere.}
#'   \item{method}{The input argument \code{method}}
#'   \item{use.mean}{The input argument \code{use.mean}}
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
  mean.lh <- mean.rh <- NULL
  gID <- getOption('bg.group')
  stopifnot(hasName(covars, gID))
  grps <- covars[, levels(as.factor(get(gID)))]
  sID <- getOption('bg.subject_id')
  if (!hasName(covars, sID)) covars[, eval(sID) := seq_len(dim(covars)[1L])]
  covars[, eval(sID) := check_sID(get(sID))]
  DT.cov <- merge(covars, dt.vol, by=sID)
  idvars <- names(covars)

  method <- match.arg(method)
  xcols <- setdiff(names(covars), exclude.cov)
  if (method == 'sep.groups') {
    xcols <- setdiff(xcols, gID)
    X <- resids <- setNames(vector('list', length(grps)), grps)
  }
  if (isTRUE(use.mean)) {
    idvars <- c(idvars, 'mean.lh', 'mean.rh')
    if (any(grepl('(\\.|_)L[0-9]*', names(dt.vol)))) {
      if (!dim(dt.vol) %in% c(148L, 162L)) {  # Hack to exclude "destrieux" atlases
        lh <- '(\\.|_)L[0-9]*'
        rh <- '(\\.|_)R[0-9]*'
      } else {
        lh <- '^l'
        rh <- '^r'
      }
    } else {
      lh <- '^l'
      rh <- '^r'
    }
    DT.cov[, mean.lh := rowMeans(.SD), .SDcols=patterns(lh)]
    DT.cov[, mean.rh := rowMeans(.SD), .SDcols=patterns(rh)]

    if (method == 'comb.groups') {
      res <- by_hemi(DT.cov, lh, rh, xcols, ...)
      X <- res$X; resids <- res$resids
    } else {
      for (g in grps) {
        res <- by_hemi(DT.cov[get(gID) == g], lh, rh, xcols, ...)
        X[[g]] <- res$X; resids[[g]] <- res$resids
      }
      resids <- do.call(rbind, resids)
    }

  } else {
    ycols <- setdiff(names(dt.vol), sID)
    if (method == 'comb.groups') {
      Y <- as.matrix(DT.cov[, ycols, with=FALSE])
      dimnames(Y)[[1L]] <- DT.cov[, get(sID)]
      X <- brainGraph_GLM_design(DT.cov[, xcols, with=FALSE], ...)
      resids <- rstudent_mat(X, Y)
    } else {
      for (g in grps) {
        Y <- as.matrix(DT.cov[get(gID) == g, ycols, with=FALSE])
        dimnames(Y)[[1L]] <- DT.cov[get(gID) == g, get(sID)]
        X[[g]] <- brainGraph_GLM_design(DT.cov[get(gID) == g, xcols, with=FALSE], ...)
        resids[[g]] <- rstudent_mat(X[[g]], Y)
      }
      resids <- do.call(rbind, resids)
    }
  }

  # Return data to "wide" format with just the residuals
  resids.all <- as.data.table(resids, keep.rownames=sID)
  resids.all <- merge(DT.cov[, c(sID, gID), with=FALSE], resids.all, by=sID)
  setkeyv(resids.all, c(gID, sID))

  out <- list(data=DT.cov, X=X, method=method, use.mean=use.mean,
              resids.all=resids.all, Group=grps, atlas=NULL)
  out$atlas <- if (is.null(atlas)) guess_atlas(dt.vol) else atlas

  class(out) <- c('brainGraph_resids', class(out))
  return(out)
}

# Do the work when including "mean.{lh,rh}" in the models
by_hemi <- function(DT, lh, rh, xcols, ...) {
  X <- resids <- setNames(vector('list', 2L), c('lh', 'rh'))
  for (hemi in names(X)) {
    Y <- as.matrix(DT[, .SD, .SDcols=patterns(get(hemi))])
    dimnames(Y)[[1L]] <- DT[, get(getOption('bg.subject_id'))]
    X[[hemi]] <- brainGraph_GLM_design(DT[, c(xcols, paste0('mean.', hemi)), with=FALSE], ...)
    resids[[hemi]] <- rstudent_mat(X[[hemi]], Y)
  }
  resids <- do.call(cbind, resids)
  list(X=X, resids=resids)
}

#' Calculate studentized residuals with matrix input
#'
#' @inheritParams fastLmBG
#' @noRd

rstudent_mat <- function(X, Y) {
  fits <- fastLmBG(X, Y)
  res <- fits$residuals

  res2 <- res^2
  var.hat <- t(colSums(res2) - t(res2)) / (fits$df.residual - 1L)
  leverage <- colSums(tcrossprod(fits$cov.unscaled, X) * t(X))
  res / sqrt(var.hat * (1 - leverage))
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
  x$Group <- x$resids.all[, levels(as.factor(get(gID)))]
  setkeyv(x$resids.all, c(gID, sID))
  x$data <- x$data[get(sID) %in% case.names(x)]
  return(x)
}

#' Print a summary of residuals for structural covariance data
#'
#' The \code{summary} method prints the number of outliers per region, and the
#' number of times a given subject was an outlier (i.e., across regions).
#'
#' @param region Character vector of region(s) to focus on; default behavior is
#'   to show summary for all regions
#' @param outlier.thresh Number indicating how many standard deviations
#'   above/below the mean indicate an outlier. Default: \code{2}
#' @export
#' @return \code{\link{summary.brainGraph_resids}} returns a list with two
#'   data tables, one of the residuals, and one of only the outlier regions
#' @rdname residuals

summary.brainGraph_resids <- function(object, region=NULL, outlier.thresh=2, ...) {
  Region <- resids <- mark <- NULL
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  if (is.null(region)) region <- region.names(object)
  DT <- melt(object$resids.all[, c(sID, gID, region), with=FALSE], id.vars=c(sID, gID),
             variable.name='Region', value.name='resids')
  setkey(DT, Region, resids)
  DT[, mark := 0]
  DT[abs(resids) > abs(mean(resids) + outlier.thresh * sd(resids)), mark := 1, by=Region]
  outliers <- DT[mark == 1, !'mark']
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
  oldwidth <- getOption('width')
  dashes <- rep_len('-', oldwidth / 4)
  options(width=max(oldwidth / 2, 80L))

  if (length(x$Region) == 1L) {
    message('# of outliers for region ', x$Region, ': ', appendLF=FALSE)
    cat(unname(x$outliers$Region), '\n')
  } else {
    message('# of outliers per region: (sorted in descending order)\n', dashes)
    print(sort(x$outliers$Region, decreasing=TRUE))
  }
  cat('\n')

  if (length(x$Region) == 1L) {
    message('Subjects that are outliers:\n', dashes)
    print(sort(names(x$outliers$subject)))
  } else {
    message('Number of times each subject (top 10%) was an outlier: (sorted in descending order)\n', dashes)
    n <- length(x$outliers$subject)
    print(sort(x$outliers$subject, decreasing=TRUE)[seq_len(floor(n / 10L))])
  }
  cat('\n')

  message('Outliers\n', dashes)
  print(x$outliers$DT[order(resids)])
  options(width=oldwidth)
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
#'
#' @return The \code{plot} method returns a \code{trellis} object or a list of
#'   \code{ggplot} objects
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

plot.brainGraph_resids <- function(x, region=NULL, outlier.thresh=2, cols=FALSE, ids=TRUE, ...) {
  Region <- ind <- mark <- resids <- NULL
  sID <- getOption('bg.subject_id')
  gID <- getOption('bg.group')
  DT <- summary(x, region=region, outlier.thresh=outlier.thresh)$DT.sum
  regions <- DT[, levels(Region)]
  setkeyv(DT, c(gID, sID))

  # 'base' plotting
  if (!requireNamespace('ggplot2', quietly=TRUE)) {
    xylabs <- paste(c('Theoretical', 'Sample'), 'Quantiles')
    if (isTRUE(cols)) {
      DT[Region %in% regions,
         qqmath(~ resids | Region, xlab=xylabs[1L], ylab=xylabs[2L], col=c('red', 'blue')[get(gID)],
                panel=function(x, ...) {
                  panel.qqmathline(x, ...)
                  panel.qqmath(x, ...)
                })]
    } else {
      DT[Region %in% regions,
         qqmath(~ resids | Region, xlab=xylabs[1L], ylab=xylabs[2L],
                panel=function(x, ...) {
                  panel.qqmathline(x, ...)
                  panel.qqmath(x, ...)
                })]
    }

  # 'ggplot2' plotting
  } else {
    rnames <- if (isTRUE(ids)) case.names(x) else as.character(seq_len(nobs(x)))
    DT[, ind := rnames, by=Region]
    setkey(DT, Region, resids)
    DT[, x := qnorm(ppoints(resids)), by=Region]
    DT[mark == 0, ind := '']
    DT[, mark := as.factor(mark)]

    textfun <- if (!requireNamespace('ggrepel', quietly=TRUE)) ggplot2::geom_text else ggrepel::geom_text_repel
    # Local function to plot for a single region
    plot_single <- function(DT.resids, cols) {
      Region <- resids <- NULL
      p <- ggplot2::ggplot(DT.resids, ggplot2::aes(x=x, y=resids, col=get(gID))) +
        textfun(ggplot2::aes(label=ind), size=3) +
        ggplot2::geom_point(ggplot2::aes(shape=mark, size=mark)) +
        ggplot2::geom_line(ggplot2::aes(x=x, y=x), col='gray50') +
        ggplot2::scale_shape_manual(values=c(20, 17)) +
        ggplot2::scale_size_manual(values=c(2, 3)) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=10, face='bold'),
              legend.position='none',
              axis.text.y=ggplot2::element_text(hjust=0.5, angle=90)) +
        ggplot2::labs(title=paste0('Normal Q-Q: ', DT.resids[, unique(Region)]),
             x='Theoretical Quantiles', y='Sample Quantiles')
      if (isFALSE(cols)) {
        p <- p + ggplot2::scale_color_manual(values=rep.int('black', DT.resids[, length(unique(get(gID)))]))
      }
      return(p)
    }

    p.all <- setNames(vector('list', length(regions)), regions)
    for (z in regions) {
      p.all[[z]] <- plot_single(DT[Region == z], cols)
    }
    return(p.all)
  }
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
