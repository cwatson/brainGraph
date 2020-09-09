#' Fit General Linear Models at each vertex of a graph
#'
#' \code{brainGraph_GLM} specifies and fits a General Linear Model (GLM) at each
#' vertex for a given vertex measure (e.g. \emph{degree}) or at the graph-level
#' (e.g., \emph{global efficiency}). Given a contrast matrix or list of
#' contrast(s), and contrast type (for t- or F-contrast(s), respectively) it
#' will calculate the associated statistic(s) for the given contrast(s).
#'
#' The \code{measure} argument will be the graph- or vertex-level measure of
#' interest. Often, this will serve as the model's \emph{outcome} (or dependent,
#' or response) variable; i.e., the variable typically denoted by \emph{y} in
#' GLMs. In other cases, you may wish to choose some other variable as the
#' outcome; e.g., IQ, age, etc. Then you could test for a direct association
#' between the network measure and outcome of interest, or test for another
#' association while adjusting for the network metric. For these applications,
#' you must provide the variable name via the \code{outcome} argument. This is
#' analogous to \code{-evperdat} in FSL's PALM and to \code{--pvr} in
#' FreeSurfer.
#'
#' @section Design matrix:
#' The GLM's \emph{design matrix} will often be identical to the \emph{model
#' matrix} associated with \code{lm} objects (if \dQuote{dummy} coding, the
#' default, is used) and is created from the input \code{data.table} and
#' arguments passed to \code{\link{brainGraph_GLM_design}}. The first column
#' should have the name of \code{getOption('bg.subject_id')} and its values
#' must match the \emph{name} graph-level attribute of the input graphs. The
#' covariates table must be supplied even if you provide your own design matrix
#' \code{X}.  If \code{level='vertex'} and \code{outcome == measure}, there will
#' be a single design for all regions but a separate model for each region
#' (since the graph measure varies by region). If \code{level='vertex'} and
#' \code{outcome != measure}, there will be a separate design (and, therefore, a
#' separate model) for each region even though the outcome is the same in all
#' models.
#'
#' @section Contrasts and statistics:
#' Either t- or F-contrasts can be calculated (specified by \code{con.type}).
#' Multiple t-contrasts can be specified by passing a multi-row \emph{matrix} to
#' \code{contrasts}. Multiple F-contrasts can be specified by passing a
#' \emph{list} of matrices; all matrices must have the same number of columns.
#' All F-contrasts are necessarily \emph{two-sided}; t-contrasts can be any
#' direction, but only one can be chosen per function call.
#' If you choose \code{con.type="f"}, the calculated effect size is represented
#' by the \code{ESS} (\dQuote{extra sum of squares}), the additional variance
#' explained for by the model parameters of interest (as determined by the
#' contrast matrix). The standard error for F-contrasts is the sum of squared
#' errors of the \emph{full model}.
#'
#' @section Non-parametric permutation tests:
#' You can calculate permutations of the data to build a null distribution of
#' the maximum statistic which corrects for multiple testing. To account for
#' complex designs, the design matrix must be \emph{partitioned} into covariates
#' of interest and nuisance; the default method is the \emph{Beckmann} method.
#' The default permutation strategy is that of Freedman & Lane (1983), and is
#' the same as that in FSL's \emph{randomise}. See \code{\link{randomise}}.
#'
#' @param g.list A \code{brainGraphList} object
#' @param covars A \code{data.table} of covariates
#' @param measure Character string of the graph measure of interest
#' @param contrasts Numeric matrix (for T statistics) or list of matrices (for F
#'   statistics) specifying the contrast(s) of interest; if only one contrast is
#'   desired, you can supply a vector (for T statistics)
#' @param con.type Character string; either \code{'t'} or \code{'f'} (for t or
#'   F-statistics). Default: \code{'t'}
#' @param outcome Character string specifying the name of the outcome variable,
#'   if it differs from the graph metric (\code{measure})
#' @param X Numeric matrix, if you wish to supply your own design matrix.
#'   Ignored if \code{outcome != measure}.
#' @param con.name Character vector of the contrast name(s); if \code{contrasts}
#'   has row/list names, those will be used for reporting results
#' @param alternative Character string, whether to do a two- or one-sided test.
#'   Default: \code{'two.sided'}
#' @param alpha Numeric; the significance level. Default: 0.05
#' @param level Character string; either \code{vertex} (default) or
#'   \code{graph}
#' @param permute Logical indicating whether or not to permute group labels.
#'   Default: \code{FALSE}
#' @param perm.method Character string indicating the permutation method.
#'   Default: \code{'freedmanLane'}
#' @param part.method Character string; the method of partitioning the design
#'   matrix into covariates of interest and nuisance. Default: \code{'beckmann'}
#' @param N Integer; number of permutations to create. Default: \code{5e3}
#' @param perms Matrix of permutations, if you would like to provide your own.
#'   Default: \code{NULL}
#' @param long Logical indicating whether or not to return all permutation
#'   results. Default: \code{FALSE}
#' @param ... Arguments passed to \code{\link{brainGraph_GLM_design}}
#' @export
#' @importFrom permute shuffleSet
#'
#' @return An object of class \code{bg_GLM} containing some input-specific
#'   variables (\code{level}, \code{outcome}, \code{measure}, \code{con.type},
#'   \code{contrasts}, \code{con.name}, \code{alt}, \code{alpha},
#'   \code{permute}, \code{perm.method}, \code{part.method}, \code{N}) in
#'   addition to:
#'   \item{DT.Xy}{A data table from which the design matrices are created and
#'     the outcome variable, for all regions.}
#'   \item{X}{A named numeric matrix or a 3D array of the design matrix.
#'     Rownames are subject IDs, column names are predictor variables, and
#'     dimnames along the 3rd dimension are region names (if applicable). This
#'     is a 3D array only if \code{outcome != measure} and \code{level ==
#'     'vertex'}.}
#'   \item{y}{A named numeric matrix of the outcome variable. Rownames are Study
#'     IDs and column names are regions. There will be multiple columns only if
#'     \code{outcome == measure} and \code{level == 'vertex'}.}
#'   \item{DT}{A data table with an entry for each vertex (region) containing
#'     statistics of interest}
#'   \item{removed.subs}{A named integer vector in which the names are subject
#'     ID's of those removed due to incomplete data (if any). The integers
#'     correspond to the row number in the input \code{covars} table.}
#'   \item{runX}{If \code{outcome != measure} and \code{level == 'vertex'}, this
#'     will be a character vector of the regions for which the design matrix is
#'     invertible. Otherwise, it is \code{NULL}.}
#'   \item{runY}{Character vector of the regions for which the outcome variable
#'     has 0 variability. For example, if \code{level='vertex'} and
#'     \code{measure='degree'}, some regions may be disconnected or have the
#'     same degree for all subjects.}
#'   \item{atlas}{Character string of the atlas used (guessed based on the
#'     vertex count).}
#'   \item{perm}{A list containing: \emph{null.dist} (the null distribution of
#'     maximum statistics), \emph{thresh} (the statistic value corresponding
#'     to the \eqn{100 \times (1 - \alpha)}th\% percentile of the null
#'     distribution)}
#'
#' @name GLM
#' @rdname glm
#' @family GLM functions
#' @family Group analysis functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#'
#' @examples
#' \dontrun{
#' conmat <- matrix(c(0, 0, 0, 1), nrow=1)
#' rownames(conmat) <- 'Control > Patient'
#'
#' res.lm <- brainGraph_GLM(g[[6]], covars=covars.all[tract == 1],
#'   measure='strength', contrasts=conmat, alt='greater', permute=TRUE, long=TRUE)
#' }
brainGraph_GLM <- function(g.list, covars, measure, contrasts, con.type=c('t', 'f'),
                           outcome=NULL, X=NULL, con.name=NULL,
                           alternative=c('two.sided', 'less', 'greater'),
                           alpha=0.05, level=c('vertex', 'graph'),
                           permute=FALSE,
                           perm.method=c('freedmanLane', 'terBraak', 'smith',
                                         'draperStoneman', 'manly', 'stillWhite'),
                           part.method=c('beckmann', 'guttman', 'ridgway'),
                           N=5e3, perms=NULL, long=FALSE, ...) {
  level <- match.arg(level)
  ctype <- match.arg(con.type)
  alt <- if (ctype == 'f') 'two.sided' else match.arg(alternative)

  # Initial GLM setup
  glmSetup <- setup_glm(g.list, level, covars, X, contrasts, ctype, con.name, measure, outcome, ...)
  X <- glmSetup$X; y <- glmSetup$y; DT.Xy <- glmSetup$DT.Xy; outcome <- glmSetup$outcome
  contrasts <- glmSetup$contrasts; nC <- glmSetup$nC

  #-------------------------------------
  # Do the model fitting/estimation
  #-------------------------------------
  dimX <- dim(X)
  namesX <- dimnames(X)
  # Handle the case when there is a different design matrix for each region
  if (length(dimX) == 3L) {
    QR <- qr(X)
    runX <- check_if_singular(QR)  # Don't fit a model for regions w/ a rank-deficient design
    runY <- namesX[[3L]]
    if (length(runX) == 0L) {
      print(head(X[, , 1L]), digits=3L)
      stop('Design matrices are rank-deficient for all regions; please check your data')
    }

    fits <- fastLmBG_3d(X, y, runX, QR[runX])
    names(fits$sigma) <- colnames(fits$coefficients) <- colnames(fits$fitted.values) <- runX
    rownames(fits$fitted.values) <- rownames(fits$residuals) <- namesX[[1L]]
    XtX <- aperm(fits$cov.unscaled, c(3L, 1L, 2L))
    fits$var.covar <- aperm(fits$sigma^2 * XtX, c(2L, 3L, 1L))
    fits$se <- sqrt(apply(fits$var.covar, 3L, diag_sq, fits$rank))
    rownames(fits$coefficients) <- rownames(fits$se) <- namesX[[2L]]

    # If any regions were skipped, add NaN/NA into the correct spot in the stats
    if (length(runX) < dimX[3L]) fits <- add_nans(fits, dimX, namesX, runX)

  # Single design matrix for all regions
  } else {
    runY <- names(which(apply(y, 2L, function(r) diff(range(r))) != 0))
    runX <- dimnames(y)[[2L]]
    fits <- fastLmBG(X, y)
    names(fits$sigma) <- colnames(fits$coefficients) <- colnames(fits$fitted.values) <- runX
    fits$var.covar <- aperm(outer(fits$sigma^2, fits$cov.unscaled), c(2L, 3L, 1L))
    fits$se <- sqrt(apply(fits$var.covar, 3L, diag_sq, fits$rank))
    rownames(fits$coefficients) <- rownames(fits$se) <- namesX[[2L]]
  }
  obs_stats <- if (ctype == 't') fastLmBG_t(fits, contrasts, alt, alpha)
    else fastLmBG_f(fits, contrasts)

  DT.lm <- rbindlist(apply(obs_stats, 3L, as.data.table, keep.rownames='region'), idcol='Contrast')

  out <- list(level=level, X=X, y=y, outcome=outcome, measure=measure, DT.Xy=DT.Xy,
              con.type=ctype, contrasts=contrasts, con.name=glmSetup$con.name, alt=alt, alpha=alpha,
              DT=DT.lm, removed.subs=glmSetup$incomp, permute=permute, runX=runX, runY=runY)
  out$atlas <- guess_atlas(g.list$graphs[[1L]])
  out <- c(out, fits)
  class(out) <- c('bg_GLM', class(out))
  if (isFALSE(permute)) return(out)

  #-------------------------------------
  # Permutation testing
  #-------------------------------------
  perm.method <- match.arg(perm.method)
  part.method <- match.arg(part.method)
  if (is.null(perms) || dim(perms)[2L] != dimX[1L]) perms <- shuffleSet(n=dimX[1L], nset=N)

  n <- dimX[1L]; p <- fits$rank; dfR <- n - p; ny <- length(runX)
  null.dist <- if (length(dimX) == 3L) randomise_3d(perm.method, part.method, N, perms, X, y,
                                                    contrasts, ctype, nC, runX, n, p, ny, dfR)
    else randomise(perm.method, part.method, N, perms, X, y,
                   contrasts, ctype, nC, n=n, p=p, ny=ny, dfR=dfR)

  myMax <- switch(alt, two.sided=colMaxAbs, less=colMin, greater=colMax)
  null.max <- apply(null.dist, 3L, myMax, ny)
  mySort <- sortfun(alt)
  index <- floor((1 - alpha) * N) + 1L
  null.thresh <- apply(null.max, 2L, mySort, index)
  names(null.thresh) <- glmSetup$con.name
  compfun <- switch(alt,
                    two.sided=function(x, y) (rowSums(abs(x) >= abs(y), na.rm=TRUE) + 1L),
                    less=function(x, y) (rowSums(x <= y, na.rm=TRUE) + 1L),
                    greater=function(x, y) (rowSums(x >= y, na.rm=TRUE) + 1L))
  p.perm <- vapply(seq_len(nC), function(x) compfun(t(matrix(null.max[, x], N, ny)), obs_stats[runX, 3L, x]),
                   numeric(ny)) / (N + 1L)
  dim(p.perm) <- c(ny, nC)
  dimnames(p.perm) <- list(runX, glmSetup$con.name)
  DT.lm[!is.nan(p), p.perm := c(p.perm[intersect(runX, runY), ])]

  perm <- list(thresh=null.thresh)
  if (isTRUE(long)) perm <- c(perm, list(null.dist=null.dist, null.max=null.max))
  out <- c(out, list(perm.order=perms, perm.method=perm.method, part.method=part.method, N=N, perm=perm))
  class(out) <- c('bg_GLM', class(out))
  return(out)
}

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Helper functions to set-up for GLM analyses
#'
#' \code{setup_glm} is used to setup the data/objects for any function that uses
#' the main GLM functionality in \code{brainGraph}.
#'
#' This function: removes unused levels from \code{covars} and \code{DT.y.m},
#' removes subjects with incomplete data, creates a design matrix (if not
#' supplied), and supplies names to the contrast matrix.
#'
#' @inheritParams GLM
#' @keywords internal
#' @name GLM helpers
#' @rdname glm_helpers

setup_glm <- function(g.list, level, covars, X, contrasts, con.type, con.name, measure, outcome, ...) {
  region <- NULL
  sID <- getOption('bg.subject_id')

  # Get the outcome variable(s) into a data.table
  if (!is.brainGraphList(g.list)) try(g.list <- as_brainGraphList(g.list))
  g.list <- g.list[]

  DT.y <- glm_data_table(g.list, level, measure)
  DT.y.m <- melt(DT.y, id.vars=sID, variable.name='region', value.name=measure)

  covars <- droplevels(copy(covars))
  if (!hasName(covars, sID)) covars[, eval(sID) := seq_len(dim(covars)[1L])]
  covars[, eval(sID) := check_sID(get(sID))]
  incomp <- covars[!complete.cases(covars), which=TRUE]
  names(incomp) <- covars[incomp, get(sID)]
  DT.Xy <- DT.y.m[covars, on=sID]
  if (length(incomp) > 0L) DT.Xy <- DT.Xy[!get(sID) %in% names(incomp)]

  regions <- DT.y.m[, levels(region)]
  # Swap the outcome and measure variables, if outcome is not a network metric
  if (!is.null(outcome)) {
    setcolorder(DT.Xy, c(sID, 'region', setdiff(names(covars), c(sID, outcome)), measure, outcome))

    # Get all design matrices into a 3-D array
    X <- setNames(vector('list', length(regions)), regions)
    for (rgn in regions) {
      X[[rgn]] <- brainGraph_GLM_design(DT.Xy[region == rgn, !c('region', outcome), with=FALSE], ...)
    }
    attrs <- attributes(X[[rgn]])[-c(1L, 2L)]  # Don't include "dim" and "dimnames"
    X <- drop(abind::abind(X, along=3L))
    attributes(X) <- c(attributes(X), attrs)
  } else {
    if (is.null(X)) X <- brainGraph_GLM_design(DT.Xy[region == regions[1L], !c('region', measure), with=FALSE], ...)
    outcome <- measure
  }

  # Create a matrix of outcome variables
  y <- as.matrix(dcast(DT.Xy, paste(sID, '~ region'), value.var=outcome), rownames=sID)
  y <- y[dimnames(X)[[1L]], , drop=FALSE]  # Since dcast changes the row order if input is out of order
  if (level == 'graph' || outcome != measure) {
    y <- y[, 1L, drop=FALSE]
    if (outcome != measure) dimnames(y)[[2L]] <- outcome
  }

  tmp <- contrast_names(contrasts, con.type, con.name, X)
  out <- list(incomp=incomp, X=X, y=y, contrasts=tmp$contrasts, con.name=tmp$con.name,
              nC=tmp$nC, DT.Xy=DT.Xy, outcome=outcome)

  return(out)
}

#' Set contrast and column names for GLM analysis
#'
#' \code{contrast_names} checks the dimensions of contrasts, generates contrast
#' names, and sets column names for GLM functions. For F-contrasts, if a
#' \code{matrix} is given, it converts it to a \code{list} to simplify processing
#' later.
#' @return \code{contrast_names} -- list containing the contrasts (matrix or
#'   list), contrast names, and number of contrasts
#' @keywords internal
#' @rdname glm_helpers

contrast_names <- function(contrasts, con.type, con.name, X) {
  if (!is.list(contrasts)) {
    if (is.vector(contrasts)) contrasts <- t(contrasts)
    contrasts <- if (con.type == 't') matrix2list(contrasts) else list(contrasts)
  } else {
    stopifnot(con.type == 'f')
  }

  ncols <- vapply(contrasts, ncol, integer(1L))
  if (any(ncols != dim(X)[2L])) {
    df <- data.frame(con.number=seq_along(contrasts), ncols=ncols)
    print(df, row.names=FALSE)
    stop('Design matrix has ', dim(X)[2L], ' columns but contrasts do not match (see above)')
  }
  nC <- length(contrasts)

  if (is.null(con.name)) {
    if (is.null(names(contrasts))) names(contrasts) <- paste('Contrast', seq_len(nC))
    con.name <- names(contrasts)
  } else {
    if (length(con.name) < nC) {
      con.name <- c(con.name, paste('Contrast', seq_len(nC)[-(seq_along(con.name))]))
    }
    names(contrasts) <- con.name
  }
  for (i in seq_len(nC)) dimnames(contrasts[[i]])[[2L]] <- dimnames(X)[[2L]]

  if (con.type == 't') contrasts <- abind::abind(contrasts, along=1L)

  return(list(contrasts=contrasts, con.name=con.name, nC=nC))
}

#' \code{check_if_singular} returns the names of the regions in which the design
#' matrix is of full rank (i.e., non-singular) and therefore is invertible.
#'
#' @param QR List of QR decompositions for each design matrix
#' @keywords internal
#' @rdname glm_helpers

check_if_singular <- function(QR) {
  stopifnot(is.list(QR), all(lapply(QR, is.qr)))
  p <- vapply(QR, function(x) dim(x$qr)[2L], integer(1L))
  rk <- vapply(QR, function(x) x$rank, integer(1L))
  names(which(p == rk))
}

################################################################################
# S3 METHODS FOR "bg_GLM"
################################################################################

#' @export
#' @rdname glm

print.bg_GLM <- function(x, ...) {
  cat('\nA', paste0(x$level, '-level'), '"bg_GLM" object\n')
  print_model_summary(x)
  print_contrast_type_summary(x)
  invisible(x)
}

#' Print a summary from brainGraph_GLM analysis
#'
#' The \code{summary} method prints the results, only for which
#' \eqn{p < \alpha}, where \code{alpha} comes from the \code{bg_GLM} object.
#' \dQuote{Simple} P-values are used by default, but you may change this to the
#' FDR-adjusted or permutation P-values via the function argument \code{p.sig}.
#' You may also choose to subset by \emph{contrast}.
#'
#' @param object,x A \code{bg_GLM} object
#' @param p.sig Character string specifying which P-value to use for displaying
#'   significant results (default: \code{p})
#' @param contrast Integer specifying the contrast to plot/summarize; defaults
#'   to showing results for all contrasts
#' @param digits Integer specifying the number of digits to display for P-values
#' @param print.head Logical indicating whether or not to print only the first
#'   and last 5 rows of the statistics tables (default: \code{TRUE})
#' @export
#' @rdname glm

summary.bg_GLM <- function(object, p.sig=c('p', 'p.fdr', 'p.perm'), contrast=NULL, alpha=object$alpha,
                           digits=max(3L, getOption('digits') - 2L), print.head=TRUE, ...) {
  stopifnot(inherits(object, 'bg_GLM'))
  threshold <- NULL
  object$p.sig <- match.arg(p.sig)
  object$printCon <- contrast
  object$digits <- digits
  object$print.head <- print.head
  DT.sum <- object$DT[get(object$p.sig) < alpha]
  if (hasName(DT.sum, 'threshold')) DT.sum[, threshold := NULL]

  # Change column order and names for `DT.sum`
  newcols <- c('Contrast', 'region', 'gamma', 'ci.low', 'ci.high',
               'se', 'stat', 'p', 'p.fdr', 'p.perm')
  oldnames <- c('region', 'gamma', 'se', 'stat', 'p', 'p.fdr', 'p.perm', 'ci.low', 'ci.high')
  newnames <- c('Region', 'Estimate', 'Std. error', 't value', 'p-value', 'p-value (FDR)', 'p-value (perm.)')
  if (isFALSE(object$permute)) {
    oldnames <- oldnames[-7L]
    newnames <- newnames[-7L]
    newcols <- newcols[-10L]
  }
  if (object$con.type == 't') {
    clp <- 100 * (1 - object$alpha)
    newnames <- c(newnames, paste0(clp, '% CI ', c('low', 'high')))
  } else if (object$con.type == 'f') {
    newcols[3L] <- oldnames[2L] <- 'ESS'
    newnames[c(2L, 4L)] <- c('Extra Sum Sq.', 'F value')
    newcols <- newcols[-c(4L, 5L)]
    oldnames <- oldnames[-grep('ci.', oldnames)]
  }
  setcolorder(DT.sum, newcols)
  setnames(DT.sum, oldnames, newnames)

  object$DT.sum <- DT.sum
  class(object) <- c('summary.bg_GLM', class(object))
  return(object)
}

#' @aliases summary.bg_GLM
#' @method print summary.bg_GLM
#' @export

print.summary.bg_GLM <- function(x, ...) {
  print_title_summary('brainGraph GLM results (', x$level, '-level)')

  print_model_summary(x)
  print_contrast_type_summary(x)
  print_subs_summary(x)

  if (isTRUE(x$permute)) {
    print_permutation_summary(x)
    cat(paste0(toupper(x$con.type), '-statistic threshold (based on the null distribution):\n'))
    print(x$perm$thresh)
  }

  # Print results for each contrast
  message('\n', 'Statistics', '\n', rep.int('-', getOption('width') / 4))
  print_contrast_stats_summary(x)

  invisible(x)
}

#' Plot GLM diagnostics for a brain network
#'
#' The \code{plot} method plots the GLM diagnostics (similar to that of
#' \code{\link[stats]{plot.lm}}). There are a total of 6 possible plots,
#' specified by the \code{which} argument; the behavior is the same as in
#' \code{\link[stats]{plot.lm}}. Please see the help for that function.
#'
#' @param region Character string specifying which region's results to
#'   plot; only relevant if \code{level='vertex'}. Default: \code{NULL}
#' @param which Integer vector indicating which of the 6 plots to print to the
#'   plot device. Default: \code{c(1:3, 5)}
#' @param ids Logical indicating whether to plot subject ID's for outliers.
#'   Otherwise plots the integer index
#' @export
#' @importFrom grid gpar grid.draw grid.newpage textGrob
#' @rdname glm
#'
#' @return The \code{plot} method returns a \emph{list} of \code{ggplot} objects
#'   (if installed) or writes the plots to a PDF in the current directory named
#'   \code{bg_GLM_diagnostics.pdf}
#' @seealso \code{\link[stats]{plot.lm}}
#'
#' @examples
#' \dontrun{
#' ## Save objects and then to multipage PDF
#' lmPlots <- plot(x)
#' ggsave('lmPlots.pdf', lmPlots)
#'
#' ## Save all the GLM sub-objects from MTPC analysis
#' res.mtpc <- mtpc(...)
#' glmPlots <- lapply(res.mtpc$res.glm, plot, which=1:6)
#' ml <- marrangeGrob(glmPlots, nrow=1, ncol=1)
#' ggsave('glmPlots.pdf', ml, width=8.5, height=11)
#' }

plot.bg_GLM <- function(x, region=NULL, which=c(1L:3L, 5L), ids=TRUE, ...) {
  cl.h <- ind <- level <- mark <- resid <- NULL
  stopifnot(inherits(x, 'bg_GLM'))
  if (!is.numeric(which) || any(which < 1L) || any(which > 6L)) stop("'which' must be in 1:6")
  show <- rep.int(FALSE, 6L)
  show[which] <- TRUE
  if (!requireNamespace('ggplot2', quietly=TRUE)) {
    mfrow <- switch(length(which),
                    c(1, 1), c(1, 2), c(1, 3),
                    c(2, 2), c(2, 3), c(2, 3))
    pdf('bg_GLM_diagnostics.pdf')
    par(mfrow=mfrow)
  } else {
    if (!requireNamespace('gridExtra', quietly=TRUE)) stop('Must install the "gridExtra" package.')
    textfun <- if (!requireNamespace('ggrepel', quietly=TRUE)) ggplot2::geom_text else ggrepel::geom_text_repel
    prows <- 1L + (length(which) > 1L)
    mytheme <- ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),
                              legend.position='none',
                              axis.text.y=ggplot2::element_text(hjust=0.5, angle=90))
    model_formula <- split_string(formula(x))
  }

  xlabs <- c('Fitted values', 'Theoretical Quantiles', 'Fitted values', 'Obs. number',
             'Leverage', expression(Leverage~h[ii]))
  ylabs <- c('Residuals', 'Sample Quantiles', expression(sqrt(Standardized~residuals)),
             'Cook\'s distance', 'Standardized residuals', 'Cook\'s distance')
  mains <- c('Residuals vs Fitted', 'Normal Q-Q', 'Scale-Location', 'Cook\'s distance',
             'Residuals vs Leverage', expression("Cook's dist vs Leverage "*h[ii] / (1 - h[ii])))

  # Local function to plot for a single region ('ggplot2')
  #-----------------------------------------------------------------------------
  plot_single <- function(dat, region) {
    leverage <- resid.std <- cook <- ind <- mark <- fit <- leverage.tr <- NULL

    diagPlots <- vector('list', length=6L)
    # 1. Resids vs fitted
    if (show[1L]) {
      diagPlots[[1L]] <- ggplot2::ggplot(dat, ggplot2::aes(x=fit, y=resid)) +
        ggplot2::geom_point(ggplot2::aes(shape=mark)) +
        textfun(ggplot2::aes(x=fit, y=resid, label=ind), size=3) +
        ggplot2::stat_smooth(method='loess', se=FALSE, span=1, col='red') +
        ggplot2::geom_hline(yintercept=0, lty=3) +
        mytheme + ggplot2::labs(title=mains[1L], xlab=xlabs[1L], ylab=ylabs[1L])
    }

    # 2. QQ-plot
    if (show[2L]) {
      dat[order(resid.std), x := qnorm(ppoints(resid.std))]
      diagPlots[[2L]] <- ggplot2::ggplot(dat, ggplot2::aes(x=x, y=resid.std)) +
        textfun(ggplot2::aes(x=x, y=resid.std, label=ind), size=3) +
        ggplot2::geom_line(ggplot2::aes(x=x, y=x), col='gray50', lty=3) +
        ggplot2::geom_point(ggplot2::aes(shape=mark)) +
        mytheme + ggplot2::labs(title=mains[2L], x=xlabs[2L], y=ylabs[2L])
    }

    # 3. Scale-Location plot
    if (show[3L]) {
      diagPlots[[3L]] <- ggplot2::ggplot(dat, ggplot2::aes(x=fit, y=sqrt(abs(resid.std)))) +
        ggplot2::geom_point(ggplot2::aes(shape=mark)) +
        textfun(ggplot2::aes(x=fit, y=sqrt(abs(resid.std)), label=ind), size=3) +
        ggplot2::stat_smooth(method='loess', se=FALSE, col='red') +
        ggplot2::ylim(c(0, NA)) +
        mytheme + ggplot2::labs(title=mains[3L], x=xlabs[3L], y=ylabs[3L])
    }

    # 4. Cook's distance
    if (show[4L]) {
      diagPlots[[4L]] <- ggplot2::ggplot(dat, ggplot2::aes(x=seq_len(dimX[1L]), y=cook)) +
        ggplot2::geom_bar(stat='identity', position='identity') +
        ggplot2::geom_text(ggplot2::aes(y=cook, label=ind), size=3, vjust='outward') +
        mytheme + ggplot2::labs(title=mains[4L], x=xlabs[4L], y=ylabs[4L])
    }

    # 5. Residual vs Leverage plot
    if (show[5L]) {
      r.hat <- dat[, range(leverage, na.rm=TRUE)]
      hh <- seq.int(min(r.hat[1L], r.hat[2L] / 100), 1L, length.out=101L)
      dt.cook <- data.table(hh=rep.int(hh, 2L), level=rep(c(0.5, 1.0), each=101L))
      dt.cook[, cl.h := sqrt(level * p * (1 - hh) / hh), by=level]

      diagPlots[[5L]] <- ggplot2::ggplot(dat, ggplot2::aes(x=leverage, y=resid.std)) +
        ggplot2::geom_point(ggplot2::aes(shape=mark)) +
        textfun(ggplot2::aes(x=leverage, y=resid.std, label=ind), size=3) +
        ggplot2::stat_smooth(method='loess', se=FALSE, col='red') +
        ggplot2::geom_vline(xintercept=0, lty=3, col='gray50') +
        ggplot2::geom_hline(yintercept=0, lty=3, col='gray50') +
        ggplot2::geom_line(data=dt.cook[level == 0.5], ggplot2::aes(x=hh, y=cl.h), col='red', lty=2, na.rm=TRUE) +
        ggplot2::geom_line(data=dt.cook[level == 1.0], ggplot2::aes(x=hh, y=cl.h), col='red', lty=2, na.rm=TRUE) +
        ggplot2::geom_line(data=dt.cook[level == 0.5], ggplot2::aes(x=hh, y=-cl.h), col='red', lty=2, na.rm=TRUE) +
        ggplot2::geom_line(data=dt.cook[level == 1.0], ggplot2::aes(x=hh, y=-cl.h), col='red', lty=2, na.rm=TRUE) +
        ggplot2::xlim(extendrange(c(0, r.hat[2L]))) +
        ggplot2::ylim(dat[, extendrange(range(resid.std), f=0.09)]) +
        mytheme + ggplot2::labs(title=mains[5L], x=xlabs[5L], y=ylabs[5L])
    }

    # 6. Cook's dist vs leverage
    if (show[6L]) {
      mytheme$plot.title$size <- 9
      diagPlots[[6L]] <- ggplot2::ggplot(dat, ggplot2::aes(x=leverage.tr, y=cook)) +
        ggplot2::geom_point(ggplot2::aes(shape=mark)) +
        textfun(ggplot2::aes(x=leverage.tr, y=cook, label=ind), size=3) +
        ggplot2::geom_abline(slope=seq.int(0, 3, 0.5), lty=2, col='gray50') +
        ggplot2::xlim(c(0, dat[, max(leverage.tr)])) +
        ggplot2::ylim(c(0, dat[, max(cook) * 1.025])) +
        mytheme +
        ggplot2::labs(title=mains[6L], x=xlabs[6L], y=ylabs[6L])
    }

    top_title <- textGrob(paste0('Outcome: ', outcome, '    Region: ', region),
                          gp=gpar(fontface='bold', cex=1.0))
    p.all <- gridExtra::arrangeGrob(grobs=diagPlots[which], nrow=prows, top=top_title, bottom=model_formula)
    return(p.all)
  }

  # Local function to plot for a single region ('base')
  #-----------------------------------------------------------------------------
  plot_single_base <- function(dat, region) {
    leverage <- fit <- leverage.tr <- NULL
    label.pos <- c(4, 2)  # Default value in 'plot.lm'
    text.id <- function(x, y, label, adj.x=TRUE) {
      labpos <- if (adj.x) label.pos[1 + as.numeric(x > mean(range(x)))] else 3
      text(x, y, label, cex=0.75, xpd=TRUE, pos=labpos, offset=0.25)
    }

    # 1. Resids vs. fitted
    if (show[1L]) {
      ylim <- dat[, extendrange(resid, f=0.08)]
      dat[mark == 0, plot(fit, resid, pch=19, ylim=ylim,
                          main=mains[1L], xlab=xlabs[1L], ylab=ylabs[1L])]
      dat[mark == 1, points(fit, resid, pch=17)]
      dat[, panel.smooth(fit, resid, iter=3, lwd=2, pch=NA)]
      dat[, abline(h=0, lty=3)]
      dat[mark == 1, text.id(fit, resid, ind)]
    }

    # 2. QQ-plot
    if (show[2L]) {
      ylim <- dat[, range(resid.std, na.rm=TRUE)]
      ylim[2L] <- extendrange(ylim, f=0.075)[2L]
      #dat[order(resid.std), x := qnorm(ppoints(resid.std))]
      qq <- dat[, qqnorm(resid.std, main=mains[2L], ylab=ylabs[2L], ylim=ylim)]
      dat[, qqline(resid.std, lty=3, col='gray50')]
      dat[mark == 1, text.id(x, resid.std, ind)]
    }

    # 3. Scale-Location plot
    if (show[3L]) {
      ylim <- dat[, c(0, max(sqrt(abs(resid.std)), na.rm=TRUE))]
      dat[mark == 0, plot(fit, sqrt(abs(resid.std)), pch=19, ylim=ylim,
                          main=mains[3L], xlab=xlabs[3L], ylab=ylabs[3L])]
      dat[mark == 1, points(fit, sqrt(abs(resid.std)), pch=17)]
      dat[, panel.smooth(fit, sqrt(abs(resid.std)), iter=3, lwd=2, pch=NA)]
      dat[mark == 1, text.id(fit, sqrt(abs(resid.std)), ind)]
    }

    # 4. Cook's distance
    if (show[4L]) {
      ylim <- dat[, c(0, max(cook, na.rm=TRUE) * 1.075)]
      dat[, plot(cook, type='h', ylim=ylim, main=mains[4L], xlab=xlabs[4L], ylab=ylabs[4L])]
      xpos <- dat[mark == 1, which=TRUE]
      dat[mark == 1, text.id(xpos, cook, ind, adj.x=FALSE)]
    }

    # 5. Residual vs Leverage plot
    if (show[5L]) {
      xlim <- dat[, c(0, max(leverage, na.rm=TRUE))]
      ylim <- dat[, extendrange(resid.std, f=0.08)]
      r.hat <- dat[, range(leverage, na.rm=TRUE)]
      hh <- seq.int(min(r.hat[1L], r.hat[2L] / 100), 1L, length.out=101L)
      dt.cook <- data.table(hh=rep.int(hh, 2L), level=rep(c(0.5, 1.0), each=101L))
      dt.cook[, cl.h := sqrt(level * p * (1 - hh) / hh), by=level]

      dat[mark == 0, plot(leverage, resid.std, xlim=xlim, ylim=ylim, main=mains[5L],
                          xlab=xlabs[5L], ylab=ylabs[5L], pch=19)]
      dat[mark == 1, points(leverage, resid.std, pch=17)]
      dat[, panel.smooth(leverage, resid.std, iter=3, lwd=2, pch=NA)]
      abline(h=0, v=0, lty=3, col='gray')
      dt.cook[, lines(hh, cl.h, lty=2, col=2), by=level]
      dt.cook[, lines(hh, -cl.h, lty=2, col=2), by=level]
      dat[mark == 1, text.id(leverage, resid.std, ind)]
    }

    # 6. Cook's dist vs leverage
    if (show[6L]) {
      xlim <- dat[, c(0, max(leverage.tr, na.rm=TRUE))]
      ylim <- dat[, c(0, max(cook, na.rm=TRUE) * 1.025)]
      dat[mark == 0, plot(leverage.tr, cook, xlim=xlim, ylim=ylim, main=mains[6L],
                          xlab=xlabs[6L], ylab=ylabs[6L], pch=19)]
      dat[mark == 1, points(leverage.tr, cook, pch=17)]
      dat[, panel.smooth(leverage.tr, cook, iter=3, lwd=2, pch=NA)]
      dat[mark == 1, text.id(leverage.tr, cook, ind)]
      bval <- dat[, pretty(sqrt(p * cook / leverage.tr), 5)]
      bi2 <- bval^2
      usr <- par("usr")
      xmax <- usr[2L]
      ymax <- usr[4L]
      for (i in seq_along(bval)) {
        if (p * ymax > bi2[i] * xmax) {
          xi <- xmax + strwidth(" ") / 3
          yi <- bi2[i] * xi / p
          abline(0, bi2[i], lty=2)
          text(xi, yi, paste(bval[i]), adj=0, xpd=TRUE)
        } else {
          yi <- ymax - 1.5 * strheight(" ")
          xi <- p * yi / bi2[i]
          lines(c(0, xi), c(0, yi), lty=2)
          text(xi, ymax - 0.8 * strheight(" "), paste(bval[i]), adj=0.5, xpd=TRUE)
        }
      }
    }

    title(paste0('Outcome: ', outcome, '    Region: ', region), line=-1, outer=TRUE)
  }

  dimX <- dim(x$X)
  p <- if (length(dimX) == 3L) qr.default(x$X[, , 1L])$rank else qr.default(x$X)$rank
  rnames <- if (isTRUE(ids)) case.names(x) else as.character(seq_len(dimX[1L]))
  outcome <- x$outcome

  # Get all GLM diagnostics
  fits <- fitted(x)
  resids <- residuals(x)
  hat <- hatvalues(x)
  hat.tr <- hat / (1 - hat)
  resid.std <- rstandard(x)
  cook <- cooks.distance(x)
  regions <- region.names(x)
  DT <- setNames(vector('list', length(regions)), regions)
  for (i in regions) {
    DT[[i]] <- data.table(fit=fits[, i], resid=resids[, i], leverage=hat[, i],
                          resid.std=resid.std[, i], cook=cook[, i], leverage.tr=hat.tr[, i])
    DT[[i]][, ind := rnames]
    DT[[i]][, mark := ifelse(abs(resid) < mean(resid) + 2 * sd(resid), 0, 1)]
    DT[[i]][mark == 0, ind := '']
    DT[[i]][, mark := as.factor(mark)]
  }

  if (x$level == 'graph') {
    if (!requireNamespace('ggplot2', quietly=TRUE)) {
      plot_single_base(DT[['graph']], region='graph-level')
      return(invisible(x))
    } else {
      p.all <- plot_single(DT[['graph']], region='graph-level')
      grid.newpage()
      grid.draw(p.all)
    }

  } else if (x$level == 'vertex') {
    if (length(x$runY) == 0) return(NULL)
    run <- if (is.null(region)) run <- intersect(intersect(x$runX, x$runY), regions) else region

    # 'base' plotting
    if (!requireNamespace('ggplot2', quietly=TRUE)) {
      for (z in run) plot_single_base(DT[[z]], z)
      return(invisible(x))

    # 'ggplot2' plotting
    } else {
      p.all <- setNames(vector('list', length(run)), run)
      for (z in run) p.all[[z]] <- plot_single(DT[[z]], z)
      p.all[vapply(p.all, is.null, logical(1L))] <- NULL
    }
  }
  return(p.all)
}
