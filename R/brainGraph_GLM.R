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
#' matrix} associated with \code{lm} objects, and is created from the input
#' \code{data.table} and arguments passed to
#' \code{\link{brainGraph_GLM_design}}. The first column
#' \emph{must} be named \emph{Study.ID}, and all graphs must have a \emph{name}
#' graph-level attribute. The covariates table must be supplied even if you
#' provide your own design matrix \code{X}.
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
#' the same as that in FSL's \emph{randomise}.
#'
#' @param g.list A \code{brainGraphList} object
#' @param covars A \code{data.table} of covariates
#' @param measure Character string of the graph measure of interest
#' @param contrasts Numeric matrix specifying the contrast(s) of interest; if
#'   only one contrast is desired, you can supply a vector
#' @param con.type Character string; either \code{'t'} or \code{'f'} (for t or
#'   F-statistics). Default: \code{'t'}
#' @param outcome Character string specifying the name of the outcome variable,
#'   if it differs from the graph metric (\code{measure})
#' @param X Numeric matrix, if you wish to supply your own design matrix
#' @param con.name Character vector of the contrast name(s); if \code{contrasts}
#'   has row names, those will be used for reporting results
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
#' matrix into covariates of interest and nuisance. Default: \code{beckmann}
#' @param N Integer; number of permutations to create. Default: \code{5e3}
#' @param perms Matrix of permutations, if you would like to provide your own.
#'   Default: \code{NULL}
#' @param long Logical indicating whether or not to return all permutation
#'   results. Default: \code{FALSE}
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM_design}}
#' @export
#' @importFrom permute shuffleSet
#'
#' @return An object of class \code{bg_GLM} containing some input-specific
#'   variables, in addition to:
#'   \item{X}{A named numeric matrix or a 3D array of the design matrix.
#'     Rownames are Study IDs, column names are predictor variables, and
#'     dimnames along the 3rd dimension are region names (if applicable). This
#'     is a 3D array only if \code{outcome != measure} and \code{level ==
#'     'vertex'}.}
#'   \item{y}{A named numeric matrix of the outcome variable. Rownames are Study
#'     IDs and column names are regions. There will be multiple columns only if
#'     \code{outcome == measure} and \code{level == 'vertex'}.}
#'   \item{DT}{A data table with an entry for each vertex (region) containing
#'     statistics of interest}
#'   \item{removed.subs}{A character vector of Study.ID's removed due to
#'     incomplete data (if any)}
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
#' @references Beckmann, C.F. and Jenkinson, M. and Smith, S.M. (2001) General
#'   multi-level linear mdoelling for group analysis in FMRI. Tech Rep.
#'   University of Oxford, Oxford.
#' @references Guttman, I. (1982) \emph{Linear Models: An Introduction}. Wiley,
#'   New York.
#' @references Ridgway, G.R. (2009) Statistical analysis for longitudinal MR
#'   imaging of dementia. PhD thesis.
#' @references Freedman, D. and Lane, D. (1983) A nonstochastic interpretation
#'   of reported significance levels. \emph{J Bus Econ Stat}, \bold{1(4)},
#'   292--298. \url{https://dx.doi.org/10.1080/07350015.1983.10509354}
#' @references Smith, S.M. and Jenkinson, M. and Beckmann, C. and Miller, K. and
#'   Woolrich, M. (2007) Meaningful design and contrast estimability in fMRI.
#'   \emph{NeuroImage}. \bold{34(1)}, 127--36.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2006.09.019}
#' @references Nichols, T.E. and Holmes, A.P. (2001) Nonparametric permutation
#'   tests for functional neuroimaging: A primer with examples. \emph{Human
#'   Brain Mapping}. \bold{15(1)}, 1--25.
#'   \url{https://dx.doi.org/10.1002/hbm.1058}
#' @references Winkler, A.M. and Ridgway, G.R. and Webster, M.A. and Smith, S.M.
#'   and Nichols, T.E. (2014) Permutation inference for the general linear
#'   model. \emph{NeuroImage}. \bold{92}, 381--397.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2014.01.060}
#'
#' @examples
#' \dontrun{
#' conmat <- matrix(c(0, 0, 0, 1), nrow=1)
#' rownames(conmat) <- 'Control > Patient'
#'
#' g.lm <- brainGraph_GLM(g[[6]], covars=covars.all[tract == 1],
#'   measure='strength', contrasts=conmat, alt='greater', permute=TRUE, long=TRUE)
#' }
brainGraph_GLM <- function(g.list, covars, measure, contrasts, con.type=c('t', 'f'),
                           outcome=NULL, X=NULL, con.name=NULL,
                           alternative=c('two.sided', 'less', 'greater'),
                           alpha=0.05, level=c('vertex', 'graph'),
                           permute=FALSE, perm.method=c('freedmanLane', 'terBraak', 'smith'),
                           part.method=c('beckmann', 'guttman', 'ridgway'),
                           N=5e3, perms=NULL, long=FALSE, ...) {
  region <- Outcome <- p.fdr <- p <- Contrast <- i <-
    stat <- p.perm <- perm <- contrast <- V1 <- NULL

  if (!inherits(g.list, 'brainGraphList')) try(g.list <- as_brainGraphList(g.list))
  g.list <- g.list[]

  # Get the outcome variable(s) into a data.table
  sID <- getOption('bg.subject_id')
  level <- match.arg(level)
  DT.y <- glm_data_table(g.list, level, measure)
  setkeyv(DT.y, sID)
  DT.y.m <- melt(DT.y, id.vars=sID, variable.name='region', value.name=measure)

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- match.arg(alternative)
  if (ctype == 'f') alt <- 'two.sided'
  glmSetup <- setup_glm(covars, X, contrasts, ctype, con.name, measure, outcome, DT.y.m, level, ...)
  X <- glmSetup$X; contrasts <- glmSetup$contrasts; con.name <- glmSetup$con.name; DT.y.m <- glmSetup$DT.y.m
  if (is.null(outcome)) outcome <- measure
  regions <- DT.y.m[, levels(region)]
  y <- matrix(DT.y.m[region %in% regions, get(outcome)], ncol=length(regions),
              dimnames=list(DT.y.m[, unique(get(sID))], regions))
  if (level != 'vertex' || outcome != measure) {
    y <- y[, 1, drop=FALSE]
    dimnames(y)[[2]] <- outcome
  }

  #-------------------------------------
  # Do the model fitting/estimation
  #-------------------------------------
  # Handle the case when there is a different design matrix for each region
  dimX <- dim(X)
  if (length(dimX) == 3 && level == 'vertex') {
    # If any columns are all 0, do not fit a model for that region
    runX <- names(which(apply(X, 3, function(p) !any(apply(p, 2, function(n) all(n == 0))))))
    if (length(runX) == 0) {
      stop('At least one covariate is equal to 0 for all regions; please check your data')
    }

    DT.lm <- setNames(vector('list', length(runX)), runX)
    for (k in runX) {
      DT.lm[[k]] <- glm_fit_helper(DT.y.m[region == k], X[, , k],
                                   ctype, contrasts, alt, outcome, 'region', alpha)
    }
    DT.lm <- rbindlist(DT.lm)
  } else {
    DT.lm <- glm_fit_helper(DT.y.m, X, ctype, contrasts, alt, outcome, 'region', alpha)
  }
  DT.lm[, Outcome := outcome]
  DT.lm[, p.fdr :=  p.adjust(p, 'fdr'), by=contrast]
  for (i in seq_along(con.name)) DT.lm[contrast == i, Contrast := con.name[i]]
  setkey(DT.lm, contrast, region)

  out <- list(level=level, covars=glmSetup$covars, X=X, y=y, outcome=outcome, measure=measure, con.type=ctype, contrasts=contrasts,
              con.name=con.name, alt=alt, alpha=alpha, DT=DT.lm, removed.subs=glmSetup$incomp, permute=permute)
  if ((outcome != measure) && level == 'vertex') out$DT.X.m <- glmSetup$DT.X.m
  out$atlas <- guess_atlas(g.list[[1]])
  class(out) <- c('bg_GLM', class(out))
  if (!isTRUE(permute)) return(out)

  #-------------------------------------
  # Permutation testing
  #-------------------------------------
  null.dist <- null.thresh <- NA
  perm.method <- match.arg(perm.method)
  part.method <- match.arg(part.method)
  if (is.null(perms) || dim(perms)[2L] != dimX[1L]) perms <- shuffleSet(n=dimX[1L], nset=N)

  myMax <- maxfun(alt)
  eqn <- if (ctype == 't') 'myMax(gamma / se)' else 'myMax(numer / (se / dfR))'
  dfR <- dimX[1L] - dimX[2L]
  runY <- DT.y.m[, which(!all(get(outcome) == 0)), by=region][, as.character(region)]
  # Different design matrix for each region
  if (length(dimX) == 3 && level == 'vertex') {
    null.dist <- setNames(vector('list', length(runX)), runX)
    for (k in intersect(runX, runY)) {
      null.dist[[k]] <- randomise(perm.method, part.method, N, perms, contrasts, ctype, glmSetup$nC, skip=NULL,
                                  DT.y.m[region == k], outcome, X[, , k], 'region')
    }
    null.dist <- rbindlist(null.dist)
  } else {
    null.dist <- randomise(perm.method, part.method, N, perms, contrasts, ctype, glmSetup$nC, skip=NULL,
                           DT.y.m[region %in% runY], outcome, X, 'region')
  }
  null.dist <- null.dist[, eval(parse(text=eqn)), by=list(perm, contrast)][, !'perm']
  mySort <- sortfun(alt)
  null.thresh <- null.dist[, mySort(V1)[floor((1 - alpha) * N) + 1], by=contrast]
  compfun <- switch(alt,
                    two.sided=function(x, y) sum(abs(x) >= abs(y), na.rm=TRUE),
                    less=function(x, y) sum(x <= y, na.rm=TRUE),
                    greater=function(x, y) sum(x >= y, na.rm=TRUE))
  for (i in seq_along(con.name)) {
    DT.lm[list(i), p.perm := (compfun(null.dist[contrast == i, V1], stat) + 1) / (N + 1), by=region]
  }

  perm <- list(thresh=null.thresh)
  if (isTRUE(long)) perm <- c(perm, list(null.dist=null.dist))
  out <- c(out, list(perm.method=perm.method, part.method=part.method, N=N, perm=perm))
  class(out) <- c('bg_GLM', class(out))
  return(out)
}

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Helper function to set-up for GLM analyses
#'
#' \code{setup_glm} is used to setup the data/objects for any function that uses
#' the main GLM functionality in \code{brainGraph}.
#'
#' This function: removes unused levels from \code{covars} and \code{DT.y.m},
#' removes subjects with incomplete data, creates a design matrix (if not
#' supplied), and supplies names to the contrast matrix.
#'
#' @param DT.y.m A \code{data.table} containing the outcome variable and
#'   \code{Study.ID}
#' @inheritParams GLM
#' @keywords internal
#' @name GLM helpers
#' @rdname glm_helpers

setup_glm <- function(covars, X, contrasts, con.type, con.name, measure, outcome, DT.y.m, level, ...) {
  region <- NULL
  sID <- getOption('bg.subject_id')
  covars <- droplevels(covars)
  if (!sID %in% names(covars)) covars[, eval(sID) := as.character(seq_len(nrow(covars)))]
  incomp <- covars[!complete.cases(covars), get(sID)]
  covars <- covars[!get(sID) %in% incomp]
  if (!is.null(DT.y.m)) DT.y.m <- DT.y.m[!get(sID) %in% incomp]   # Not called for NBS
  setkeyv(covars, sID)

  # Swap the outcome and measure variables, if outcome is not a network metric
  if (!is.null(outcome)) {
    DT.y.m[, eval(outcome) := covars[, get(outcome)], by=region]
    covars[, eval(outcome) := NULL]

    # Graph-level case is simple; only 1 design matrix
    if (level == 'graph') {
      covars[, eval(measure) := DT.y.m[, get(measure)]]

    # Vertex-level has 1 design matrix per region, with the measure changing for each
    } else if (level == 'vertex') {
      DT.X.m <- merge(DT.y.m, covars, by=sID)
      setcolorder(DT.X.m, c(sID, 'region', names(covars[, !get(sID)]), measure))
      DT.X.m[, eval(outcome) := NULL]

      # Get all design matrices into a 3-D array
      X <- setNames(vector('list', DT.X.m[, nlevels(region)]), DT.X.m[, levels(region)])
      for (rgn in DT.X.m[, levels(region)]) {
        X[[rgn]] <- brainGraph_GLM_design(DT.X.m[region == rgn, !'region', with=FALSE], ...)
      }
      X <- abind::abind(X, along=3)
    }
    DT.y.m[, eval(measure) := NULL]
  }
  if (is.null(X)) X <- brainGraph_GLM_design(covars, ...)
  rownames(X) <- covars[, get(sID)]

  tmp <- contrast_names(contrasts, con.type, con.name, X)
  out <- list(covars=covars, incomp=incomp, X=X, contrasts=tmp$contrasts, con.name=tmp$con.name,
              nC=tmp$nC, DT.y.m=DT.y.m)
  if (!is.null(outcome)) if (level == 'vertex') out$DT.X.m <- DT.X.m
  return(out)
}

#' Set contrast and column names for GLM analysis
#'
#' Simple helper function to check the dimensions of contrasts, generate
#' contrast names, and set column names for GLM functions. For F-contrasts, if a
#' \code{matrix} is given, convert it to a \code{list} to simplify processing
#' later.
#' @keywords internal
#' @rdname glm_helpers

contrast_names <- function(contrasts, con.type, con.name, X) {
  if (!is.list(contrasts)) {
    if (is.vector(contrasts)) contrasts <- t(contrasts)
    contrasts <- if (con.type == 't') matrix2list(contrasts) else list(contrasts)
  } else {
    stopifnot(con.type == 'f')
  }

  stopifnot(all(vapply(contrasts, ncol, integer(1)) == dim(X)[2L]))
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
  for (i in seq_len(nC)) dimnames(contrasts[[i]])[[2]] <- dimnames(X)[[2]]

  if (con.type == 't') contrasts <- abind::abind(contrasts, along=1)

  return(list(contrasts=contrasts, con.name=con.name, nC=nC))
}

################################################################################
# MODEL FITTING FUNCTIONS
################################################################################

#' Helper function for GLM fitting
#'
#' @param DT A data.table with all the necessary data; namely \code{region}
#'   (equals \code{graph} if \code{level='graph'}), and the outcome measure(s)
#' @param mykey The \code{key} to key by (to differentiate NBS and other GLM
#'   analyses). For GLM, it is \code{'region'}; for NBS, it is
#'   \code{'Var1,Var2'}.
#' @inheritParams GLM
#' @keywords internal
#' @rdname glm_helpers

glm_fit_helper <- function(DT, X, con.type, contrasts, alt, outcome, mykey, alpha=NULL) {
  ESS <- numer <- stat <- se <- p <- contrast <- ci.low <- ci.high <- NULL
  dimX <- dim(X)
  dfR <- dimX[1L] - dimX[2L]

  XtX <- solve(crossprod(X))
  if (con.type == 'f') {
    DT.lm <- vector('list', length(contrasts))
    CXtX <- lapply(contrasts, function(x) solve(x %*% tcrossprod(XtX, x)))
    rkC <- unlist(lapply(contrasts, function(x) qr(x)$rank))
    for (i in seq_along(contrasts)) {
      DT.lm[[i]] <- DT[, brainGraph_GLM_fit_f(X, get(outcome), dfR, contrasts[[i]], rkC[i], CXtX[[i]]), by=mykey]
    }
    DT.lm <- rbindlist(DT.lm, idcol='contrast')
    DT.lm[, ESS := numer * rkC]
    DT.lm[, stat := (numer / (se / dfR))]
    DT.lm[, numer := NULL]
    DT.lm[, p := pf(stat, rkC, dfR, lower.tail=FALSE)]
  } else if (con.type == 't') {
    pfun <- switch(alt,
                   two.sided=function(stat, df) 2 * pt(abs(stat), df, lower.tail=FALSE),
                   less=function(stat, df) pt(stat, df),
                   greater=function(stat, df) pt(stat, df, lower.tail=FALSE))

    DT.lm <- vector('list', dim(contrasts)[1L])
    for (i in seq_along(DT.lm)) {
      DT.lm[[i]] <- DT[, brainGraph_GLM_fit_t(X, get(outcome), XtX, contrasts[i, , drop=FALSE]), by=mykey]
    }
    DT.lm <- rbindlist(DT.lm, idcol='contrast')
    DT.lm[, stat := gamma / se]
    DT.lm[, p := pfun(stat, dfR)]
    if (!is.null(alpha)) {
      DT.lm[, ci.low := gamma + qt(alpha / 2, dfR) * se]
      DT.lm[, ci.high := gamma + qt(1 - (alpha / 2), dfR) * se]
    }
  }
  return(DT.lm)
}

#' Fit linear models for t and F contrasts
#'
#' \code{brainGraph_GLM_fit_t} fits a linear model for t-contrasts (i.e.,
#' uni-dimensional contrasts) and returns the contrasts of parameter estimates,
#' standard errors, t-statistics, and P-values. If a contrast matrix is
#' supplied, it will return the above values for each row of the matrix.
#'
#' For speed purposes (if it is called from \code{\link{brainGraph_GLM}} and
#' permutation testing is done), this function does not do argument checking.
#'
#' @param y Numeric vector; the outcome variable
#' @param XtX Numeric matrix
#' @param contrast A numeric matrix (1 row) for a single contrast
#' @inheritParams GLM
#' @importFrom RcppEigen fastLmPure
#'
#' @name GLM fits
#' @rdname glm_fit
#'
#' @return \code{brainGraph_GLM_fit_t} - A list containing:
#'   \item{gamma}{The contrast of parameter estimates}
#'   \item{se}{The standard error}
#'
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

brainGraph_GLM_fit_t <- function(X, y, XtX, contrast) {
  est <- fastLmPure(X, y, method=2)
  b <- est$coefficients
  gamma <- contrast %*% b
  sigma.squared <- est$s^2
  var.covar <- sigma.squared * XtX
  se <- sqrt(diag(contrast %*% tcrossprod(var.covar, contrast)))

  list(gamma=as.numeric(gamma), se=se)
}

#' Fit linear models for F contrasts
#'
#' \code{brainGraph_GLM_fit_f} fits a linear model for F contrasts (i.e.,
#' multi-dimensional contrasts) and returns the \emph{extra sum of squares} due
#' to the full model, the sum of squared errors of the full model, the
#' F-statistic, and associated P-value.
#'
#' @param dfR Integer; residual degrees of freedom
#' @param rkC Integer; rank of the contrast matrix
#' @param CXtX Numeric matrix
#' @importFrom RcppEigen fastLmPure
#'
#' @return \code{brainGraph_GLM_fit_f} - A list containing:
#'   \item{numer}{The extra sum of squares due to the full model divided by the
#'     rank of the contrast matrix}
#'   \item{se}{The sum of squared errors of the full model}
#'   \item{contrast}{The contrast number; defaults to \code{1}}
#'
#' @rdname glm_fit

brainGraph_GLM_fit_f <- function(X, y, dfR, contrast, rkC, CXtX) {
  est <- fastLmPure(X, y, method=2)
  b <- as.matrix(est$coefficients)
  gamma <- contrast %*% b
  SSEF <- as.numeric(crossprod(est$residuals))

  numer <- as.numeric(crossprod(gamma, CXtX) %*% gamma / rkC)
  list(numer=numer, se=SSEF)
}

################################################################################
# S3 METHODS FOR "bg_GLM"
################################################################################

#' @method print bg_GLM
#' @keywords internal

print.bg_GLM <- function(x, ...) {
  cat('\nA bg_GLM object at the', x$level, 'level with model:\n\n', formula(x))
  cat('\n\n')
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
  Outcome <- threshold <- NULL
  object$p.sig <- match.arg(p.sig)
  object$printCon <- contrast
  object$digits <- digits
  object$print.head <- print.head
  DT.sum <- object$DT[get(object$p.sig) < alpha]
  if (object$outcome == object$measure) DT.sum[, Outcome := NULL]
  if ('threshold' %in% names(DT.sum)) DT.sum[, threshold := NULL]

  # Change column order and names for `DT.sum`
  newcols <- c('Contrast', 'region', 'gamma', 'ci.low', 'ci.high', 'se',
               'stat', 'p', 'p.fdr', 'p.perm', 'contrast')
  oldnames <- c('region', 'gamma', 'se', 'stat', 'p', 'p.fdr', 'p.perm', 'ci.low', 'ci.high')
  newnames <- c('Region', 'Estimate', 'Std. error', 't value', 'p-value', 'p-value (FDR)', 'p-value (perm.)')
  if (!isTRUE(object$permute)) {
    oldnames <- oldnames[-7]
    newnames <- newnames[-7]
    newcols <- newcols[-10]
  }
  if (object$con.type == 't') {
    clp <- 100 * (1 - object$alpha)
    newnames <- c(newnames, paste0(clp, '% CI ', c('low', 'high')))
  } else if (object$con.type == 'f') {
    newcols[3] <- oldnames[2] <- 'ESS'
    newnames[c(2, 4)] <- c('Extra Sum Sq.', 'F value')
    newcols <- newcols[-c(4, 5)]
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
#' @keywords internal

print.summary.bg_GLM <- function(x, ...) {
  print_title_summary('brainGraph GLM results')
  cat('Level: ', x$level, '\n')

  print_measure_summary(x)
  print_contrast_type_summary(x)
  print_subs_summary(x)

  if (isTRUE(x$permute)) {
    print_permutation_summary(x)
    cat(paste0(toupper(x$con.type), '-statistic threshold (based on the null distribution):\n'))
    print(x$perm$thresh)
  }

  # Print results for each contrast
  message('\n', 'Statistics', '\n', rep('-', getOption('width') / 4))
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
#' @param ids Logical indicating whether to plot Study ID's for outliers.
#'   Otherwise plots the integer index
#' @export
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid gpar grid.draw grid.newpage textGrob
#' @rdname glm
#'
#' @return The \code{plot} method returns a \emph{list} of
#'   \code{\link[ggplot2]{ggplot}} objects
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
  cl.h <- ind <- level <- mark <- resid <- ymax <- NULL
  stopifnot(inherits(x, 'bg_GLM'))
  if (!requireNamespace('gridExtra', quietly=TRUE)) {
    stop('Must install the "gridExtra" package.')
  } else {
    requireNamespace('gridExtra')
  }
  if (!is.numeric(which) || any(which < 1) || any(which > 6)) stop("'which' must be in 1:6")
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  prows <- 1L + (length(which) > 1L)

  model_formula <- split_string(formula(x))
  mytheme <- theme(plot.title=element_text(hjust=0.5),
                   legend.position='none',
                   axis.text.y=element_text(hjust=0.5, angle=90))

  # Local function to plot for a single region
  plot_single <- function(dat, region) {
    leverage <- resid.std <- cook <- ind <- mark <- fit <- leverage.tr <- NULL

    diagPlots <- vector('list', length=6)
    # 1. Resids vs fitted
    if (show[1L]) {
      diagPlots[[1L]] <- ggplot(dat, aes(x=fit, y=resid)) +
        geom_point(aes(shape=mark)) +
        geom_text_repel(aes(x=fit, y=resid, label=ind), size=3) +
        stat_smooth(method='loess', se=FALSE, span=1, col='red') +
        geom_hline(yintercept=0, lty=3) +
        mytheme + labs(title='Residuals vs Fitted', x='Fitted values', y='Residuals')
    }

    # 2. QQ-plot
    if (show[2L]) {
      dat[order(resid.std), x := qnorm(ppoints(resid.std))]
      diagPlots[[2L]] <- ggplot(dat, aes(x=x, y=resid.std)) +
        geom_text_repel(aes(x=x, y=resid.std, label=ind), size=3) +
        geom_line(aes(x=x, y=x), col='gray50', lty=3) +
        geom_point(aes(shape=mark)) +
        mytheme + labs(title='Normal Q-Q', x='Theoretical Quantiles', y='Sample Quantiles')
    }

    # 3. Scale-Location plot
    if (show[3L]) {
      diagPlots[[3L]] <- ggplot(dat, aes(x=fit, y=sqrt(abs(resid.std)))) +
        geom_point(aes(shape=mark)) +
        geom_text_repel(aes(x=fit, y=sqrt(abs(resid.std)), label=ind), size=3) +
        stat_smooth(method='loess', se=FALSE, col='red') +
        ylim(c(0, NA)) +
        mytheme + labs(title='Scale-Location', x='Fitted values', y=expression(sqrt(Standardized~residuals)))
    }

    # 4. Cook's distance
    if (show[4L]) {
      diagPlots[[4L]] <- ggplot(dat, aes(x=seq_len(dimX[1L]), y=cook)) +
        geom_bar(stat='identity', position='identity') +
        geom_text(aes(y=cook, label=ind), size=3, vjust='outward') +
        mytheme + labs(title='Cook\'s distance', x='Obs. number', y='Cook\'s distance')
    }

    # 5. Residual vs Leverage plot
    if (show[5L]) {
      r.hat <- dat[, range(leverage, na.rm=TRUE)]
      hh <- seq.int(min(r.hat[1L], r.hat[2L] / 100), 1, length.out=101)
      dt.cook <- data.table(hh=rep(hh, 2), level=rep(c(0.5, 1.0), each=101))
      dt.cook[, cl.h := sqrt(level * p * (1 - hh) / hh), by=level]
      xmax <- dat[, round(max(leverage), 2)]
      dt.cook[, ymax := .SD[which(xmax == round(hh, 2)), cl.h], by=level]

      diagPlots[[5L]] <- ggplot(dat, aes(x=leverage, y=resid.std)) +
        geom_point(aes(shape=mark)) +
        geom_text_repel(aes(x=leverage, y=resid.std, label=ind), size=3) +
        stat_smooth(method='loess', se=FALSE, col='red') +
        geom_vline(xintercept=0, lty=3, col='gray50') +
        geom_hline(yintercept=0, lty=3, col='gray50') +
        geom_line(data=dt.cook[level == 0.5], aes(x=hh, y=cl.h), col='red', lty=2) +
        geom_line(data=dt.cook[level == 1.0], aes(x=hh, y=cl.h), col='red', lty=2) +
        geom_line(data=dt.cook[level == 0.5], aes(x=hh, y=-cl.h), col='red', lty=2) +
        geom_line(data=dt.cook[level == 1.0], aes(x=hh, y=-cl.h), col='red', lty=2) +
        xlim(extendrange(c(0, xmax))) +
        ylim(dat[, extendrange(range(resid.std), f=0.09)]) +
        mytheme + labs(title='Residuals vs Leverage', x='Leverage', y='Standardized residuals')
    }

    # 6. Cook's dist vs leverage
    if (show[6L]) {
      mytheme$plot.title$size <- 9
      diagPlots[[6L]] <- ggplot(dat, aes(x=leverage.tr, y=cook)) +
        geom_point(aes(shape=mark)) +
        geom_text_repel(aes(x=leverage.tr, y=cook, label=ind), size=3) +
        geom_abline(slope=seq(0, 3, 0.5), lty=2, col='gray50') +
        xlim(c(0, dat[, max(leverage.tr)])) +
        ylim(c(0, dat[, max(cook) * 1.025])) +
        mytheme +
        labs(title=expression("Cook's dist vs Leverage "*h[ii] / (1 - h[ii])),
             x=expression(Leverage~h[ii]), y='Cook\'s distance')
    }

    top_title <- textGrob(paste0('Outcome: ', outcome, '    Region: ', region),
                          gp=gpar(fontface='bold', cex=1.0))
    p.all <- gridExtra::arrangeGrob(grobs=diagPlots[which], nrow=prows, top=top_title, bottom=model_formula)
    return(p.all)
  }

  dimX <- dim(x$X)
  p <- if (length(dimX) == 3L) qr(x$X[, , 1])$rank else qr(x$X)$rank
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
    p.all <- plot_single(DT[['graph']], region='graph-level')
    grid.newpage()
    grid.draw(p.all)

  } else if (x$level == 'vertex') {
    if (all(x$y == 0)) return(NULL)
    if (is.null(region)) region <- regions
    p.all <- setNames(vector('list', length(region)), region)
    runY <- if (outcome == x$measure) names(which(colSums(x$y[, region, drop=FALSE]) != 0)) else region
    for (z in runY) p.all[[z]] <- plot_single(DT[[z]], z)
    p.all[vapply(p.all, is.null, logical(1))] <- NULL
  }
  return(p.all)
}
