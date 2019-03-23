#' Fit General Linear Models at each vertex of a graph
#'
#' \code{brainGraph_GLM} specifies and fits a General Linear Model (GLM) at each
#' vertex for a given vertex measure (e.g. \emph{degree}) or at the graph-level
#' (e.g., \emph{global efficiency}). Given a contrast matrix and type (i.e., t-
#' or F-contrast), it will calculate the associated statistic(s).
#'
#' The \code{measure} argument will be the graph- or vertex-level measure of
#' interest. Often, this will serve as the model's \emph{outcome} (or dependent,
#' or response) variable; i.e., the variable typically denoted by \emph{y} in
#' GLMs. In other cases, you may wish to choose some other variable as the
#' outcome; e.g., IQ, age, etc. Then you could test for a direct association
#' between the network measure and outcome of interest, or test for another
#' association while adjusting for the network metric. For these applications,
#' you must provide the variable name via the \code{outcome} argument.
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
#' @section Statistics:
#' Both t- and F-contrasts are allowed by supplying a \emph{matrix} to the
#' argument \code{con.mat}. If you supply a multi-row matrix and you choose
#' \code{con.type="t"}, then statistics will be calculated for each contrast
#' individually. If you choose \code{con.type="f"}, in the return object's data
#' table the calculated effect size is represented by the \code{ESS}
#' (\dQuote{extra sum of squares}), the additional variance explained for by the
#' model parameters of interest (as determined by the contrast matrix). The
#' standard error for F-contrasts is the sum of squared errors of the \emph{full
#' model}.
#'
#' @section Non-parametric permutation tests:
#' You can calculate permutations of the data to build a null distribution of
#' the maximum statistic which corrects for multiple testing. To account for
#' complex designs, the design matrix must be \emph{partitioned} into covariates
#' of interest and nuisance; the default method is the \emph{Beckmann} method.
#' The default permutation strategy is that of Freedman & Lane (1983), and is
#' the same as that in FSL's \emph{randomise}.
#'
#' @param g.list A list of \code{igraph} graph objects for all subjects
#' @param covars A \code{data.table} of covariates
#' @param measure Character string of the graph measure of interest
#' @param con.mat Numeric matrix specifying the contrast(s) of interest; if
#'   only one contrast is desired, you can supply a vector
#' @param con.type Character string; either \code{'t'} or \code{'f'} (for t or
#'   F-statistics). Default: \code{'t'}
#' @param outcome Character string specifying the name of the outcome variable,
#'   if it differs from the graph metric (\code{measure})
#' @param X Numeric matrix, if you wish to supply your own design matrix
#' @param con.name Character vector of the contrast name(s); if \code{con.mat}
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
#'   \item{X}{A numeric matrix; a copy of the \emph{design matrix}}
#'   \item{y}{A numeric vector or matrix of the outcome variable}
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
#' @aliases brainGraph_GLM
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
#' ## Note that I concatenate the graphs from each group's 6th threshold
#' g.lm <- brainGraph_GLM(g.list=do.call(Map, c(c, g))[[6]],
#'   covars=covars.all[tract == 1],
#'   measure='strength', con.mat=conmat, alt='greater',
#'   permute=TRUE, long=TRUE)
#' }
brainGraph_GLM <- function(g.list, covars, measure, con.mat, con.type=c('t', 'f'),
                           outcome=NULL, X=NULL, con.name=NULL,
                           alternative=c('two.sided', 'less', 'greater'),
                           alpha=0.05, level=c('vertex', 'graph'),
                           permute=FALSE, perm.method=c('freedmanLane', 'terBraak', 'smith'),
                           part.method=c('beckmann', 'guttman', 'ridgway'),
                           N=5e3, perms=NULL, long=FALSE, ...) {
  Study.ID <- region <- Outcome <- p.fdr <- p <- Contrast <- i <-
    stat <- p.perm <- perm <- contrast <- V1 <- NULL

  # Get the outcome variable(s) into a data.table
  DT.y <- data.table(Study.ID=vapply(g.list, graph_attr, character(1), 'name'))
  level <- match.arg(level)
  if (level == 'vertex') {
    y <- t(vapply(g.list, vertex_attr, numeric(vcount(g.list[[1]])), measure))
    DT.y <- cbind(DT.y, y)
    setnames(DT.y, 2:ncol(DT.y), V(g.list[[1]])$name)
  } else if (level == 'graph') {
    DT.y[, graph := vapply(g.list, graph_attr, numeric(1), measure)]
  }
  setkey(DT.y, Study.ID)
  DT.y.m <- melt(DT.y, id.vars='Study.ID', variable.name='region', value.name=measure)

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- match.arg(alternative)
  if (ctype == 'f') alt <- 'two.sided'
  glmSetup <- setup_glm(covars, X, con.mat, ctype, con.name, measure, outcome, DT.y.m, level, ...)
  covars <- glmSetup$covars; X <- glmSetup$X; incomp <- glmSetup$incomp
  con.mat <- glmSetup$con.mat; con.name <- glmSetup$con.name; nC <- glmSetup$nC
  DT.y.m <- glmSetup$DT.y.m; outcome <- glmSetup$outcome
  if (level == 'vertex') {
    y <- as.matrix(DT.y[!Study.ID %in% incomp, -1])
  } else if (level == 'graph') {
    y <- DT.y.m[, get(outcome)]
  }

  #-------------------------------------
  # Do the model fitting/estimation
  #-------------------------------------
  # Handle the case when there is a different design matrix for each region
  if (length(dim(X)) == 3 && level == 'vertex') {
    regions <- dimnames(X)[[3]]
    run <- rep(1, dim(X)[3])
    for (k in seq_along(run)) {
      if (any(colSums(X[, , k]) == 0)) run[k] <- 0
    }
    run <- which(run == 1)
    DT.lm <- setNames(vector('list', length(run)), regions[run])
    for (k in run) {
      DT.lm[[regions[k]]] <- glm_fit_helper(DT.y.m[region == regions[k]], X[, , k],
                                            ctype, con.mat, alt, outcome, 'region', alpha)
    }
    DT.lm <- rbindlist(DT.lm)
  } else {
    DT.lm <- glm_fit_helper(DT.y.m, X, ctype, con.mat, alt, outcome, 'region', alpha)
  }
  DT.lm[, Outcome := outcome]
  DT.lm[, p.fdr :=  p.adjust(p, 'fdr'), by=contrast]
  for (i in seq_along(con.name)) DT.lm[contrast == i, Contrast := con.name[i]]
  setkey(DT.lm, contrast, region)

  #-------------------------------------
  # Permutation testing
  #-------------------------------------
  null.dist <- null.thresh <- NA
  if (isTRUE(permute)) {
    perm.method <- match.arg(perm.method)
    part.method <- match.arg(part.method)
    if (is.null(perms) || ncol(perms) != nrow(X)) {
      perms <- shuffleSet(n=nrow(X), nset=N)
    }

    zeroregions <- DT.y.m[, .SD[all(get(outcome) == 0)], by=region][, unique(as.character(region))]
    # Different design matrix for each region
    if (length(dim(X)) == 3 && level == 'vertex') {
      maxfun <- switch(alt,
                       two.sided=function(x) max(abs(x), na.rm=TRUE),
                       less=function(x) min(x, na.rm=TRUE),
                       greater=function(x) max(x, na.rm=TRUE))
      null.dist <- setNames(vector('list', length(run)), regions[run])
      for (k in run) {
        if (!regions[k] %in% zeroregions) {
          null.dist[[regions[k]]] <- randomise(perm.method, part.method, ctype, N, perms, DT.y.m[region == regions[k]],
                                               nC, outcome, X[, , k], con.mat, alt)
        }
      }
      null.dist <- rbindlist(null.dist, idcol='region')
      null.dist <- cbind(null.dist, data.table(perm=rep(seq_len(N), times=length(run))))
      if (ctype == 't') {
        null.dist <- null.dist[, maxfun(V1), by=list(perm, contrast)][, !'perm']
      } else if (ctype == 'f') {
        null.dist <- null.dist[, max(V1), by=list(perm, contrast)][, !'perm']
      }
    } else {
      null.dist <- randomise(perm.method, part.method, ctype, N, perms, DT.y.m[!region %in% zeroregions], nC, outcome, X, con.mat, alt)
    }
    sortfun <- switch(alt,
                      two.sided=function(x) sort(abs(x)),
                      less=function(x) sort(x, decreasing=TRUE),
                      greater=function(x) sort(x))
    null.thresh <- null.dist[, sortfun(V1)[floor((1 - alpha) * N) + 1], by=contrast]
    compfun <- switch(alt,
                      two.sided=function(x, y) sum(abs(x) >= abs(y), na.rm=TRUE),
                      less=function(x, y) sum(x <= y, na.rm=TRUE),
                      greater=function(x, y) sum(x >= y, na.rm=TRUE))
    for (i in seq_along(con.name)) {
      DT.lm[list(i), p.perm := (compfun(null.dist[contrast == i, V1], stat) + 1) / (N + 1), by=region]
    }
  }

  perm <- list(thresh=null.thresh)
  if (isTRUE(long)) perm <- c(perm, list(null.dist=null.dist))
  out <- list(level=level, X=X, y=y, outcome=outcome, measure=measure, con.type=ctype, con.mat=con.mat,
              con.name=con.name, alt=alt, alpha=alpha, DT=DT.lm, removed.subs=incomp,
              permute=permute, perm.method=perm.method, part.method=part.method, N=N, perm=perm)
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
#' @name GLMhelpers
#' @aliases setup_glm
#' @rdname glm_helpers

setup_glm <- function(covars, X, con.mat, con.type, con.name, measure, outcome, DT.y.m, level, ...) {
  Study.ID <- region <- NULL
  covars <- droplevels(covars)
  if (!'Study.ID' %in% names(covars)) covars$Study.ID <- as.character(seq_len(nrow(covars)))
  incomp <- covars[!complete.cases(covars), Study.ID]
  covars <- covars[!Study.ID %in% incomp]
  if (!is.null(DT.y.m)) DT.y.m <- DT.y.m[!Study.ID %in% incomp]   # Not called for NBS
  setkey(covars, Study.ID)

  # Swap the outcome and measure variables, if outcome is not a network metric
  if (!is.null(outcome)) {
    DT.y.m[, eval(outcome) := covars[, get(outcome)], by=region]
    covars[, eval(outcome) := NULL]

    # Graph-level case is simple; only 1 design matrix
    if (level == 'graph') {
      covars[, eval(measure) := DT.y.m[, get(measure)]]

    # Vertex-level has 1 design matrix per region, with the measure changing for each
    } else if (level == 'vertex') {
      DT.X.m <- merge(DT.y.m, covars, by='Study.ID')
      setcolorder(DT.X.m, c('Study.ID', 'region', names(covars)[-1], measure))
      DT.X.m[, eval(outcome) := NULL]

      # Get all design matrices into a 3-D array
      X <- setNames(vector('list', DT.X.m[, nlevels(region)]), DT.X.m[, levels(region)])
      for (rgn in DT.X.m[, levels(region)]) {
        X[[rgn]] <- brainGraph_GLM_design(DT.X.m[region == rgn, !'region', with=F], ...)
      }
      X <- abind::abind(X, along=3)
    }
    DT.y.m[, eval(measure) := NULL]
  } else {
    outcome <- measure
  }

  # Get design matrix and contrast + column names
  if (is.null(X)) X <- brainGraph_GLM_design(covars, ...)
  if (is.vector(con.mat)) con.mat <- t(con.mat)
  stopifnot(ncol(X) == ncol(con.mat))
  tmp <- contrast_names(con.mat, con.type, con.name)
  con.mat <- tmp$con.mat; con.name <- tmp$con.name; nC <- tmp$nC
  if (is.null(colnames(con.mat))) colnames(con.mat) <- colnames(X)

  return(list(covars=covars, incomp=incomp, X=X, con.mat=con.mat, con.name=con.name,
              nC=nC, DT.y.m=DT.y.m, outcome=outcome))
}

#' Create contrast names for GLM analysis
#'
#' Simple helper function to generate contrast names for GLM functions.
#' @keywords internal
#' @rdname glm_helpers

contrast_names <- function(con.mat, con.type, con.name) {
  nC <- ifelse(con.type == 't', nrow(con.mat), 1)
  if (is.null(con.name)) {
    if (!is.null(rownames(con.mat))) {
      con.name <- rownames(con.mat)
      if (con.type == 'f') con.name <- rownames(con.mat)[1]
    } else {
      con.name <- rownames(con.mat) <- paste('Contrast', seq_len(nC))
    }
  } else {
    if (con.type == 't' & length(con.name) < nC) {
      con.name <- c(con.name, paste('Contrast', seq_len(nC)[-(1:length(con.name))]))
      rownames(con.mat) <- con.name
    } else if (con.type == 'f') {
      con.name <- con.name[1]
      rownames(con.mat) <- c(con.name, rep('', nrow(con.mat) - 1))
    }
  }
  return(list(con.mat=con.mat, con.name=con.name, nC=nC))
}

################################################################################
# MODEL FITTING FUNCTIONS
################################################################################

#' Helper function for GLM fitting
#'
#' @param DT A data.table with all the necessary data; namely, \code{Study.ID},
#'   \code{region} (which is just \code{graph} if \code{level='graph'}), and the
#'   outcome measure(s)
#' @param mykey The \code{key} to key by (to differentiate NBS and other GLM
#'   analyses)
#' @inheritParams GLM
#' @keywords internal
#' @aliases glm_fit_helper
#' @rdname glm_helpers

glm_fit_helper <- function(DT, X, con.type, con.mat, alt, outcome, mykey, alpha=NULL) {
  ESS <- numer <- stat <- se <- p <- contrast <- ci.low <- ci.high <- NULL
  dfR <- nrow(X) - ncol(X)

  XtX <- solve(crossprod(X))
  if (con.type == 'f') {
    CXtX <- solve(con.mat %*% XtX %*% t(con.mat))
    rkC <- qr(con.mat)$rank
    DT.lm <- DT[, brainGraph_GLM_fit_f(X, get(outcome), dfR, con.mat, rkC, CXtX), by=mykey]
    DT.lm[, ESS := numer * rkC]
    DT.lm[, stat := (numer / (se / dfR))]
    DT.lm[, numer := NULL]
    DT.lm[, p := pf(stat, rkC, dfR, lower.tail=FALSE)]
  } else if (con.type == 't') {
    pfun <- switch(alt,
                   two.sided=function(stat, df) 2 * pt(abs(stat), df, lower.tail=FALSE),
                   less=function(stat, df) pt(stat, df),
                   greater=function(stat, df) pt(stat, df, lower.tail=FALSE))

    DT.lm <- vector('list', nrow(con.mat))
    for (i in seq_along(DT.lm)) {
      DT.lm[[i]] <- DT[, brainGraph_GLM_fit_t(X, get(outcome), XtX, con.mat[i, , drop=FALSE]), by=mykey]
    }
    DT.lm <- rbindlist(DT.lm, idcol='contrast')
    DT.lm[, stat := gamma / se]
    DT.lm[, p := pfun(stat, dfR)]
    if (!is.null(alpha)) {
      DT.lm[, ci.low := gamma - qt(alpha / 2, dfR, lower.tail=F) * se]
      DT.lm[, ci.high := gamma + qt(1 - (alpha / 2), dfR) * se]
    }
  }
  return(DT.lm)
}

#' Fit linear models for t contrasts
#'
#' \code{brainGraph_GLM_fit_t} fits a linear model for t-contrasts (i.e.,
#' uni-dimensional contrasts) and returns the contrasts of parameter estimates,
#' standard errors, t-statistics, and p-values. If a contrast matrix is
#' supplied, it will return the above values for each row of the matrix.
#'
#' For speed purposes (if it is called from \code{\link{brainGraph_GLM}} and
#' permutation testing is done), this function does not do argument checking.
#'
#' @param y Numeric vector; the outcome variable
#' @param XtX Numeric matrix
#' @inheritParams GLM
#' @importFrom RcppEigen fastLmPure
#'
#' @name GLMfit
#' @aliases brainGraph_GLM_fit_t
#' @rdname brainGraph_GLM_fit
#'
#' @return \code{brainGraph_GLM_fit_t} - A list containing:
#'   \item{gamma}{The contrast of parameter estimates}
#'   \item{se}{The standard error}
#'
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

brainGraph_GLM_fit_t <- function(X, y, XtX, con.mat) {
  est <- fastLmPure(X, y, method=2)
  b <- est$coefficients
  gamma <- con.mat %*% b
  sigma.squared <- est$s^2
  var.covar <- sigma.squared * XtX
  se <- sqrt(diag(con.mat %*% tcrossprod(var.covar, con.mat)))

  list(gamma=as.numeric(gamma), se=se)
}

#' Fit linear models for f contrasts
#'
#' \code{brainGraph_GLM_fit_f} fits a linear model for f-contrasts (i.e.,
#' multi-dimensional contrasts) and returns the \emph{extra sum of squares} due
#' to the full model, the sum of squared errors of the full model, the
#' f-statistic, and associated p-value.
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
#' @name GLMfit
#' @aliases brainGraph_GLM_fit_f
#' @rdname brainGraph_GLM_fit

brainGraph_GLM_fit_f <- function(X, y, dfR, con.mat, rkC, CXtX) {
  est <- fastLmPure(X, y, method=2)
  b <- as.matrix(est$coefficients)
  gamma <- con.mat %*% b
  SSEF <- as.numeric(crossprod(est$residuals))

  numer <- as.numeric(t(gamma) %*% CXtX %*% gamma / rkC)
  list(numer=numer, se=SSEF, contrast=1)
}

################################################################################
# S3 METHODS FOR "bg_GLM"
################################################################################

#' Print a summary from brainGraph_GLM analysis
#'
#' The \code{summary} method prints the results, only for which
#' \eqn{p < \alpha}, where \code{alpha} comes from the \code{bg_GLM} object.
#' \dQuote{Simple} P-values are used by default, but you may change this to the
#' FDR-adjusted or permutation P-values via the function argument \code{p.sig}.
#' You may also choose to subset by \emph{contrast}.
#'
#' @param object,x A \code{bg_GLM} object
#' @param p.sig Character string specifying which p-value to use for displaying
#'   significant results (default: \code{p})
#' @param contrast Integer specifying the contrast to plot/summarize; defaults
#'   to showing results for all contrasts
#' @param digits Integer specifying the number of digits to display for p-values
#' @param print.head Logical indicating whether or not to print only the first
#'   and last 5 rows of the statistics tables (default: \code{TRUE})
#' @export
#' @method summary bg_GLM
#' @rdname glm

summary.bg_GLM <- function(object, p.sig=c('p', 'p.fdr', 'p.perm'), contrast=NULL,
                           digits=max(3L, getOption('digits') - 2L), print.head=TRUE, ...) {
  stopifnot(inherits(object, 'bg_GLM'))
  Outcome <- threshold <- NULL
  object$p.sig <- match.arg(p.sig)
  object$contrast <- contrast
  object$digits <- digits
  object$print.head <- print.head
  DT.sum <- object$DT[get(object$p.sig) < object$alpha]
  DT.sum[, Outcome := NULL]
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
  Outcome <- threshold <- contrast <- NULL
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
#'   plot; only relevant if \code{level='vertex'} (default: \code{NULL})
#' @param which Integer vector indicating which of the 6 plots to print to the
#'   plot device (default: \code{c(1:3, 5)})
#' @export
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra arrangeGrob
#' @importFrom gridExtra grid.arrange
#' @method plot bg_GLM
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

plot.bg_GLM <- function(x, region=NULL, which=c(1L:3L, 5L), ...) {
  stopifnot(inherits(x, 'bg_GLM'))

  # Local function to plot for a single region
  plot_single <- function(X, y, region, outcome) {
    leverage <- resid.std <- cook <- ind <- mark <- fit <- leverage.tr <- NULL
    est <- fastLmPure(X, y, method=2)
    dt.p1 <- data.table(fit=est$fitted.values, resid=est$residuals)
    H <- X %*% solve(crossprod(X)) %*% t(X)
    dt.p1[, leverage := diag(H)]
    dt.p1[, resid.std := resid / (est$s * sqrt(1 - leverage))]
    dt.p1[, cook := (resid.std^2 / ncol(X)) * (leverage / (1 - leverage))]
    dt.p1[, ind := as.character(.I)]
    dt.p1[, mark := ifelse(abs(resid) < mean(resid) + 2 * sd(resid), 0, 1)]
    dt.p1[mark == 0, ind := '']
    dt.p1[, mark := as.factor(mark)]

    diagPlots <- vector('list', length=6)
    # 1. Resids vs fitted
    diagPlots[[1]] <- ggplot(dt.p1, aes(x=fit, y=resid)) +
      geom_point(aes(shape=mark)) +
      geom_text_repel(aes(x=fit, y=resid, label=ind), size=3) +
      stat_smooth(method='loess', se=FALSE, span=1, col='red') +
      geom_hline(yintercept=0, lty=3) +
      theme(plot.title=element_text(hjust=0.5),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title='Residuals vs Fitted', x='Fitted values', y='Residuals')

    # 2. QQ-plot
    dt.p2 <- copy(dt.p1)
    setkey(dt.p2, resid.std)
    dt.p2[, x := qnorm(ppoints(resid.std))]
    diagPlots[[2]] <- ggplot(dt.p2, aes(x=x, y=resid.std)) +
      geom_text_repel(aes(x=x, y=resid.std, label=ind), size=3) +
      geom_line(aes(x=x, y=x), col='gray50', lty=3) +
      geom_point(aes(shape=mark)) +
      theme(plot.title=element_text(hjust=0.5),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title='Normal Q-Q', x='Theoretical Quantiles', y='Sample Quantiles')

    # 3. Scale-Location plot
    diagPlots[[3]] <- ggplot(dt.p1, aes(x=fit, y=sqrt(abs(resid.std)))) +
      geom_point(aes(shape=mark)) +
      geom_text_repel(aes(x=fit, y=sqrt(abs(resid.std)), label=ind), size=3) +
      stat_smooth(method='loess', se=FALSE, col='red') +
      ylim(c(0, NA)) +
      theme(plot.title=element_text(hjust=0.5),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title='Scale-Location', x='Fitted values', y=expression(sqrt(Standardized~residuals)))

    # 4. Cook's distance
    diagPlots[[4]] <- ggplot(dt.p1, aes(x=seq_len(nrow(X)), y=cook)) +
      geom_bar(stat='identity', position='identity') +
      geom_text_repel(aes(y=cook, label=ind), size=3) +
      theme(plot.title=element_text(hjust=0.5),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title='Cook\'s distance', x='Obs. number', y='Cook\'s distance')

    # 5. Residual vs Leverage plot
    diagPlots[[5]] <- ggplot(dt.p1, aes(x=leverage, y=resid.std)) +
      geom_point(aes(shape=mark)) +
      geom_text_repel(aes(x=leverage, y=resid.std, label=ind), size=3) +
      stat_smooth(method='loess', se=FALSE, col='red') +
      geom_vline(xintercept=0, lty=3, col='gray50') +
      geom_hline(yintercept=0, lty=3, col='gray50') +
      xlim(c(0, NA)) +
      theme(plot.title=element_text(hjust=0.5),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title='Residuals vs Leverage', x='Leverage', y='Standardized residuals')

    # 6. Cook's dist vs leverage
    dt.p1[, leverage.tr := leverage / (1 - leverage)]
    diagPlots[[6]] <- ggplot(dt.p1, aes(x=leverage.tr, y=cook)) +
      geom_point(aes(shape=mark)) +
      geom_text_repel(aes(x=leverage.tr, y=cook, label=ind), size=3) +
      geom_abline(slope=seq(0, 3, 0.5), lty=2, col='gray50') +
      xlim(c(0, NA)) +
      theme(plot.title=element_text(hjust=0.5, size=9),
            legend.position='none',
            axis.text.y=element_text(hjust=0.5, angle=90)) +
      labs(title=expression("Cook's dist vs Leverage "*h[ii] / (1 - h[ii])),
           x=expression(Leverage~h[ii]), y='Cook\'s distance')

    prows <- ifelse(length(which) == 1, 1, 2)
    p.all <- arrangeGrob(grobs=diagPlots[which], nrow=prows,
                         top=paste0('Outcome: ', outcome, '    Region: ', region),
                         bottom=paste0('y ~ ', paste(colnames(X), collapse=' + ')))
    return(p.all)
  }

  X <- x$X
  if (x$level == 'graph') {
    p.all <- plot_single(X, x$y, region='graph-level', x$outcome)
    return(p.all)
  } else if (x$level == 'vertex') {
    if (is.null(region)) {
      region <- x$DT[, levels(region)]
    }
    p.all <- sapply(region, function(x) NULL)
    for (z in region) {
      if (all(x$y[, z] == 0)) next
      p.all[[z]] <- plot_single(X, x$y[, z], z, x$outcome)
    }
    return(p.all)
  }
}
