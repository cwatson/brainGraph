#' Fit linear models at each vertex of a graph
#'
#' \code{brainGraph_GLM} specifies and fits a linear model at each vertex for a
#' given vertex measure (e.g. \emph{degree}) or at the graph-level (e.g.,
#' \emph{global efficiency}). Given a contrast matrix, it will calculate the
#' associated statistics.
#'
#' The input list of graphs \code{g.list} must not be nested; i.e., if you have
#' multiple groups, they will have to be combined into one list. See the code in
#' the Examples below.
#'
#' A \code{data.table} of covariates is required input; the first column
#' \emph{must} be named \emph{Study.ID}. Additionally, all graphs must
#' have a \emph{name} attribute (at the graph level) which matches the
#' \emph{Study.ID} for a given subject. If you create the design matrix
#' \code{X} yourself, you still must supply the covariates table so that
#' subjects can be correctly matched with their network data.
#'
#' Both t- and F-contrasts are allowed. You may supply a \emph{matrix} to the
#' argument \code{con.mat}. If you supply a multi-row matrix and you choose
#' \code{con.type="t"}, then statistics will be calculated for each contrast
#' individually. If you choose \code{con.type="f"}, in the result data table,
#' \code{ESS} stands for "extra sum of squares", the additional variance
#' explained for by the model parameters of interest (as determined by the
#' contrast matrix). Finally, the standard error in these tables is the sum of
#' squared errors of the \emph{full model}.
#'
#' Finally, you can calculate permutations of the data to build a null
#' distribution of the maximum statistic, to provide control over false
#' positives. The permutation strategy is that of Freedman & Lane (1983), and
#' is the same as that in FSL's \emph{randomise}.
#'
#' @param g.list A list of \code{igraph} graph objects for all subjects
#' @param covars A \code{data.table} of covariates
#' @param measure Character string of the graph measure of interest
#' @param con.mat Numeric matrix specifying the contrast(s) of interest; if
#'   only one contrast is desired, you can supply a vector
#' @param con.type Character string; either \code{'t'} or \code{'f'} (for t or
#'   F-statistics). Default: \code{'t'}
#' @param X Numeric matrix, if you wish to supply your own design matrix
#'   (default: \code{NULL})
#' @param con.name Character vector of the contrast name(s); if \code{con.mat}
#'   has row names, those will be used for reporting results (default:
#'   \code{NULL})
#' @param alternative Character string, whether to do a two- or one-sided test
#'   (default: \code{'two.sided'})
#' @param alpha Numeric; the significance level (default: 0.05)
#' @param level Character string; either \code{vertex} (default) or
#'   \code{graph}
#' @param permute Logical indicating whether or not to permute group labels
#'   (default: \code{FALSE})
#' @param N Integer; number of permutations to create (default: 5e3)
#' @param perms Matrix of permutations, if you would like to provide your own
#'   (default: \code{NULL})
#' @param long Logical indicating whether or not to return all permutation
#'   results (default: \code{FALSE})
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
#'   \item{removed}{A character vector of Study.ID's removed due to incomplete
#'     data (if any)}
#'   \item{perm}{A list containing: \emph{null.dist} (the null distribution of
#'     maximum statistics), \emph{thresh} (the statistic value corresponding
#'     to the \eqn{100 \times (1 - \alpha)}th\% percentile of the null
#'     distribution)}.
#'
#' @family GLM functions
#' @family Group analysis functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Freedman D & Lane D (1983). \emph{A nonstochastic interpretation
#'   of reported significance levels}. J Bus Econ Stat, 1(4):292-298.
#' @references Nichols TE & Holmes AP (2001). \emph{Nonparametric permutation
#'   tests for functional neuroimaging: A primer with examples.} Human Brain
#'   Mapping, 15(1):1-25.
#' @references Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE (2014).
#'   \emph{Permutation inference for the general linear model}. NeuroImage,
#'   92:381-397.
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
                           X=NULL, con.name=NULL,
                           alternative=c('two.sided', 'less', 'greater'),
                           alpha=0.05, level=c('vertex', 'graph'),
                           permute=FALSE, N=5e3, perms=NULL, long=FALSE, ...) {
  Study.ID <- region <- Outcome <- Covariate <- p.fdr <- p <- Contrast <- i <-
    stat <- p.perm <- perm <- contrast <- V1 <- incomp <- ci.low <- ci.high <- se <- ESS <- numer <- NULL

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- match.arg(alternative)
  if (ctype == 'f') alt <- 'two.sided'
  glmSetup <- setup_glm(covars, X, con.mat, ctype, con.name, ...)
  covars <- glmSetup$covars; X <- glmSetup$X; incomp <- glmSetup$incomp
  con.mat <- glmSetup$con.mat; con.name <- glmSetup$con.name; nC <- glmSetup$nC

  # Get the outcome variable(s) into a data.table
  level <- match.arg(level)
  if (level == 'vertex') {
    y <- t(vapply(g.list, vertex_attr, numeric(vcount(g.list[[1]])), measure))
    DT <- data.table(Study.ID=vapply(g.list, graph_attr, character(1), 'name'))
    DT <- cbind(DT, y)
    setnames(DT, 2:ncol(DT), V(g.list[[1]])$name)
    y <- as.matrix(DT[!Study.ID %in% incomp, -1])
  } else if (level == 'graph') {
    DT <- data.table(Study.ID=vapply(g.list, graph_attr, character(1), 'name'),
                     graph=vapply(g.list, graph_attr, numeric(1), measure))
    y <- DT[!Study.ID %in% incomp, graph]
  }
  DT <- DT[!Study.ID %in% incomp]
  setkey(DT, Study.ID)
  DT.m <- melt(DT, id.vars='Study.ID', variable.name='region', value.name=measure)

  # Do the model fitting/estimation
  glmFits <- glm_fit_helper(DT.m, X, ctype, con.mat, alt, measure, 'region')
  DT.lm <- glmFits$DT.lm; dfR <- glmFits$dfR

  if (ctype == 't') {
    DT.lm[, ci.low := gamma - qt(alpha / 2, dfR, lower.tail=F) * se]
    DT.lm[, ci.high := gamma + qt(1 - (alpha / 2), dfR) * se]
  }
  DT.lm[, Outcome := measure]
  DT.lm[, p.fdr :=  p.adjust(p, 'fdr'), by=contrast]
  for (i in seq_along(con.name)) DT.lm[contrast == i, Contrast := con.name[i]]
  setkey(DT.lm, contrast, region)

  # Permutation testing
  null.dist <- null.thresh <- NA
  if (isTRUE(permute)) {
    if (is.null(perms)) perms <- shuffleSet(n=nrow(X), nset=N)

    zeroregions <- DT.m[, .SD[all(get(measure) == 0)], by=region][, unique(as.character(region))]
    null.dist <- randomise(ctype, N, perms, DT.m[!region %in% zeroregions], nC, measure, X, con.mat, alt)
    sortfun <- switch(alt,
                      two.sided=function(x) {sort(abs(x))},
                      less=function(x) {sort(x, decreasing=TRUE)},
                      greater=function(x) {sort(x)})
    null.thresh <- null.dist[, sortfun(V1)[floor((1 - alpha) * N) + 1], by=contrast]
    compfun <- switch(alt,
                      two.sided=function(x, y) {sum(abs(x) >= abs(y), na.rm=TRUE)},
                      less=function(x, y) {sum(x <= y, na.rm=TRUE)},
                      greater=function(x, y) {sum(x >= y, na.rm=TRUE)})
    for (i in seq_along(con.name)) {
      DT.lm[list(i), p.perm := (compfun(null.dist[contrast == i, V1], stat) + 1) / (N + 1), by=region]
    }
  }

  if (isTRUE(long)) perm <- list(null.dist=null.dist, thresh=null.thresh)
  out <- list(level=level, X=X, y=y, outcome=measure, con.type=ctype, con.mat=con.mat,
              con.name=con.name, alt=alt, alpha=alpha, DT=DT.lm, removed=incomp,
              permute=permute, N=N, perm=perm)
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
#' This function performs several tasks: it removes unused levels from
#' \code{covars}, removes subjects with incomplete data, creates a design matrix
#' (if not supplied), and supplies names to the contrast matrix.
#'
#' @inheritParams brainGraph_GLM
#' @keywords internal
#' @name GLMhelpers
#' @aliases setup_glm
#' @rdname glm_helpers

setup_glm <- function(covars, X, con.mat, con.type, con.name, ...) {
  Study.ID <- NULL
  covars <- droplevels(covars)
  if (!'Study.ID' %in% names(covars)) covars$Study.ID <- as.character(seq_len(nrow(covars)))
  incomp <- covars[!complete.cases(covars), Study.ID]
  covars <- covars[!Study.ID %in% incomp]
  setkey(covars, Study.ID)

  if (is.null(X)) X <- brainGraph_GLM_design(covars, ...)
  if (is.vector(con.mat)) con.mat <- t(con.mat)
  stopifnot(ncol(X) == ncol(con.mat))
  if (is.null(con.name)) {
    if (!is.null(rownames(con.mat))) {
      con.name <- rownames(con.mat)
      if (con.type == 'f') con.name <- rownames(con.mat)[1]
    } else {
      con.name <- rownames(con.mat) <- paste('Contrast', seq_len(nrow(con.mat)))
    }
  } else {
    if (con.type == 't' & length(con.name) < nrow(con.mat)) {
      con.name <- c(con.name, paste('Contrast', seq_len(nrow(con.mat))[-(1:length(con.name))]))
      rownames(con.mat) <- con.name
    } else if (con.type == 'f') {
      con.name <- con.name[1]
      rownames(con.mat) <- c(con.name, rep('', nrow(con.mat) - 1))
    }
  }
  if (is.null(colnames(con.mat))) colnames(con.mat) <- colnames(X)

  nC <- ifelse(con.type == 't', nrow(con.mat), 1)
  return(list(covars=covars, incomp=incomp, X=X, con.mat=con.mat, con.name=con.name, nC=nC))
}

#' Partition a design matrix into columns of interest and nuisance
#'
#' @param M Numeric matrix; the full design matrix
#' @param con.mat Numeric matrix; the contrast matrix
#' @param method Character string; the method of partitioning (default:
#'   \code{beckmann})
#' @importFrom MASS Null
#' @keywords internal
#'
#' @return A list containing:
#'   \item{X}{Numeric matrix for the covariates of interest}
#'   \item{Z}{Numeric matrix for the nuisance covariates}
#'   \item{eCm}{The \emph{effective contrast}, equivalent to the original, for
#'     the partitioned model \code{[X, Z]} and considering all covariates}
#'   \item{eCx}{Same as \code{eCx}, but considering only \code{X}}
#' @references Guttman I. Linear Models: An Introduction. Wiley, New York, 1982.
#' @references Smith SM, Jenkinson M, Beckmann C, Miller K, Woolrich M (2007).
#'   \emph{Meaningful design and contrast estimability in fMRI.} NeuroImage,
#'   34(1):127-36.

partition <- function(M, con.mat, method=c('beckmann', 'guttman')) {
  method <- match.arg(method)
  if (method == 'guttman') {
    idx <- which(con.mat != 0, arr.ind=TRUE)[, 2]
    X <- M[, idx, drop=FALSE]
    Z <- M[, -idx, drop=FALSE]
    eCm <- cbind(con.mat[, idx, drop=FALSE], con.mat[, -idx, drop=FALSE])
  } else if (method == 'beckmann') {
    Q <- solve(crossprod(M))
    cdc <- solve(con.mat %*% Q %*% t(con.mat))
    X <- M %*% Q %*% t(con.mat) %*% cdc

    Cu <- MASS::Null(t(con.mat))
    Cv <- Cu - (t(con.mat) %*% cdc %*% con.mat %*% Q %*% Cu)
    Z <- M %*% Q %*% Cv %*% solve(t(Cv) %*% Q %*% Cv)
    eCm <- cbind(diag(ncol(X)), matrix(0, nrow=ncol(X), ncol=ncol(Z)))
  }
  eCx <- eCm[, 1:ncol(X), drop=FALSE]
  return(list(X=X, Z=Z, eCm=eCm, eCx=eCx))
}

#' Helper function to setup for randomise
#'
#' \code{setup_randomise} is used to setup the data/objects for any function
#' that does permutations for GLM-based analysis.
#'
#' The tasks performed by this function are: separate the design matrix into the
#' independent variables of interest and nuisance variables (based on the
#' contrast(s)), and calculate the new contrast(s) based on this new design
#' matrix.
#'
#' @param nC Integer; the number of contrasts
#' @inheritParams brainGraph_GLM
#' @keywords internal
#' @rdname randomise
#' @return A list containing:
#'   \item{Mp}{The full partitioned model, joined}
#'   \item{Rz}{The residual-forming matrix}
#'   \item{MtM}{The inverse of the cross product of the full model}
#'   \item{eC}{The \emph{effective contrast}, equivalent to the original, for
#'     the partitioned model \code{[X, Z]} and considering all covariates}

setup_randomise <- function(X, con.mat, nC) {
  n <- nrow(X)
  Mp <- MtM <- Rz <- eC <- vector('list', length=nC)

  for (j in seq_len(nC)) {
    if (nC == 1) {
      parts <- partition(X, con.mat, 'beckmann')
    } else {
      parts <- partition(X, con.mat[j, , drop=FALSE], 'beckmann')
    }
    Mp[[j]] <- with(parts, cbind(X, Z))
    MtM[[j]] <- solve(crossprod(Mp[[j]]))
    Hz <- with(parts, Z %*% solve(crossprod(Z)) %*% t(Z))
    Rz[[j]] <- diag(n) - Hz
    eC[[j]] <- parts$eCm
  }

  return(list(Mp=Mp, Rz=Rz, MtM=MtM, eC=eC))
}

#' Randomize and fit a model and find the maximum statistic
#'
#' @param DT \code{data.table} with outcome variables
#' @inheritParams brainGraph_GLM
#' @keywords internal

randomise <- function(ctype, N, perms, DT, nC, measure, X, con.mat, alternative) {
  i <- region <- numer <- se <- perm <- contrast <- NULL
  randMats <- setup_randomise(X, con.mat, nC)
  Mp <- randMats$Mp; Rz <- randMats$Rz; MtM <- randMats$MtM; eC <- randMats$eC
  dfR <- nrow(Mp[[1]]) - qr(randMats$Mp[[1]])$rank

  if (ctype == 'f') {
    CMtM <- solve(randMats$eC[[1]] %*% randMats$MtM[[1]] %*% t(randMats$eC[[1]]))
    rkC <- qr(randMats$eC[[1]])$rank
  }
  maxfun <- switch(alternative,
                   two.sided=function(x) {max(abs(x), na.rm=TRUE)},
                   less=function(x) {min(x, na.rm=TRUE)},
                   greater=function(x) {max(x, na.rm=TRUE)})
  null.dist <- vector('list', length=nC)
  perm.order <- rep(seq_len(N), each=DT[, length(unique(region))])

  for (j in seq_len(nC)) {
    # T-contrasts
    if (ctype == 't') {
      null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        DT[, brainGraph_GLM_fit_t(Mp[[j]], Rz[[j]][perms[i, ], ] %*% get(measure), MtM[[j]], eC[[j]]), by=region]
      }
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][, maxfun(gamma / se), by=perm][, !'perm']

    # F-contrasts
    } else if (ctype == 'f') {
      null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        DT[, brainGraph_GLM_fit_f(Mp[[j]], Rz[[j]][perms[i, ], ] %*% get(measure), dfR, eC[[j]], rkC, CMtM), by=region]
      }
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][, max(numer / (se / dfR), na.rm=TRUE), by=perm][, !'perm']
    }
  }

  null.dist <- rbindlist(null.dist, idcol='contrast')
  return(null.dist)
}

################################################################################
# MODEL FITTING FUNCTIONS
################################################################################

#' Helper function for GLM fitting
#' @param DT A data.table with all the necessary data
#' @param mykey The \code{key} to key by (to differentiate NBS and other GLM
#'   analyses)
#' @inheritParams brainGraph_GLM
#' @keywords internal
#' @aliases glm_fit_helper
#' @rdname glm_helpers

glm_fit_helper <- function(DT, X, con.type, con.mat, alt, measure, mykey) {
  ESS <- numer <- stat <- se <- p <- contrast <- NULL
  dfR <- nrow(X) - ncol(X)
  XtX <- solve(crossprod(X))
  if (con.type == 'f') {
    CXtX <- solve(con.mat %*% XtX %*% t(con.mat))
    rkC <- qr(con.mat)$rank
    DT.lm <- DT[, brainGraph_GLM_fit_f(X, get(measure), dfR, con.mat, rkC, CXtX), by=mykey]
    DT.lm[, ESS := numer * rkC]
    DT.lm[, stat := (numer / (se / dfR))]
    DT.lm[, numer := NULL]
    DT.lm[, p := pf(stat, rkC, dfR, lower.tail=FALSE)]
  } else if (con.type == 't') {
    pfun <- switch(alt,
                   two.sided=function(stat, df) {2 * pt(abs(stat), df, lower.tail=FALSE)},
                   less=function(stat, df) {pt(stat, df)},
                   greater=function(stat, df) {pt(stat, df, lower.tail=FALSE)})

    DT.lm <- vector('list', nrow(con.mat))
    for (i in seq_along(DT.lm)) {
      DT.lm[[i]] <- DT[, brainGraph_GLM_fit_t(X, get(measure), XtX, con.mat[i, , drop=FALSE]), by=mykey]
    }
    DT.lm <- rbindlist(DT.lm, idcol='contrast')
    DT.lm[, stat := gamma / se]
    DT.lm[, p := pfun(stat, dfR)]
  }
  return(list(DT.lm=DT.lm, dfR=dfR))
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
#' @inheritParams brainGraph_GLM
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
  list(numer=numer, se=SSEF)
}

#' Create a design matrix for linear model analysis
#'
#' \code{brainGraph_GLM_design} takes a \code{data.table} of covariates and
#' returns a \emph{design matrix} to be used in linear model analysis.
#'
#' There are three different ways to code factors: \emph{dummy}, \emph{effects},
#' or \emph{cell-means} (chosen by the argument \code{coding}). To understand
#' the difference, see Chapter 8 of the User Guide.
#'
#' Importantly, the default behavior (as of v2.1.0) is to convert all character
#' columns (excluding the Study ID column and any that you list in the
#' \code{binarize} argument) to factor variables. To change this, set
#' \code{factorize=FALSE}. So, if your covariates include multiple character
#' columns, but you want to convert \emph{Scanner} to binary instead of a
#' factor, you may still specify \code{binarize='Scanner'} and get the expected
#' result. \code{binarize} will convert the given
#' factor variable(s) into numeric variable(s), which is performed \emph{before}
#' mean-centering.
#'
#' The argument \code{mean.center} will mean-center (i.e., subtract the mean of
#' the entire dataset from each variable) any non-factor variables (including
#' any dummy/indicator covariates). This is done \emph{after} "factorizing" and
#' "binarizing".
#'
#' \code{int} specifies which variables should interact with one another. This
#' argument accepts both numeric (e.g., \emph{Age}) and factor variables (e.g.,
#' \emph{Sex}). All interaction combinations will be generated: if you supply 3
#' variables, all two-way and the single three-way interaction will be
#' generated. This variable \emph{must} have at least two elements.
#'
#' @param covars A \code{data.table} of covariates
#' @param coding Character string indicating how factor variables will be coded
#'   (default: \code{'dummy'})
#' @param factorize Logical indicating whether to convert \emph{character}
#'   columns into \emph{factor} (default: \code{TRUE})
#' @param mean.center Logical indicating whether to mean center non-factor
#'   variables (default: \code{FALSE})
#' @param binarize Character vector specifying the column name(s) of the
#'   covariate(s) to be converted from type \code{factor} to \code{numeric}
#'   (default: \code{NULL})
#' @param int Character vector specifying the column name(s) of the
#'   covariate(s) to test for an interaction (default: \code{NULL})
#' @export
#'
#' @return A numeric matrix
#'
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

brainGraph_GLM_design <- function(covars, coding=c('dummy', 'effects', 'cell.means'),
                                  factorize=TRUE, mean.center=FALSE, binarize=NULL, int=NULL) {

  covars <- copy(covars)
  covars[, Study.ID := as.character(Study.ID)]
  X <- matrix(1, nrow=nrow(covars), ncol=1)
  colnames(X) <- 'Intercept'

  if (isTRUE(factorize)) {
    cols <- names(which(sapply(covars, is.character)))
    cols <- cols[!is.element(cols, 'Study.ID')]
    cols <- setdiff(cols, binarize)
    for (z in cols) covars[, eval(z) := as.factor(get(z))]
  }
  if (!is.null(binarize)) {
    stopifnot(all(binarize %in% names(covars)))
    for (b in binarize) covars[, eval(b) := as.numeric(get(b)) - 1]
  }

  nums <- which(sapply(covars, is.numeric))
  if (isTRUE(mean.center)) {
    for (n in nums) covars[[n]] <- covars[[n]] - mean(covars[[n]])
  }
  if (length(nums) > 0) X <- cbind(X, as.matrix(covars[, nums, with=F]))

  factors <- which(sapply(covars, class) == 'factor')
  coding <- match.arg(coding)
  for (f in factors) {
    cov.name <- names(covars)[f]
    cov.levels <- covars[, levels(get(cov.name))]
    cov.vec <- covars[, as.numeric(get(cov.name))]

    if (coding == 'cell.means') {
      for (i in 1:max(cov.vec)) {
        cov.vec.sub <- ifelse(cov.vec == i, 1, 0)
        X <- cbind(X, cov.vec.sub)
      }
      if (all(X[, 1] == 1)) X <- X[, -1]  # Remove intercept term
      p <- ncol(X)
      colnames(X)[(p - max(cov.vec) + 1):p] <- paste0(cov.name, cov.levels)

    } else {
      code <- ifelse(coding == 'dummy', 0, -1)
      for (i in 2:max(cov.vec)) {
        cov.vec.sub <- ifelse(cov.vec == i, 1, 0)
        cov.vec.sub <- ifelse(cov.vec == 1, code, cov.vec.sub)
        X <- cbind(X, cov.vec.sub)
      }
      p <- ncol(X)
      colnames(X)[(p - (max(cov.vec) - 1) + 1):p] <- paste0(cov.name, cov.levels[-1])
    }
  }

  if (!is.null(int) & length(int) > 1) {
    stopifnot(all(int %in% names(covars)))
    intcomb <- combn(int, 2, simplify=FALSE)
    if (length(int) == 3) intcomb <- c(intcomb, combn(int, 3, simplify=FALSE))
    for (x in intcomb) X <- get_int(X, coding, names(factors), x)
  }

  return(X)
}

get_int <- function(X, coding, factors, int) {
  get_colnames <- function(string, X, factors) {
    cnames <- colnames(X)[grep(string, colnames(X))]
    if ( (!string %in% factors) && length(cnames) > 1) {
      cnames <- colnames(X)[grep(paste0('^', string, '$'), colnames(X))]
    }
    if (any(grepl(':', cnames))) cnames <- cnames[-grep(':', cnames)]
    return(cnames)
  }

  p <- ncol(X)
  intnames <- lapply(int, get_colnames, X, factors)

  if (length(int) == 3) {
    X <- cbind(X, X[, intnames[[1]]] * X[, intnames[[2]]] * X[, intnames[[3]]])
    colnames(X)[(p + 1):ncol(X)] <- paste(intnames[[1]], intnames[[2]], intnames[[3]], sep=':')
  } else {
    if (int[1] %in% factors && int[2] %in% factors && coding == 'cell.means') {
      for (j in seq_along(intnames[[1]])) {
        X <- cbind(X, X[, intnames[[1]][j]] * X[, intnames[[2]]])
        p2 <- ncol(X)
        colnames(X)[(p + 1):p2] <- paste(intnames[[1]][j], intnames[[2]], sep=':')
        p <- ncol(X)
      }
      X <- X[, -which(colnames(X) %in% intnames[[1]])]
      X <- X[, -which(colnames(X) %in% intnames[[2]])]
    } else {  # One of the int. terms is numeric
      X <- cbind(X, X[, intnames[[1]]] * X[, intnames[[2]]])
      p2 <- ncol(X)
      colnames(X)[(p + 1):p2] <- paste(intnames[[1]], intnames[[2]], sep=':')
      if (coding == 'cell.means') {
        if (!int[1] %in% factors) {
          X <- X[, -which(colnames(X) %in% intnames[[1]])]
        } else if (!int[2] %in% names(factors)) {
          X <- X[, -which(colnames(X) %in% intnames[[2]])]
        }
      }
    }
  }
  return(X)
}

################################################################################
# S3 METHODS FOR "bg_GLM"
################################################################################

#' Print a summary from brainGraph_GLM analysis
#'
#' \code{summary} prints the results from a \code{\link{brainGraph_GLM}}
#' analysis. It will only print results for which \eqn{p < \alpha}; you may
#' change this to the FDR-adjusted or permutation p-values via the function
#' argument \code{p.sig}.
#'
#' @param object A \code{bg_GLM} object
#' @param p.sig Character string specifying which p-value to use for displaying
#'   significant results (default: \code{p})
#' @param contrast Integer specifying the contrast to summarize; defaults to
#'   showing results for all contrasts
#' @param digits Integer specifying the number of digits to display for p-values
#' @param print.head Logical indicating whether or not to print only the first
#'   and last 5 rows of the statistics tables (default: \code{TRUE})
#' @param ... Unused
#' @export
#' @method summary bg_GLM

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

  if (isTRUE(object$permute)) {
    if (object$con.type == 't') {
      setcolorder(DT.sum, c('Contrast', 'region', 'gamma', 'ci.low', 'ci.high', 'se',
                            'stat', 'p', 'p.fdr', 'p.perm', 'contrast'))
    } else {
      setcolorder(DT.sum, c('Contrast', 'region', 'ESS', 'se', 'stat', 'p', 'p.fdr',
                            'p.perm', 'contrast'))
    }
    setnames(DT.sum, 'p.perm', 'p-value (perm.)')
  } else {
    if (object$con.type == 't') {
      setcolorder(DT.sum, c('Contrast', 'region', 'gamma', 'ci.low', 'ci.high', 'se',
                            'stat', 'p', 'p.fdr', 'contrast'))
    } else {
      setcolorder(DT.sum, c('Contrast', 'region', 'ESS', 'se', 'stat', 'p', 'p.fdr',
                            'contrast'))
    }
  }
  if (object$con.type == 't') {
    clp <- 100 * (1 - object$alpha)
    setnames(DT.sum, c('region', 'gamma', 'se', 'stat', 'p', 'p.fdr'),
             c('Region', 'Estimate', 'Std. error', 't value', 'p-value', 'p-value (FDR)'))
    setnames(DT.sum, c('ci.low', 'ci.high'), paste0(clp, '% CI ', c('low', 'high')))
  } else {
    setnames(DT.sum, c('region', 'ESS', 'se', 'stat', 'p', 'p.fdr'),
             c('Region', 'Extra Sum Sq.', 'Std. error', 'F value', 'p-value', 'p-value (FDR)'))
  }

  object$DT.sum <- DT.sum
  class(object) <- c('summary.bg_GLM', class(object))
  return(object)
}

#' @aliases summary.bg_GLM
#' @method print summary.bg_GLM
#' @keywords internal

print.summary.bg_GLM <- function(x, ...) {
  Outcome <- threshold <- contrast <- NULL
  title <- paste('brainGraph GLM results')
  message('\n', title, '\n', rep('-', getOption('width') / 2))
  cat('Level: ', x$level, '\n')
  cat('Graph metric of interest: ', x$outcome, '\n\n')

  cat('Contrast type: ', paste(toupper(x$con.type), 'contrast'), '\n')
  alt <- switch(x$alt,
                two.sided='C != 0',
                greater='C > 0',
                less='C < 0')
  cat('Alternative hypothesis: ', alt, '\n')
  cat('Contrast matrix: ', '\n')
  print(x$con.mat)

  if (length(x$removed) != 0) cat('\nSubjects removed due to incomplete data:\n', x$removed, '\n')
  if (isTRUE(x$permute)) {
    message('\n', 'Permutation analysis', '\n', rep('-', getOption('width') / 4))
    cat('# of permutations: ', prettyNum(x$N, ','), '\n')
    cat(paste0(toupper(x$con.type), '-statistic threshold (based on the null distribution):\n'))
    print(x$perm$thresh)
  }

  # Print results for each contrast
  message('\n', 'Statistics', '\n', rep('-', getOption('width') / 4))
  if (is.null(x$contrast)) {
    contrast <- x$DT[, unique(contrast)]
  } else {
    contrast <- x$contrast
  }
  for (i in contrast) {
    message(x$con.name[i])
    if (nrow(x$DT.sum[contrast == i]) == 0) {
      message('\tNo signficant results!\n')
    } else {
      if (isTRUE(x$print.head)) {
        print(x$DT.sum[contrast == i, !c('Contrast', 'contrast')], topn=5, nrows=10, digits=x$digits)
      } else {
        print(x$DT.sum[contrast == i, !c('Contrast', 'contrast')], digits=x$digits)
      }
      cat('\n')
    }
  }
  invisible(x)
}

#' Plot GLM diagnostics for a brain network
#'
#' Plots the GLM diagnostics (similar to that of
#' \code{\link[stats]{plot.lm}}) for the output of \code{\link{brainGraph_GLM}}.
#' There are a total of 6 possible plots, specified by the \code{which}
#' argument; the behavior is the same as in \code{\link[stats]{plot.lm}}. Please
#' see the help for that function.
#'
#' @param x A \code{bg_GLM} object
#' @param region Character string specifying which region's results to
#'   plot; only relevant if \code{level='vertex'} (default: \code{NULL})
#' @param which Integer vector indicating which of the 6 plots to print to the
#'   plot device (default: \code{c(1:3, 5)})
#' @param ... Unused
#' @export
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra arrangeGrob
#' @importFrom gridExtra grid.arrange
#' @method plot bg_GLM
#'
#' @return A list of \code{\link[ggplot2]{ggplot}} objects
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
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
