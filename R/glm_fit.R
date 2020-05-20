# There is a lot of re-used code across these functions to avoid excessive
# branching when being repeatedly run in a loop (i.e., avoiding code such as
# "if (perm.method == ...) {FOO} else if (...) "). This is because, for a
# given function call, many of the inputs (in some cases, all) do not change and
# should not be checked in every permutation.

#===============================================================================
# Model fitting functions
#===============================================================================

#' Fit design matrices to one or multiple outcomes
#'
#' These are the \dQuote{base} model-fitting functions that solve the
#' \emph{least squares problem} to estimate model coefficients, residuals, etc.
#' for brain network data.
#'
#' @section Parameter estimation:
#' These functions use the \emph{QR} decomposition to calculate the least
#' squares solution which is the same as the base \code{\link[stats]{lm}}
#' function. If we substitute \eqn{X = QR} in the standard normal equations, the
#' equation to be solved reduces to
#' \deqn{X^T X \hat{\beta} = X^T y \Rightarrow R \hat{\beta} = Q^T y}
#'
#' Since \code{R} is an \emph{upper-triangular} matrix, we can use the
#' \code{\link{backsolve}} function which is a bit faster than
#' \code{\link{solve}}. In some cases, the \code{fastLmBG*} functions are about
#' as fast or faster (particularly when \code{X} is not permuted) as one in
#' which the normal equations are solved directly; additionally, using the
#' \emph{QR} method affords greater numerical stability.
#'
#' @section Different scenarios:
#' There are a few different scenarios for fitting models of the data, with a
#' separate function for each:
#' \describe{
#'   \item{fastLmBG}{The main function for when there is a single design matrix
#'     \eqn{X} and any number of outcome variables \eqn{Y}.}
#'   \item{fastLmBG_3d}{Fits models when there is a different design matrix
#'     \eqn{X} for each region and a single outcome variable \eqn{Y}, which in
#'     this case will be a column matrix.}
#'   \item{fastLmBG_3dY}{Fits models when there is both a different design
#'     matrix \eqn{X} and outcome variable \eqn{Y} for each region. Occurs under
#'     permutation for the Freedman-Lane, ter Braak, and Still-White methods.}
#'   \item{fastLmBG_3dY_1p}{Fits models when there is both a different design
#'     and outcome variable for each region, and also when \eqn{X} is a rank-1
#'     matrix (i.e., it has 1 column). Only occurs under permutation with the
#'     Still-White method if there is a single regressor of interest.}
#' }
#'
#' In the last case above, model coefficients are calculated by simple (i.e.,
#' non-matrix) algebra.
#'
#' @section Improving speed/efficiency:
#' Speed/efficiency gains will be vast for analyses in which there is a single
#' design matrix \eqn{X} for all regions, there are multiple outcome variables
#' (i.e., vertex-level analysis), and the permutation method chosen does
#' not permute \eqn{X}. Specifically, these are \emph{Freedman-Lane}, \emph{ter
#' Braak}, and \emph{Manly} methods. Therefore, the QR decomposition, the
#' \eqn{Q} and \eqn{R} matrices, and the \dQuote{unscaled covariance matrix}
#' (which is \eqn{(X^T X)^{-1}}) only need to be calculated once for the entire
#' analysis. Other functions (e.g., \code{lm.fit}) would recalculate these for
#' each permutation.
#'
#' Furthermore, this (and the other model fitting functions in the package) will
#' likely only work in models with full rank. I sacrifice proper error checking
#' in favor of speed, but hopefully any issues with the model will be identified
#' prior to the permutation step. Finally, the number of observations, model
#' rank, number of outcome variables, and degrees of freedom will not change and
#' therefore do not need to be recalculated (although these probably amount to a
#' negligible speed boost).
#'
#' In case there are multiple design matrices, or the permutation method
#' permutes the design, then the QR decomposition will need to be calculated
#' each time anyway. For these cases, I use more simplified functions
#' \code{qr_Q2} and \code{qr_R2} to calculate the \eqn{Q} and \eqn{R} matrices,
#' and then the fitted values, residuals, and residual standard deviation are
#' calculated at the same time (whereas \code{lm.fit} and others would calculate
#' these each time).
#'
#' @param X Design matrix or 3D array of design matrices
#' @param Y Numeric matrix; there should be 1 column for each outcome variable
#'   (so that in a graph-level analysis, this is a column matrix)
#' @param QR,Q,R The QR decomposition(s) and Q and R matrix(es) of the design
#'   matrix(es). If \code{X} is a 3D array, these should be \emph{lists}
#' @param n,p,ny,dfR Integers; the number of observations, model \emph{rank},
#'   number of regions/outcome variables, and residual degrees of freedom
#' @param XtXinv Numeric matrix or array; the inverse of the cross-product of
#'   the design matrix(es)
#' @export
#'
#' @return A list with elements
#'   \item{coefficients}{Parameter estimates}
#'   \item{rank}{Model rank}
#'   \item{df.residual}{Residual degrees of freedom}
#'   \item{residuals}{Model residuals}
#'   \item{sigma}{The residual standard deviation, or \emph{root mean square
#'     error (RMSE)}}
#'   \item{fitted.values}{Model fitted values}
#'   \item{qr}{The design matrix QR decomposition(s)}
#'   \item{cov.unscaled}{The \dQuote{unscaled covariance matrix}}
#' @name GLM fits
#' @rdname glm_fit
#' @family GLM functions
#' @seealso randomise
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

fastLmBG <- function(X, Y, QR=qr.default(X), Q=qr_Q2(QR, n=n, p=p), R=qr_R2(QR, p),
                     n=dim(X)[1L], p=QR$rank, ny=dim(Y)[2L], dfR=n-p, XtXinv=inv(QR)) {

  beta <- backsolve(R, crossprod(Q, Y), p)
  yhat <- X %*% beta
  ehat <- Y - yhat
  s <- if (ny == 1L) sum(ehat^2) else .colSums(ehat^2, n, ny)

  list(coefficients=beta, rank=p, df.residual=dfR, residuals=ehat,
       sigma=sqrt(s / dfR), fitted.values=yhat, qr=QR, cov.unscaled=XtXinv)
}

#' @param runX Character vector of the regions for which the design matrix is
#'   not singular
#' @export
#' @rdname glm_fit

fastLmBG_3d <- function(X, Y, runX, QR=qr(X[, , runX, drop=FALSE]),
                        Q=lapply(QR, qr_Q2, n=n, p=p), R=lapply(QR, qr_R2, p),
                        n=dim(X)[1L], p=QR[[1L]]$rank, ny=length(runX), dfR=n-p,
                        XtXinv=inv(QR)) {

  beta <- vapply(runX, function(r) backsolve(R[[r]], crossprod(Q[[r]], Y), p), numeric(p))
  yhat <- vapply(runX, function(r) X[, , r] %*% beta[, r], numeric(n))
  ehat <- c(Y) - yhat
  s <- .colSums(ehat^2, n, ny) / dfR

  list(coefficients=beta, rank=p, df.residual=dfR, residuals=ehat,
       sigma=sqrt(s), fitted.values=yhat, qr=QR, cov.unscaled=XtXinv)
}

#' @export
#' @rdname glm_fit

fastLmBG_3dY <- function(X, Y, runX, QR=qr(X[, , runX, drop=FALSE]),
                         Q=lapply(QR, qr_Q2, n=n, p=p), R=lapply(QR, qr_R2, p),
                         n=dim(X)[1L], p=QR[[1L]]$rank, ny=length(runX), dfR=n-p,
                         XtXinv=inv(QR)) {

  beta <- vapply(runX, function(r) backsolve(R[[r]], crossprod(Q[[r]], Y[, r]), p), numeric(p))
  yhat <- vapply(runX, function(r) X[, , r] %*% beta[, r], numeric(n))
  ehat <- Y - yhat
  s <- .colSums(ehat^2, n, ny) / dfR

  list(coefficients=beta, rank=p, df.residual=dfR, residuals=ehat,
       sigma=sqrt(s), fitted.values=yhat, qr=QR, cov.unscaled=XtXinv)
}

#' @export
#' @rdname glm_fit

fastLmBG_3dY_1p <- function(X, Y, runX, QR=qr(X[, , runX, drop=FALSE]),
                            Q=lapply(QR, qr_Q2, diag(1L, n, 1L), n, 1L),
                            R=lapply(QR, function(r) r$qr[1L]), n=dim(X)[1L], p=1L,
                            ny=length(runX), dfR=n-1L, XtXinv=inv(QR)) {

  Qarr <- sapply(Q, function(x) x, simplify='array')
  beta <- colSums(Qarr[, 1L, ] * Y) / unlist(R, recursive=FALSE, use.names=FALSE)
  yhat <- X[, 1L, ] %*% diag(beta)
  ehat <- Y - yhat
  dim(beta) <- c(1L, ny)
  dimnames(beta)[[2L]] <- dimnames(yhat)[[2L]] <- dimnames(ehat)[[2L]] <- runX
  s <- .colSums(ehat^2, n, ny) / dfR

  list(coefficients=beta, rank=p, df.residual=dfR, residuals=ehat,
       sigma=sqrt(s), fitted.values=yhat, qr=QR, cov.unscaled=XtXinv)
}

#===============================================================================
# Contrast-based statistics
#===============================================================================

#' Calculate GLM statistics for T and F contrasts
#'
#' \code{fastLmBG_t} and \code{fastLmBG_f} calculate contrast-based statistics
#' for T or F contrasts, respectively. It accepts any number of \emph{contrasts}
#' (i.e., a multi-row contrast matrix).
#'
#' @section Contrast-based statistics:
#' The \emph{contrast of parameter estimates}, \eqn{\gamma}, for T contrasts is
#' \deqn{\gamma = C \hat{\beta}}
#' where \eqn{C} is the contrast matrix with size \eqn{k \times p} (where
#' \eqn{k} is the number of contrasts) and \eqn{\hat{\beta}} is the matrix of
#' parameter estimates with size \eqn{p \times r} (where \eqn{r} is the number
#' of regions). For F contrasts, the effect size is the \emph{extra sum of
#' squares} and is calculated as
#' \deqn{\gamma (C (X^T X)^{-1} C^T)^{-1} \gamma^T}
#' The \emph{standard error} of a T contrast is
#' \deqn{\sqrt{\hat{\sigma} (X^T X)^{-1}}}
#' where \eqn{\hat{\sigma}} is the \emph{residual standard deviation} of the
#' model and the second term is the unscaled covariance matrix. The standard
#' error for F contrasts is simply the \emph{residual sum of squares}. P-values
#' and FDR-adjusted P-values (across regions) are also calculated. Finally, if
#' \eqn{\alpha} is provided for T contrasts, confidence limits are calculated.
#'
#' @param fits List object output by one of the model fitting functions (e.g.,
#'   \code{fastLmBG})
#' @inheritParams GLM
#' @export
#' @return \code{fastLmBG_t} -- A multidimensional array with the third
#'   dimension equalling the number of contrasts; each matrix contains the
#'   contrast of parameter estimates, standard error of the contrast,
#'   T-statistics, P-values, FDR-adjusted P-values, and confidence intervals (if
#'   \code{alpha} is given)
#' @rdname glm_fit

fastLmBG_t <- function(fits, contrasts, alternative=c('two.sided', 'less', 'greater'), alpha=NULL) {
  is3d <- inherits(fits$qr, 'list')
  dfR <- fits$df.residual
  XtX <- fits$cov.unscaled

  gamma <- t(contrasts %*% fits$coefficients)

  if (is3d) {
    CXtX <- cxtxfun_3d('t', transpose=FALSE)(contrasts, xtx=XtX, rkC=qr.default(contrasts)$rank)
    se <- sqrt(fits$sigma^2 * t(CXtX))
  } else {
    CXtX <- contrasts %*% tcrossprod(XtX, contrasts)
    se <- sqrt(outer_vec(fits$sigma^2, diag(CXtX)))
  }

  stat <- gamma / se
  alt <- match.arg(alternative)
  pfun <- switch(alt,
                 two.sided=function(stat, df) 2 * pt(abs(stat), df, lower.tail=FALSE),
                 less=function(stat, df) pt(stat, df),
                 greater=function(stat, df) pt(stat, df, lower.tail=FALSE))
  p <- pfun(stat, dfR)
  pfdr <- apply(p, 2L, p.adjust, 'fdr')
  dim(pfdr) <- dim(p)

  res <- abind::abind(gamma=gamma, se=se, stat=stat, p=p, p.fdr=pfdr, along=1.5)
  if (!is.null(alpha)) {
    t.limits <- c(1, -1) * qt(alpha / 2, dfR)
    cilow <- gamma + t.limits[1L] * se
    cihigh <- gamma + t.limits[2L] * se
    res <- abind::abind(res, ci.low=cilow, ci.high=cihigh, along=2L)
  }
  res
}

#' @param rkC,nC Integers; the rank of the contrast matrix and number of
#'   contrasts, respectively (for F contrasts)
#' @export
#' @return \code{fastLmBG_f} -- A numeric matrix with columns for the effect
#'   size, standard error, F statistic, P-values, and FDR-adjusted P-values
#' @rdname glm_fit

fastLmBG_f <- function(fits, contrasts, rkC=NULL, nC=length(contrasts)) {
  is3d <- inherits(fits$qr, 'list')
  dfR <- fits$df.residual
  XtX <- fits$cov.unscaled
  if (is.null(rkC)) rkC <- vapply(contrasts, function(x) qr.default(x)$rank, integer(1L))

  # This will be a list of "nC" matrices w/ dims "nregions X rkC"
  gamma <- lapply(contrasts, function(x) t(x %*% fits$coefficients))

  if (is3d) {
    cxtxfun <- cxtxfun_3d('f')
    rgns <- dimnames(XtX)[[3L]]
    CXtX <- vector('list', nC)
    for (x in seq_len(nC)) {
      CXtX[[x]] <- cxtxfun(contrasts[[x]], XtX, rkC[x], dim(XtX)[3L])
      dimnames(CXtX[[x]])[[3L]] <- rgns
    }
    ESS <- mapply(function(x, y)
                  vapply(rgns, function(r) rowSums(x[r, ] %*% y[, , r] * x[r, ]), numeric(1L)),
                  gamma, CXtX)
  } else {
    CXtX <- lapply(contrasts, function(x) inv(x %*% XtX, x, transpose=TRUE))
    ESS <- mapply(function(x, y) rowSums(x %*% y * x), gamma, CXtX)
  }

  SSEF <- matrix(fits$sigma^2 * dfR, dim(ESS)[1L], nC)
  stat <- t(t(ESS) / rkC) / (SSEF / dfR)
  p <- t(pf(t(stat), rkC, dfR, lower.tail=FALSE))
  pfdr <- apply(p, 2L, p.adjust, 'fdr')
  dim(pfdr) <- dim(p)
  abind::abind(ESS=ESS, se=SSEF, stat=stat, p=p, p.fdr=pfdr, along=1.5)
}

#===============================================================================
# Statistics for "randomise" and "randomise_3d"
#===============================================================================
#-----------------------------------------------------------
# T statistic functions
#-----------------------------------------------------------

#' Calculate T and F statistics when called from randomise
#'
#' These functions only call \code{fastLmBG} (the main model-fitting functions)
#' and calculates the T- or F-statistics.
#'
#' @noRd

fastLmBG_t_rand <- function(X, Y, contrasts, QR, Q, R, n, p, ny, dfR, XtXinv, CXtX, rkC) {
  all_fits <- fastLmBG(X, Y, QR, Q, R, n, p, ny, dfR, XtXinv)
  gamma <- t(contrasts %*% all_fits$coefficients)
  se <- sqrt(tcrossprod(all_fits$sigma^2, CXtX))  # Only works if CXtX is the  diagonal
  gamma / se
}

#' @param fitfun Function that will calculate model coefficients, residuals,
#'   etc. Can be one of: \code{fastLmBG_3d}, \code{fastLmBG_3dY},
#'   \code{fastLmBG_3dY_1p}
#' @noRd

fastLmBG_t_rand_3d <- function(fitfun, X, Y, contrasts, runX, QR, Q, R, n, p, ny, dfR, XtXinv, CXtX, rkC) {
  all_fits <- fitfun(X, Y, runX, QR, Q, R, n, p, ny, dfR, XtXinv)
  gamma <- t(contrasts %*% all_fits$coefficients)
  se <- sqrt(all_fits$sigma^2 * t(CXtX))
  gamma / se
}

#-----------------------------------------------------------
# F statistic functions
#-----------------------------------------------------------

#' @noRd
fastLmBG_f_rand <- function(X, Y, contrast, QR, Q, R, n, p, ny, dfR, XtXinv, CXtX, rkC) {
  all_fits <- fastLmBG(X, Y, QR, Q, R, n, p, ny, dfR, XtXinv)
  gamma <- t(contrast %*% all_fits$coefficients)
  ESS <- rowSums(gamma %*% CXtX * gamma)
  (ESS / rkC) / all_fits$sigma^2
}

#' @noRd
fastLmBG_f_rand_3d <- function(fitfun, X, Y, contrasts, runX, QR, Q, R, n, p, ny, dfR, XtXinv, CXtX, rkC) {
  all_fits <- fitfun(X, Y, runX, QR, Q, R, n, p, ny, dfR, XtXinv)
  gamma <- t(contrasts %*% all_fits$coefficients)
  ESS <- vapply(seq_len(ny), function(r) rowSums(gamma[r, ] %*% CXtX[, , r] * gamma[r, ]), numeric(1L))
  (ESS / rkC) / all_fits$sigma^2
}

#-----------------------------------------------------------
# Helper functions
#-----------------------------------------------------------

#' Add rows/columns for regions not fitted
#'
#' \code{add_nans} adds rows/columns (or higher dimensions) to model fit data
#' for regions which were skipped (due to having a singular design matrix,
#' usually).
#'
#' @param fits List object output by one of the model fitting functions (e.g.,
#'   \code{\link{fastLmBG}})
#' @param dimX Integer vector containing the dimensions of the original design
#'   matrix/array (including singular designs)
#' @param namesX List of character vectors containing the dimension names from
#'   the original design matrix/array
#' @param runX Character vector of regions for which models were fit
#'
#' @return \code{add_nans} -- the original \code{fits} object with \code{NaN} or
#'   \code{NA} inserted
#' @keywords internal
#' @rdname glm_helpers

add_nans <- function(fits, dimX, namesX, runX=names(fits$qr)) {
  # Initialize the full matrices/arrays
  coeffs <- se <- matrix(NaN, fits$rank, dimX[3L], dimnames=namesX[c(2L, 3L)])
  resids <- fvals <- matrix(NaN, dimX[1L], dimX[3L], dimnames=namesX[c(1L, 3L)])
  sig <- setNames(rep_len(NaN, dimX[3L]), namesX[[3L]])
  QRx <- setNames(vector('list', dimX[3L]), namesX[[3L]])
  XtX <- vcv <- array(NaN, dim=dimX[c(2L, 2L, 3L)], dimnames=namesX[c(2L, 2L, 3L)])

  # Fill in and re-assign
  coeffs[, runX] <- fits$coefficients
  fits$coefficients <- coeffs
  resids[, runX] <- fits$residuals
  fits$residuals <- resids
  fvals[, runX] <- fits$fitted.values
  fits$fitted.values <- fvals
  names(fits$sigma) <- runX
  sig[runX] <- fits$sigma
  fits$sigma <- sig
  QRx[runX] <- fits$qr
  fits$qr <- QRx
  XtX[, , runX] <- fits$cov.unscaled
  fits$cov.unscaled <- XtX
  if (!is.null(fits$var.covar)) {
    vcv[, , runX] <- fits$var.covar
    fits$var.covar <- vcv
  }
  se[, runX] <- fits$se
  fits$se <- se

  fits
}

#' Return a function to calculate the CXtX matrix/array for multiple designs
#'
#' \code{cxtxfun_3d} returns a function that calculates the \dQuote{CXtX}
#' matrix/array, used to calculate the standard error of a contrast. The
#' function signature will be \code{f(contrast, xtx, rkC, ny)}.
#'
#' The \dQuote{CXtX} matrix/array for T-contrasts is the diagonal of
#' \deqn{C (X^T X)^{-1} C^T}
#' For F-contrasts, it is the \emph{inverse} of the matrix:
#' \deqn{(C (X^T X)^{-1} C^T)^{-1}}
#' where in both cases, \eqn{C} is the contrast matrix and \eqn{X} is the design
#' matrix/array.
#'
#' For T-contrasts, the function will return a numeric vector/matrix with
#' dimensions \eqn{k \times r} where \eqn{k} is the \emph{number} of contrasts
#' (i.e., the number of rows in the contrast matrix) and \eqn{r} is the number
#' of regions. For F-contrasts, the function would return a numeric array with
#' dimensions \eqn{k \times k \times r}, where \eqn{k} is the rank of the
#' contrast matrix and \eqn{r} is the number of regions. If there is a single
#' design for all regions, it will be a \eqn{k \times k} matrix.
#'
#' @param con.type Either \code{'t'} or \code{'f'}
#' @param transpose Logical indicating whether to transpose the output of the
#'   selected function. Ignored for F-contrasts. Should be \code{FALSE} for
#'   T-contrasts when the \dQuote{Guttman} partition method is used. Default:
#'   \code{TRUE}
#' @return \code{cxtxfun_3d} -- A function with arguments for the contrast
#'   (numeric matrix), the \dQuote{unscaled covariance} array, the (row) rank of
#'   the contrast, and the number of regions in the analysis (only used for
#'   F-contrasts)
#' @keywords internal
#' @rdname glm_helpers

cxtxfun_3d <- function(con.type, transpose=TRUE) {
    if (con.type == 't') {
      if (isTRUE(transpose)) {
        fun <- function(contrast, xtx, rkC, ny) {
          t(apply(xtx, 3L, function(r) diag_sq(tcrossprod(contrast %*% r, contrast), rkC)))
        }
      } else {
        fun <- function(contrast, xtx, rkC, ny) {
          apply(xtx, 3L, function(r) diag_sq(tcrossprod(contrast %*% r, contrast), rkC))
        }
      }

    # The F-contrast function uses the last argument (ny)
    } else if (con.type == 'f') {
      fun <- function(contrast, xtx, rkC, ny) {
        res <- apply(xtx, 3L, function(r) inv(contrast %*% r, contrast, transpose=TRUE))
        dim(res) <- c(rkC, rkC, ny)
        res
      }
    }

  fun
}
