#' Extract model fit statistics from a bg_GLM object
#'
#' These functions extract or calculate model fit statistics of a
#' \code{bg_GLM} object. These can be found in the output from
#' \code{\link[stats]{summary.lm}}.
#'
#' These mimic the same functions that operate on \code{\link{lm}} objects, and
#' include:
#' \describe{
#'   \item{coef}{Regression coefficients (parameter estimates)}
#'   \item{confint}{Confidence intervals (by default, 95\%) for parameter
#'     estimates}
#'   \item{fitted}{Fitted (mean) values; i.e., the design matrix multiplied by
#'     the parameter estimates, \eqn{X \hat{\beta}}}
#'   \item{residuals}{Model residuals; i.e., the response/outcome variable minus
#'     the \emph{fitted} values. Partial residuals can also be calculated}
#'   \item{deviance}{Model deviance, or the \emph{residual sum of squares}}
#'   \item{coeff_determ}{Calculate the \emph{coefficient of determination} (or
#'     \eqn{R^2}), adjusted or unadjusted}
#'   \item{df.residual}{Residual degrees of freedom}
#'   \item{sigma}{Residual standard deviation, sometimes called the \emph{root
#'     mean squared error (RMSE)}}
#'   \item{vcov}{Variance-covariance matrix of the model parameters}
#' }
#'
#' \code{coeff_table} returns model coefficients, standard errors, T-statistics,
#' and P-values for all model terms and regions in a \code{bg_GLM} object. This
#' is the same as running \code{summary(x)$coefficients} for a \code{lm} object.
#'
#' @note \code{sigma} -- The denominator is \emph{not} the number of
#'   observations, but rather the model's \emph{residual degrees of freedom}.
#' @note When calculating \emph{partial residuals}, the parameter estimates are
#'   \emph{not} re-calculated after removing one of the model terms.
#'
#' @param object A \code{bg_GLM} object
#' @param ... Unused
#' @export
#' @return A named numeric vector, matrix, or array, depending on the function:
#'   \item{coef}{Matrix in which rownames are parameter names and column names
#'     are regions}
#'   \item{fitted,residuals}{Matrix in which rownames are Study ID's and column
#'     names are regions. If \code{type='partial'}, an array is returned in
#'     which columns are \emph{terms} and the 3rd dimension are regions}
#'   \item{deviance,coeff_determ,sigma}{Numeric vector with elements for each
#'     region}
#'   \item{df.residual}{Single integer; the degrees of freedom}
#'   \item{confint,vcov,coeff_table}{Numeric array; the extent of the third
#'     dimension equals the number of regions}
#'
#' @name GLM statistics
#' @rdname glm_stats
#' @seealso \code{\link{GLM}}, \code{\link[car]{Anova}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

coef.bg_GLM <- function(object, ...) {
  if (!is.null(object$coefficients)) return(object$coefficients)

  X <- object$X
  if (is.null(p <- object$rank)) p <- dim(X)[2L]
  if (length(dim(X)) == 3L) {
    dimX <- dim(X)
    runX <- object$runX
    if (is.null(QR <- object$qr)) QR <- qr(X[, , runX])
    Q <- lapply(QR, qr_Q2, n=dimX[1L], p=p)
    R <- lapply(QR, qr_R2, p=p)
    coeffs <- matrix(NaN, dimX[2L], dimX[3L], dimnames=dimnames(X)[c(2L, 3L)])
    coeffs[, runX] <- vapply(runX, function(r) backsolve(R[[r]], crossprod(Q[[r]], object$y), p),
                             numeric(p))
  } else {
    if (is.null(QR <- object$qr)) QR <- qr.default(X)
    Q <- qr_Q2(QR, p=p)
    R <- qr_R2(QR, p)
    coeffs <- backsolve(R, crossprod(Q, object$y), p)
  }

  return(coeffs)
}

# Helper function to calculate standard errors of coefficients
coef_se <- function(object, pnames=variable.names(object)) {
  if (!is.null(object$se)) return(object$se)
  dimX <- dim(object$X)
  if (length(dimX) == 3L) {
    std_err <- apply(vcov(object), 3L, function(x) sqrt(diag(x)))[pnames, , drop=FALSE]
  } else {
    XtX <- object$cov.unscaled
    if (is.null(XtX)) XtX <- inv(object$X)[pnames, pnames, drop=FALSE]
    sig <- sigma(object)^2
    std_err <- sqrt(outer_vec(diag(XtX), sig))
  }
  std_err
}

#' @param parm Vector of parameters to calculate confidence intervals for.
#'   Default is to use all parameters
#' @param level The confidence level. Default: \code{0.95}
#' @export
#' @rdname glm_stats

confint.bg_GLM <- function(object, parm, level=0.95, ...) {
  pnames <- variable.names(object)
  dfR <- object$df.residual <- df.residual(object)
  object$coefficients <- coef(object)
  dimC <- dim(object$coefficients)
  std_err <- coef_se(object, pnames)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  t.crit <- qt(a, dfR)
  ci <- array(object$coefficients, dim=c(dimC, 2L)) + outer(std_err, t.crit)
  ci <- aperm(ci, c(1L, 3L, 2L))
  ci <- ci[parm, , , drop=FALSE]
  dimnames(ci)[[2L]] <- paste(format(100 * a, trim=TRUE, scientific=FALSE, digits=3L), '%')
  return(ci)
}

#' @export
#' @rdname glm_stats

fitted.bg_GLM <- function(object, ...) {
  if (!is.null(object$fitted.values)) return(object$fitted.values)
  coeffs <- coef(object)

  # outcome != measure && level == 'vertex' (multiple design mats)
  dimX <- dim(object$X)
  if (length(dimX) == 3L) {
    runX <- object$runX
    fits <- matrix(NaN, dimX[1L], dimX[3L], dimnames=dimnames(object$X)[c(1L, 3L)])
    fits[, runX] <- vapply(runX, function(x) object$X[, , x] %*% coeffs[, x],
                           numeric(dimX[1L]))
  } else {
    fits <- object$X %*% coeffs
  }
  return(fits)
}

#' @param type Character string specifying the type of residuals to return.
#'   Default: \code{'response'}
#' @export
#' @rdname glm_stats

residuals.bg_GLM <- function(object, type=c('response', 'partial'), ...) {
  cf <- object$coefficients <- coef(object)
  type <- match.arg(type)
  if (type == 'partial') {
    object$residuals <- object$fitted.values <- NULL
    regions <- region.names(object)
    tm <- terms(object)
    if (any(grepl('Intercept', names(tm)))) tm <- tm[-grep('Intercept', names(tm))]
    resids <- array(0, dim=c(nobs(object), length(regions), length(tm)))
    for (i in seq_along(tm)) {
      object$coefficients <- cf[-tm[[i]], , drop=FALSE]
      resids[, , i] <- residuals(object[, -tm[[i]]])
    }
    resids <- aperm(resids, c(1L, 3L, 2L))
    dimnames(resids) <- list(case.names(object), names(tm), regions)

  } else {
    if (!is.null(object$residuals)) return(object$residuals)
    fits <- fitted(object)
    y <- if (dim(object$y)[2L] == 1L) c(object$y) else object$y
    resids <- y - fits
  }
  return(resids)
}

#' @export
#' @rdname glm_stats
deviance.bg_GLM <- function(object, ...) colSums(residuals(object)^2)

#' @param adjusted Logical indicating whether to calculate the adjusted
#'   R-squared. Default: \code{FALSE}
#' @export
#' @rdname glm_stats

coeff_determ <- function(object, adjusted=FALSE) {
  stopifnot(inherits(object, 'bg_GLM'))
  SSR <- deviance(object)
  if (dim(object$y)[2L] > 1L) {
    SStot <- diag(tcrossprod(t(object$y) - colMeans(object$y)))
  } else {
    SStot <- c(crossprod(object$y - mean(object$y)))
  }
  numer <- if (isTRUE(adjusted)) SSR / df.residual(object) else SSR
  denom <- if (isTRUE(adjusted)) SStot / (nobs(object) - 1L) else SStot
  1 - numer / denom
}

#' @export
#' @method df.residual bg_GLM
#' @rdname glm_stats

df.residual.bg_GLM <- function(object, ...) {
  if (!is.null(object$df.residual)) return(object$df.residual)
  dimX <- dim(object$X)
  dimX[1L] - dimX[2L]
}

#' @export
#' @rdname glm_stats

sigma.bg_GLM <- function(object, ...) {
  if (!is.null(object$sigma)) return(object$sigma)
  sqrt(deviance(object) / df.residual(object))
}

#' @export
#' @rdname glm_stats

vcov.bg_GLM <- function(object, ...) {
  if (!is.null(object$var.covar)) return(object$var.covar)
  sig <- sigma(object)^2
  X <- object$X
  dimX <- dim(X)

  XtX <- object$cov.unscaled
  if (is.null(XtX)) {
    if (length(dimX) == 3L) {
      namesX <- dimnames(X)
      dimV <- dimX[c(2L, 2L, 3L)]
      namesV <- namesX[c(2L, 2L, 3L)]
      runX <- object$runX
      XtX <- array(NaN, dim=dimV, dimnames=namesV)
      XtX[, , runX] <- inv(X[, , runX, drop=FALSE])
    } else {
      XtX <- inv(X)
    }
  }
  if (length(dimX) == 3L) {
    XtX <- aperm(XtX, c(3L, 1L, 2L))
    vcv <- sig * XtX
  } else {
    vcv <- outer(sig, XtX)
  }
  aperm(vcv, c(2L, 3L, 1L))
}

#' @param CI Logical indicating whether to include confidence intervals of
#'   parameter estimates in the coefficient summary table. Default: \code{FALSE}
#' @export
#' @importFrom abind abind
#' @rdname glm_stats

coeff_table <- function(object, CI=FALSE, level=0.95) {
  stopifnot(inherits(object, 'bg_GLM'))
  dfR <- object$df.residual <- df.residual(object)
  cf <- object$coefficients <- coef(object)
  se <- object$se <- coef_se(object)
  tab <- abind(cf, se, along=1.5)
  cnames <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
  if (isTRUE(CI)) {
    ci <- confint(object, level=level)
    tab <- abind(tab, ci, along=2L)
    cnames <- append(cnames, dimnames(ci)[[2L]], after=2L)
  }
  t.stat <- tab[, 1L, , drop=FALSE] / tab[, 2L, , drop=FALSE]
  tab <- abind(tab, t.stat, along=2L)
  tab <- abind(tab, 2 * pt(abs(t.stat), dfR, lower.tail=FALSE), along=2L)
  dimnames(tab)[2L:3L] <- list(cnames, region.names(object))
  return(tab)
}

# Convenience function to calculate leave-one-out resid SD
sigma_loo <- function(model) {
  if (!is.null(model$sigma.loo)) return(model$sigma.loo)
  dfR <- df.residual(model)
  resids2 <- residuals(model)^2
  var.hat <- t(colSums(resids2) - t(resids2))
  sqrt(var.hat / (dfR - 1L))
}

#' Influence measures for a bg_GLM object
#'
#' These functions compute common (leave-one-out) diagnostics for the models in
#' a \code{bg_GLM} object.
#'
#' The \code{influence} method calculates all diagnostics present in
#' \code{\link[stats]{lm.influence}} and
#' \code{\link[stats]{influence.measures}}, consisting of the following
#' functions:
#' \describe{
#'   \item{rstandard}{Standardized residuals. Choosing \code{type='predictive'}
#'     returns leave-one-out cross validation residuals. The \dQuote{PRESS}
#'     statistic can be calculated as \code{colSums(resids.p^2)}}
#'   \item{rstudent}{Studentized residuals}
#'   \item{hatvalues}{The \emph{leverage}, or the diagonal of the
#'     \emph{hat/projection matrix}}
#'   \item{cooks.distance}{Cook's distance}
#'   \item{dffits.bg_GLM}{The change in fitted values when deleting
#'     observations}
#'   \item{dfbeta}{The change in parameter estimates (coefficients) when
#'     deleting observations}
#'   \item{dfbetas}{The \emph{scaled} change in parameter estimates}
#'   \item{covratio.bg_GLM}{The covariance ratios, or the change in the
#'     determinant of the covariance matrix of parameter estimates when deleting
#'     observations}
#' }
#'
#' @section Outlier values:
#' Each variable has a different criterion for determining outliers. In the
#' following: \code{x} is the influence variable (for \code{DFBETA}, the
#' criterion applies to all DFBETAs); \code{k} is the number of columns of the
#' design matrix; \code{dfR} is the residual degrees of freedom; and \code{n} is
#' the number of observations.
#' \describe{
#'   \item{DFBETAs}{If \eqn{|x| > 1}}
#'   \item{DFFITs}{If \eqn{|x| > 3 \sqrt{k / dfR}}}
#'   \item{covratio}{If \eqn{|1 - x| > (3k / dfR)}}
#'   \item{cook}{If \eqn{F_{k, dfR}(x) > 0.5}}
#'   \item{hat}{If \eqn{x > 3k / n}}
#' }
#' The return object of \code{influence} has a \code{print} method which will
#' list the subjects/variables/regions for which an outlier was detected.
#'
#' @param model A \code{bg_GLM} object
#' @param type The type of standardized residuals. Default: \code{'sd.1'}
#' @param ... Unused
#' @export
#' @return Most influence functions return a numeric matrix in which rownames
#'   are Study ID's and column names are regions. \code{dfbeta} and
#'   \code{dfbetas} return a numeric array in which each column is a parameter
#'   estimate and the 3rd dimension is for each region. \code{influence} returns
#'   a list with class \code{infl.bg_GLM} and elements:
#'   \item{infmat}{Numeric array (like \code{dfbeta}) with DFBETAs, DFFITs,
#'     covratios, Cook's distance, and hat values}
#'   \item{is.inf}{Logical array of the same data as \code{infmat}; values of
#'     \code{TRUE} indicate the subject-variable-region combination is an
#'     outlier value}
#'   \item{f}{The model \emph{formula}}
#'   \item{sigma}{The leave-one-out residual standard deviation}
#'   \item{wt.res}{Model residuals}
#'
#' @name GLM influence measures
#' @rdname glm_influence
#' @seealso \code{\link{GLM}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

rstandard.bg_GLM <- function(model, type=c('sd.1', 'predictive'), ...) {
  model$residuals <- residuals(model)
  hat <- hatvalues(model)
  type <- match.arg(type)
  denom <- if (type == 'sd.1') t(sigma(model) * sqrt(1 - t(hat))) else 1 - hat
  model$residuals / denom
}

#' @export
#' @rdname glm_influence

rstudent.bg_GLM <- function(model, ...) {
  res <- model$residuals <- residuals(model)
  hat <- hatvalues(model)
  var.hat <- sigma_loo(model)
  res / (var.hat * sqrt(1 - hat))
}

#' @export
#' @rdname glm_influence

hatvalues.bg_GLM <- function(model, ...) {
  if (!is.null(model$hatvalues)) return(model$hatvalues)
  X <- model$X
  dimX <- dim(X)
  namesX <- dimnames(X)
  if (!is.null(model$cov.unscaled)) {
    XtX <- model$cov.unscaled
    if (length(dimX) == 3L) {
      hv <- matrix(NaN, dimX[1L], dimX[3L], dimnames=namesX[c(1L, 3L)])
      for (r in model$runX) {
        hv[, r] <- colSums(tcrossprod(XtX[, , r], X[, , r]) * t(X[, , r]))
      }
    } else {
      hv <- colSums(tcrossprod(XtX, X) * t(X))
      rgn <- region.names(model)
      hv <- matrix(hv, dimX[1L], length(rgn), dimnames=list(namesX[[1L]], rgn))
    }
    return(hv)
  }
  if (length(dimX) == 3L) {
    runX <- model$runX
    A <- U <- array(NaN, dim=dimX[c(2L, 2L, 3L)], dimnames=namesX[c(2L, 2L, 3L)])
    A[, , runX] <- apply(X[, , runX, drop=FALSE], 3L, crossprod)
    U[, , runX] <- apply(A[, , runX, drop=FALSE], 3L, chol.default)
    Z <- array(0, dim=dimX[c(2L, 1L, 3L)])
    for (i in seq_len(dimX[3L])) {
      Z[, , i] <- forwardsolve(t(U[, , i]), t(X[, , i]))
    }
    res <- colSums(Z^2)
    dimnames(res) <- c(namesX[1L], namesX[3L])
  } else {
    A <- crossprod(X)
    U <- chol.default(A)
    Z <- forwardsolve(t(U), t(X))
    res <- colSums(Z^2)
    rgn <- region.names(model)
    res <- matrix(res, dimX[1L], length(rgn), dimnames=list(namesX[[1L]], rgn))
  }
  return(res)
}

#' @export
#' @rdname glm_influence

cooks.distance.bg_GLM <- function(model, ...) {
  hat <- hatvalues(model)
  res <- model$residuals <- residuals(model)
  sd <- sigma(model)
  hat * (res / t(sd * (1 - t(hat))))^2 / model$rank
}

#' @export
#' @rdname glm_influence

dffits.bg_GLM <- function(model) {
  res <- model$residuals <- residuals(model)
  hat <- hatvalues(model)
  sig <- sigma_loo(model)
  res * sqrt(hat) / (sig * (1 - hat))
}

#' @export
#' @rdname glm_influence

dfbeta.bg_GLM <- function(model, ...) {
  X <- model$X
  dimX <- dim(X)
  namesX <- dimnames(X)
  regions <- region.names(model)
  Nv <- length(regions)

  resids <- residuals(model)
  hat <- hatvalues(model)
  mult <- resids / (1 - hat)
  multArr <- array(mult, dim=c(dimX[1L], Nv, dimX[2L]))
  multArr <- aperm(multArr, c(1L, 3L, 2L))

  if (length(dimX) == 3L) {
    runX <- model$runX
    numer <- array(NaN, dim=dimX, dimnames=namesX)
    XtX <- model$cov.unscaled
    if (is.null(XtX)) {
      XtX <- array(NaN, dim=dimX[c(2L, 2L, 3L)], dimnames=namesX[c(2L, 2L, 3L)])
      XtX[, , runX] <- inv(X[, , runX, drop=FALSE])
    }
    numer[, , runX] <- vapply(runX, function(r) t(tcrossprod(XtX[, , r], X[, , r])),
                              numeric(dimX[1L] * dimX[2L]))
  } else {
    XtX <- model$cov.unscaled
    if (is.null(XtX)) XtX <- inv(X)
    numer <- t(tcrossprod(XtX, X))
    numer <- array(numer, dim=c(dimX, Nv), dimnames=c(namesX, list(regions)))
  }
  numer * multArr
}

#' @export
#' @rdname glm_influence

dfbetas.bg_GLM <- function(model, ...) {
  resids <- model$residuals <- residuals(model)
  sig <- sigma_loo(model)

  X <- model$X
  dimX <- dim(X)
  XtX <- model$cov.unscaled
  if (length(dimX) == 3L) {
    namesX <- dimnames(X)
    runX <- model$runX
    if (is.null(XtX)) {
      XtX <- array(NaN, dim=dimX[c(2L, 2L, 3L)], dimnames=namesX[c(2L, 2L, 3L)])
      XtX[, , runX] <- inv(model$X[, , runX, drop=FALSE])
      model$cov.unscaled <- XtX
    }
    xxi <- apply(XtX, 3L, function(x) sqrt(diag(x)))
    dfb <- if (is.null(model$dfbeta)) dfbeta(model) else model$dfbeta

    dfbs <- array(NaN, dim=dimX, dimnames=dimnames(dfb))
    for (i in runX) {
      dfbs[, , i] <- dfb[, , i] / tcrossprod(sig[, i], xxi[, i])
    }
  } else {
    if (is.null(XtX)) XtX <- model$cov.unscaled <- inv(X)
    xxi <- sqrt(diag(XtX))
    dfb <- if (is.null(model$dfbeta)) dfbeta(model) else model$dfbeta
    dfbs <- dfb / aperm(outer(sig, xxi), c(1L, 3L, 2L))
  }
  return(dfbs)
}

#' @export
#' @rdname glm_influence

covratio.bg_GLM <- function(model) {
  dfR <- model$df.residual <- df.residual(model)
  p <- model$rank
  if (is.null(p)) p <- dim(model$X)[2L]
  omh <- 1 - hatvalues(model)
  res <- model$residuals <- residuals(model)
  sig <- sigma_loo(model)
  e.star <- res / (sig * sqrt(omh))
  1 / (omh * ((e.star^2 + dfR - 1) / dfR)^p)
}

#' @param do.coef Logical indicating whether to calculate \code{dfbeta}
#' @param region Character string of the region(s) to return results for.
#'   Default is to calculate for all regions
#' @export
#' @rdname glm_influence
#' @importFrom abind abind

influence.bg_GLM <- function(model, do.coef=TRUE, region=NULL, ...) {
  n <- nobs(model)
  hat <- model$hatvalues <- hatvalues(model)
  model$residuals <- residuals(model)
  sig <- model$sigma.loo <- sigma_loo(model)
  dff <- dffits.bg_GLM(model)
  covr <- covratio.bg_GLM(model)
  cook <- cooks.distance(model)

  # Create an array and determine which observations are influential based on any metric
  vnames <- variable.names(model)
  k <- length(vnames)
  dfR <- n - k

  cnames <- c('dffit', 'cov.r', 'cook.d', 'hat')
  infmat <- abind(dff, covr, cook, hat, along=1.5)
  infl <- abind(abs(dff) > 3 * sqrt(k / dfR),
                abs(1 - covr) > (3 * k) / dfR,
                pf(cook, k, dfR) > 0.5,
                hat > (3 * k) / n,
                along=1.5)
  if (isTRUE(do.coef)) {
    coeffs <- dfbetas(model)
    cnames <- c(paste0('dfb.', vnames), cnames)
    infmat <- abind(coeffs, infmat, along=2L)
    infl <- abind(abs(coeffs) > 1, infl, along=2L)
  }
  dimnames(infmat)[[2L]] <- dimnames(infl)[[2L]] <- cnames
  if (!is.null(region)) {
    infmat <- infmat[, , region, drop=FALSE]
    infl <- infl[, , region, drop=FALSE]
    sig <- sig[, region, drop=FALSE]
    model$residuals <- model$residuals[, region, drop=FALSE]
  }
  ans <- list(infmat=infmat, is.inf=infl, f=formula(model), sigma=sig, wt.res=model$residuals)
  class(ans) <- c('infl.bg_GLM', class(ans))
  return(ans)
}

#' @method print infl.bg_GLM
#' @export

print.infl.bg_GLM <- function(x, region=NULL, ...) {
  sID <- getOption('bg.subject_id')
  Region <- total <- value <- NULL
  message('\nInfluence measures and outliers for a bg_GLM model with formula:')
  cat('  ', x$f, '\n\n')
  DT <- as.data.table(x$is.inf, key='V1')
  setnames(DT, c(sID, 'variable', 'Region', 'value'))
  if (!is.null(region)) DT <- DT[Region %in% region]
  DT <- DT[value == TRUE]
  DT.split <- split(DT, by='variable', keep.by=FALSE)

  if (dim(x$infmat)[3L] > 1L) {
    DT.split.wide <- lapply(DT.split, dcast, paste(sID, '~ Region'))
    for (nm in names(DT.split.wide)) {
      DT.split.wide[[nm]][, total := DT.split[[nm]][, .N, by=sID]$N]
      message('Variable: ', nm)
      print(DT.split.wide[[nm]])
      cat('\n')
    }
  } else {  # single region
    for (nm in names(DT.split)) {
      message('Variable: ', nm)
      cat('  ', paste(DT.split[[nm]][, as.character(get(sID))], collapse=', '), '\n')
      cat('\n')
    }
  }
  invisible(x)
}

#' Variance inflation factors for \code{bg_GLM} objects
#'
#' @param mod A \code{bg_GLM} object
#' @param ... Unused
#' @export
#' @return A named array of VIFs; names of the 3rd dimension are regions

vif.bg_GLM <- function(mod, ...) {
  stopifnot(inherits(mod, 'bg_GLM'))
  X <- mod$X
  dimX <- dim(X)
  cols <- terms(mod)
  tlabels <- names(cols)
  if ((hasInt <- any(grepl('Intercept', tlabels)))) {
    int <- grep('Intercept', tlabels)
    cols <- lapply(cols[-int], `-`, 1L)
  } else {
    warning('No intercept; VIFs may not be sensible.')
  }

  regions <- region.names(mod)
  kNumRegions <- length(regions)
  n.terms <- length(cols)
  if (n.terms < 2L) stop('Model contains fewer than 2 terms!')
  df <- lengths(cols, use.names=FALSE)
  terms_single <- which(df == 1L)
  terms_mult <- which(df > 1L)
  v <- mod$cov.unscaled
  if (length(dimX) == 3L) {
    runX <- mod$runX
    #TODO: to remove
    if (is.null(v)) {
      v <- array(NaN, dim=dimX[c(2L, 2L, 3L)], dimnames=dimnames(X)[c(2L, 2L, 3L)])
      v[, , runX] <- inv(X[, , runX, drop=FALSE])
    }
    R <- array(NaN, dim=dim(v), dimnames=dimnames(v))
    if (isTRUE(hasInt)) {
      v <- v[-int, -int, , drop=FALSE]
      R <- R[-int, -int, , drop=FALSE]
    }
    R[, , runX] <- apply(v[, , runX, drop=FALSE], 3L, cov2cor)
    detR <- setNames(rep(NaN, kNumRegions), regions)
    detR[runX] <- apply(R[, , runX, drop=FALSE], 3L, det)
    res <- array(NaN, dim=c(n.terms, 3L, kNumRegions), dimnames=list(NULL, NULL, regions))

    if (length(terms_single) > 0L) {
      for (rgn in runX) {
        res[terms_single, 1L, rgn] <- vapply(cols[terms_single], function(j)
                                             det(as.matrix(R[-j, -j, rgn])),
                                             numeric(1L))
      }
    }
    if (length(terms_mult) > 0L) {
      for (rgn in runX) {
        res[terms_mult, 1L, rgn] <- vapply(cols[terms_mult], function(j)
                                           det(R[j, j, rgn]) * det(as.matrix(R[-j, -j, rgn])),
                                           numeric(1L))
      }
    }
    res[, 1L, ] <- t(t(res[, 1L, ]) / detR)
  } else {
    if (is.null(v)) v <- inv(X)
    if (hasInt) v <- v[-int, -int, drop=FALSE]
    R <- cov2cor(v)
    detR <- det(R)
    res <- matrix(0, nrow=n.terms, ncol=3L)
    if (length(terms_single) > 0L) {
      res[terms_single, 1L] <- vapply(cols[terms_single], function(j) det(R[-j, -j, drop=FALSE]), numeric(1L))
    }
    if (length(terms_mult) > 0L) {
      res[terms_mult, 1L] <- vapply(cols[terms_mult], function(j)
                                    det(R[j, j, drop=FALSE]) * det(R[-j, -j, drop=FALSE]),
                                    numeric(1L))
    }
    res[, 1L] <- res[, 1L] / detR
    res <- array(res, dim=c(dim(res), kNumRegions))
  }
  dimnames(res) <- list(names(cols), c('GVIF', 'Df', 'GVIF^(1/(2*Df))'), regions)
  res[, 2L, ] <- df
  if (any(df != 1L)) res[, 3L, ] <- res[, 1L, ]^(1 / (2L * df))
  return(res)
}

#' @section ANOVA tables:
#' The \code{anova} method calculates the so-called \emph{Type III} test
#' statistics for a \code{bg_GLM} object. These standard ANOVA statistics
#' include: sum of squares, mean squares, degrees of freedom, F statistics, and
#' P-values. Additional statistics calculated are: \eqn{\eta^2}, partial
#' \eqn{\eta^2}, \eqn{\omega^2}, and partial \eqn{\omega^2} as measures of
#' \emph{effect size}.
#'
#' @param region Character vector indicating the region(s) to calculate ANOVA
#'   statistics for. Default: \code{NULL} (use all regions)
#' @export
#' @return \code{anova} returns a \emph{list} of tables of class \code{anova}
#' @rdname glm_stats

anova.bg_GLM <- function(object, region=NULL, ...) {
  regions <- if (is.null(region)) region.names(object) else region
  kNumRegions <- length(regions)
  RSS <- deviance(object)[regions]

  cols <- terms(object)
  object$coefficients <- object$fitted.values <- object$residuals <- NULL
  SSt <- vapply(cols, function(x) deviance(object[, -x])[regions], numeric(kNumRegions))
  ss <- cbind(SSt - RSS, Residuals=RSS)
  df <- c(lengths(cols), df.residual(object))
  ms <- t(t(ss) / df)

  nrows <- length(df)      # Number of terms + 1 (residuals)
  SSTot <- .rowSums(ss, kNumRegions, nrows)
  eta2 <- ss / SSTot
  eta2.part <- ss / (ss + RSS)
  MSE <- RSS / df[nrows]
  Ndf <- nobs(object) - df
  MS <- tcrossprod(MSE, df)
  omega2 <- (ss - MS) / (SSTot + MSE)
  omega2.part <- (ss - MS) / (ss + tcrossprod(MSE, Ndf))
  f <- ms / MSE
  p <- matrix(nrow=kNumRegions, ncol=nrows, dimnames=dimnames(ss))
  p[, -nrows] <- t(pf(t(f[, -nrows]), df[-nrows], df[nrows], lower.tail=FALSE))

  dfMat <- matrix(df, kNumRegions, nrows, byrow=TRUE)
  arr <- abind(ss, ms, dfMat, f, eta2, eta2.part, omega2, omega2.part, p, along=1.5)
  arr <- aperm(arr, c(3L, 2L, 1L))
  dimnames(arr)[[2L]] <- c('Sum Sq', 'Mean Sq', 'Df', 'F value', 'eta^2', 'Partial eta^2',
                           'omega^2', 'Partial omega^2', 'Pr(>F)')
  arr[nrows, 4L:8L, ] <- NA
  tab <- apply(arr, 3L, as.data.frame)
  for (i in regions) {
    attr(tab[[i]], 'heading') <- c('Anova Table (Type III tests)\n',
                                   paste('Response:', object$outcome))
    class(tab[[i]]) <- c('anova', class(tab[[i]]))
  }
  return(tab)
}

#' Model selection for bg_GLM objects
#'
#' These functions compute the log-likelihood and Akaike's \emph{An Information
#' Criterion (AIC)} of a \code{bg_GLM} object. See
#' \code{\link[stats]{logLik.lm}} and \code{\link[stats]{extractAIC}} for
#' details.
#'
#' The functions \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}} will also
#' work for \code{bg_GLM} objects because they each call \code{logLik}.
#'
#' @param object,fit A \code{bg_GLM} object
#' @param REML Logical indicating whether to return the \emph{restricted}
#'   log-likelihood. Default: \code{FALSE}
#' @param ... Unused
#' @export
#' @return \code{logLik} returns an object of class \code{logLik} with several
#'   attributes. \code{extractAIC} returns a numeric vector in which the first
#'   element is the \emph{equivalent degrees of freedom} and the remaining are
#'   the AIC's for each region
#'
#' @name GLM model selection
#' @rdname glm_model_select

logLik.bg_GLM <- function(object, REML=FALSE, ...) {
  X <- object$X
  p <- object$rank
  N0 <- N <- nobs(object)

  if (isTRUE(REML)) {
    dimX <- dim(X)
    N <- N - p
    if (length(dimX) == 3L) {
      runX <- object$runX
      QR <- setNames(rep(NaN, dimX[3L]), dimnames(X)[[3L]])
      QR[runX] <- colSums(vapply(object$qr[runX], function(x)
                                 log(abs(diag(x$qr))), numeric(dimX[2L])))
    } else {
      QR <- sum(log(abs(diag(object$qr$qr))))
    }
  }
  val <- -0.5 * N * (log(2 * pi) + 1 - log(N) + log(deviance(object)))
  if (isTRUE(REML)) val <- val - QR
  structure(val, nall=N0, nobs=N, df=p+1L, class='logLik')
}

#' @param scale Should be left at its default
#' @param k Numeric; the weight of the equivalent degrees of freedom
#' @export
#' @rdname glm_model_select

extractAIC.bg_GLM <- function(fit, scale=0, k=2, ...) {
  n <- nobs(fit)
  edf <- n - df.residual(fit)
  RSS <- deviance(fit)
  dev <- if (scale > 0) RSS / scale - n else n * log(RSS / n)
  c(edf, dev + k * edf)
}
