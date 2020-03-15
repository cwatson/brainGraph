# Helper function to convert design matrix to an array
design2array <- function(object) {
  X <- object$X
  dimX <- dim(X)
  if (length(dimX) == 2L) {
    regions <- region.names(object)
    X <- array(X, dim=c(dimX, length(regions)), dimnames=c(dimnames(X), list(regions)))
  }
  return(X)
}

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
#'   \item{deviance,coeff_determ,df.residual,sigma}{Numeric vector with elements
#'     for each region}
#'   \item{confint,vcov}{Numeric array; the extent of the third dimension equals
#'     the number of regions}
#'
#' @importFrom RcppEigen fastLmPure
#'
#' @name GLM statistics
#' @rdname glm_stats
#' @seealso \code{\link{GLM}}, \code{\link[car]{Anova}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

coef.bg_GLM <- function(object, ...) {
  if (!is.null(object$coefficients)) return(object$coefficients)

  # outcome != measure && level == 'vertex' (multiple design matrices)
  if (length(dim(object$X)) == 3L) {
    coeffs <- apply(object$X, 3L, function(x) fastLmPure(x, object$y, method=2L)$coefficients)
  } else {
    coeffs <- apply(object$y, 2L, function(y) fastLmPure(object$X, y, method=2L)$coefficients)
  }

  if (object$level == 'graph') dimnames(coeffs)[[2L]] <- 'graph'
  return(coeffs)
}

#' @param parm Vector of parameters to calculate confidence intervals for.
#'   Default is to use all parameters
#' @param level The confidence level. Default: \code{0.95}
#' @export
#' @rdname glm_stats

confint.bg_GLM <- function(object, parm, level=0.95, ...) {
  object$coefficients <- coef(object)
  dimC <- dim(object$coefficients)
  std_err <- apply(vcov(object), 3L, function(x) sqrt(diag(x)))
  pnames <- variable.names(object)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  t.crit <- qt(a, df.residual(object))
  pct <- paste(format(100 * a, trim=TRUE, scientific=FALSE, digits=3L), '%')
  ci <- array(object$coefficients, dim=c(dimC, 2L)) + outer(std_err, t.crit)
  ci <- aperm(ci, c(1L, 3L, 2L))
  ci <- ci[parm, , , drop=FALSE]
  dimnames(ci)[[2L]] <- pct
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
    fits <- vapply(seq_len(dimX[3L]), function(x) object$X[, , x] %*% coeffs[, x],
                   numeric(dimX[1L]))
    dimnames(fits) <- list(case.names(object), dimnames(coeffs)[[2L]])
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
  if (!is.null(object$residuals)) return(object$residuals)
  cf <- object$coefficients <- coef(object)
  type <- match.arg(type)
  if (type == 'partial') {
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
    fits <- fitted(object)
    y <- if (dim(object$y)[2L] == 1L) c(object$y) else object$y
    resids <- y - fits
  }
  return(resids)
}

#' @export
#' @rdname glm_stats
deviance.bg_GLM <- function(object, ...) apply(residuals(object), 2L, crossprod)

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
  df <- object$df.residual
  if (is.null(df)) df <- -diff(dim(object$X)[1L:2L])
  return(df)
}

#' @export
#' @rdname glm_stats
sigma.bg_GLM <- function(object, ...) sqrt(deviance(object) / df.residual(object))

#' @export
#' @rdname glm_stats

vcov.bg_GLM <- function(object, ...) {
  X <- design2array(object)
  dimV <- dim(X)[c(2L, 2L, 3L)]
  namesV <- dimnames(X)[c(2L, 2L, 3L)]
  runX <- if (is.null(object$runX)) namesV[[3L]] else object$runX
  XtX <- array(NaN, dim=dimV, dimnames=namesV)
  XtX[, , runX] <- apply(X[, , runX, drop=FALSE], 3L, function(x) chol2inv(chol.default(crossprod(x))))
  sig <- sigma(object)^2
  XtX <- aperm(XtX, c(3L, 1L, 2L))
  aperm(sig * XtX, c(2L, 3L, 1L))
}

# Convenience function to calculate leave-one-out resid SD
sigma_loo <- function(resids, dfR) {
  var.hat <- apply(resids, 2L, function(x) c(crossprod(x)) - x^2)
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
  model$coefficients <- coef(model)
  model$residuals <- residuals(model)
  hat <- hatvalues(model)
  type <- match.arg(type)
  denom <- if (type == 'sd.1') t(sigma(model) * sqrt(1 - t(hat))) else 1 - hat
  model$residuals / denom
}

#' @export
#' @rdname glm_influence

rstudent.bg_GLM <- function(model, ...) {
  dfR <- df.residual(model)
  res <- residuals(model)
  hat <- hatvalues(model)
  var.hat <- sigma_loo(res, dfR)
  res / (var.hat * sqrt(1 - hat))
}

#' @export
#' @rdname glm_influence

hatvalues.bg_GLM <- function(model, ...) {
  if (!is.null(model$hatvalues)) return(model$hatvalues)
  X <- design2array(model)
  dimX <- dim(X)
  namesX <- dimnames(X)
  runX <- if (is.null(model$runX)) namesX[[3L]] else model$runX
  A <- U <- array(NaN, dim=dimX[c(2L, 2L, 3L)], dimnames=c(namesX[2L], namesX[2L], namesX[3L]))
  A[, , runX] <- apply(X[, , runX, drop=FALSE], 3L, crossprod)
  U[, , runX] <- apply(A[, , runX, drop=FALSE], 3L, chol.default)
  Z <- array(0, dim=dimX[c(2L, 1L, 3L)])
  for (i in seq_len(dimX[3L])) {
    Z[, , i] <- forwardsolve(t(U[, , i]), t(X[, , i]))
  }
  res <- apply(Z, 3L, function(x) colSums(x^2))
  dimnames(res) <- c(namesX[1L], namesX[3L])
  return(res)
}

#' @export
#' @rdname glm_influence

cooks.distance.bg_GLM <- function(model, ...) {
  X <- design2array(model)
  p <- qr(X[, , 1L])$rank
  hat <- hatvalues(model)
  res <- model$residuals <- residuals(model)
  sd <- sigma(model)
  hat * (res / t(sd * (1 - t(hat))))^2 / p
}

#' @export
#' @rdname glm_influence

dffits.bg_GLM <- function(model) {
  res <- residuals(model)
  hat <- hatvalues(model)
  sig <- sigma_loo(res, df.residual(model))
  res * sqrt(hat) / (sig * (1 - hat))
}

#' @export
#' @rdname glm_influence

dfbeta.bg_GLM <- function(model, ...) {
  Nobs <- nobs(model)
  coeffs <- coef(model)
  dimC <- dim(coeffs)
  dfb <- array(0, dim=c(dimC, Nobs), dimnames=c(dimnames(coeffs), dimnames(model$X)[1L]))
  for (i in seq_len(Nobs)) dfb[, , i] <- coeffs - coef(model[-i])
  aperm(dfb, c(3L, 1L, 2L))
}

#' @export
#' @rdname glm_influence

dfbetas.bg_GLM <- function(model, ...) {
  X <- design2array(model)
  dimX <- dim(X)
  namesX <- dimnames(X)

  runX <- if (is.null(model$runX)) namesX[[3L]] else model$runX
  QR <- apply(X, 3L, qr)
  xxi <- matrix(NaN, dimX[2L], dimX[3L], dimnames=c(namesX[2L], namesX[3L]))
  xxi[, runX] <- vapply(QR[runX], function(x) sqrt(diag(chol2inv(x$qr, x$rank))), numeric(dimX[2L]))
  sig <- sigma_loo(residuals(model), df.residual(model))
  dfb <- if (is.null(model$dfbeta)) dfbeta(model) else model$dfbeta

  dfbs <- array(0, dim=dimX, dimnames=dimnames(dfb))
  for (i in seq_len(dimX[3L])) {
    dfbs[, , i] <- dfb[, , i] / tcrossprod(sig[, i], xxi[, i])
  }
  return(dfbs)
}

#' @export
#' @rdname glm_influence

covratio.bg_GLM <- function(model) {
  dfR <- df.residual(model)
  omh <- 1 - hatvalues(model)
  res <- residuals(model)
  sig <- sigma_loo(res, dfR)
  e.star <- res / (sig * sqrt(omh))
  1 / (omh * ((e.star^2 + dfR - 1) / dfR)^dim(model$X)[2L])
}

#' @param do.coef Logical indicating whether to calculate \code{dfbeta}
#' @export
#' @rdname glm_influence
#' @importFrom abind abind

influence.bg_GLM <- function(model, do.coef=TRUE, ...) {
  n <- nobs(model)
  hat <- model$hatvalues <- hatvalues(model)
  model$residuals <- residuals(model)
  sig <- sigma_loo(model$residuals, df.residual(model))
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
  cols <- terms(mod)
  tlabels <- labels(mod)
  v <- vcov(mod)
  if (any(grepl('Intercept', tlabels))) {
    int <- grep('Intercept', tlabels)
    v <- v[-int, -int, , drop=FALSE]
    tlabels <- tlabels[-int]
    cols <- lapply(cols[-int], `-`, 1L)
  } else {
    warning('No intercept; VIFs may not be sensible.')
  }

  regions <- region.names(mod)
  n.terms <- length(cols)
  if (n.terms < 2L) stop('Model contains fewer than 2 terms!')
  R <- array(apply(v, 3L, cov2cor), dim=dim(v), dimnames=dimnames(v))
  detR <- apply(R, 3L, det)
  res <- array(0, dim=c(n.terms, 3L, dim(v)[3L]))
  dimnames(res) <- list(tlabels, c('GVIF', 'Df', 'GVIF^(1/(2*Df))'), regions)

  for (rgn in regions) {
    res[, 1L, rgn] <- vapply(cols, function(j)
                            det(as.matrix(R[j, j, rgn])) * det(as.matrix(R[-j, -j, rgn])),
                            numeric(1L))
  }
  res[, 1L, ] <- t(t(res[, 1L, ]) / detR)
  df <- lengths(cols)
  res[, 2L, ] <- df
  if (any(df != 1L)) res[, 3L, ] <- res[, 1L, ]^(1 / (2L * df))
  return(res)
}

#' Calculate ANOVA table for a bg_GLM object
#'
#' \code{anova} calculates ANOVA tables for a \code{bg_GLM} object. The tests
#' performed are so-called \emph{Type III} tests.
#'
#' In addition to the standard ANOVA statistics (sum of squares, mean squares,
#' degrees of freedom, F statistics, and P-values), the output tables include:
#' \eqn{\eta^2}, partial \eqn{\eta^2}, \eqn{\omega^2}, and partial
#' \eqn{\omega^2} as measures of \emph{effect size}.
#'
#' @param region Character vector indicating the region(s) to calculate ANOVA
#'   statistics for. Default: \code{NULL} (use all regions)
#' @export
#' @return \code{anova} returns a \emph{list} of tables of class \code{anova}
#' @rdname glm_stats

anova.bg_GLM <- function(object, region=NULL, ...) {
  dfR <- df.residual(object)
  regions <- if (is.null(region)) region.names(object) else region
  RSS <- deviance(object)[regions]

  cols <- terms(object)
  tlabels <- labels(object)
  SSt <- matrix(0, length(RSS), length(cols), dimnames=list(names(RSS), tlabels))
  for (i in seq_along(cols)) SSt[, i] <- deviance(object[, -cols[[i]]])[regions]
  ss <- cbind(SSt - RSS, RSS)
  SSTot <- rowSums(ss)
  df_terms <- lengths(cols)
  df <- c(df_terms, dfR)
  ms <- t(t(ss) / df)
  eta2 <- ss / SSTot
  eta2.part <- ss / (ss + RSS)
  MSE <- ms[, dim(ms)[2L]]
  Ndf <- nobs(object) - df
  omega2 <- (ss - (outer(MSE, df))) / (SSTot + MSE)
  omega2.part <- (ss - (outer(MSE, df))) / (ss + outer(MSE, Ndf))
  f <- ms / (RSS / dfR)
  p <- t(apply(f, 1L, pf, df_terms, dfR, lower.tail=FALSE))

  tab <- setNames(vector('list', length(regions)), regions)
  for (i in regions) {
    tab[[i]] <- data.table(`Sum Sq`=ss[i, ], `Mean Sq`=ms[i, ], Df=df,
                           `F value`=f[i, ], `eta^2`=eta2[i, ], `Partial eta^2`=eta2.part[i, ],
                           `omega^2`=omega2[i, ], `Partial omega^2`=omega2.part[i, ],
                           `Pr(>F)`=p[i, ])
    tab[[i]][.N, c('F value', 'Pr(>F)', 'eta^2', 'Partial eta^2', 'omega^2', 'Partial omega^2') := NA]
    setDF(tab[[i]], rownames=c(tlabels, 'Residuals'))
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
  X <- design2array(object)
  p <- qr(X[, , 1L])$rank
  N0 <- N <- nobs(object)

  if (isTRUE(REML)) {
    N <- N - p
    QR <- colSums(apply(X, 3L, function(x) log(abs(diag(qr(x)$qr)[1L:p]))))
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
