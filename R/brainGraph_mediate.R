#' Mediation analysis with brain graph measures as mediator variables
#'
#' \code{brainGraph_mediate} performs simple mediation analyses in which a given
#' graph- or vertex-level measure (e.g., \emph{weighted global efficiency}) is
#' the mediator \emph{M}. The outcome (or dependent/response) variable \emph{Y}
#' can be a neuropsychological measure (e.g., \emph{IQ}) or can be a
#' disease-specific metric (e.g., recovery time).
#'
#' This code was adapted closely from \code{\link[mediation]{mediate}} in the
#' \code{mediation} package, and the procedure is exactly the same as theirs
#' (see the references listed below). If you use this function, please cite
#' their work.
#'
#' @note As of \code{brainGraph v2.0.0}, this function has been tested only for
#' a treatment (independent) variable \emph{X} being a \emph{factor} (e.g.,
#' disease group, old vs. young, etc.). If your treatment variable has more
#' than 2 levels, then you must explicitly specify the levels you would like to
#' compare; otherwise, the baseline and first levels are taken to be the
#' control and treatment values, respectively. Be aware that these are \emph{0}
#' indexed; that is, if you have 3 groups and you would like the treatment
#' group to be the 3rd, you should specify as either the group's character
#' string or as \code{treat.value=2}.
#'
#' @note Allowing for treatment-mediator interaction (setting \code{int=TRUE})
#' currently will only work properly if the mediator is a continuous variable;
#' since the mediator is always a graph metric, this should always be the case.
#'
#' @param covars A data table containing covariates of interest. It must include
#'   columns for \code{getOption('bg.subject_id')}, \code{treat},
#'   \code{outcome}, and \code{covar.names}.
#' @param mediator Character string; the name of the graph measure acting as
#'   the \emph{mediating} variable
#' @param treat Character string; the \emph{treatment} variable (e.g.,
#'   \emph{Group})
#' @param outcome Character string; the name of the outcome variable of interest
#' @param covar.names Character vector of the column name(s) in \code{covars} to
#'   include in the models as pre-treatment covariate(s).
#' @param control.value Value of \code{treat} to be used as the control
#'   condition. Default: \code{0}
#' @param treat.value Value of \code{treat} to be used as the treatment
#'   condition. Default: \code{1}
#' @param int Logical indicating whether or not to include an interaction of the
#'   mediator and treatment. Default: \code{FALSE}
#' @param boot Logical indicating whether or not to perform bootstrapping. This
#'   should always be done. Default: \code{TRUE}
#' @param boot.ci.type Character string; which type of CI's to calculate.
#'   Default: \code{perc}
#' @param N Integer; the number of bootstrap samples to run. Default:
#'   \code{1e3}
#' @param conf.level Numeric between 0 and 1; the level of the CI's to
#'   calculate. Default: \code{0.95} for the 2.5 and 97.5 percentiles)
#' @param long Logical indicating whether or not to return all bootstrap
#'   samples. Default: \code{FALSE}
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM_design}}
#'   (e.g., \code{binarize}) (unused in the \code{summary} method)
#' @inheritParams GLM
#' @export
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#'
#' @return An object of class \code{bg_mediate} with elements:
#'   \item{level}{Either \code{graph} or \code{vertex}.}
#'   \item{removed.subs}{A character vector of Study.ID's removed due to
#'     incomplete data}
#'   \item{X.m, X.y}{Design matrix and numeric array for the model with the
#'     mediator as the outcome variable (\code{X.m}) and for the model with the
#'     mediator as an additional predictor (\code{X.y}), respectively}
#'   \item{y.m, y.y}{Outcome variables for the associated design matrices above.
#'     \code{y.m} will be a matrix of size \emph{# subj. X # regions}}
#'   \item{res.obs}{A \code{data.table} of the observed values of the point
#'     estimates.}
#'   \item{res.ci}{A \code{data.table} of the confidence intervals for the
#'     effect estimates.}
#'   \item{res.p}{A \code{data.table} of the two-sided p-values for the effect
#'     estimates}
#'   \item{boot}{Logical, the \code{boot} argument.}
#'   \item{boot.ci.type}{Character string indicating which type of bootstrap
#'     confidence intervals were calculated.}
#'   \item{res.boot}{A \code{data.table} with \code{N} rows of the bootstrap
#'     results for all effects.}
#'   \item{treat}{Character string of the treatment variable.}
#'   \item{mediator}{Character string of the mediator variable.}
#'   \item{outcome}{Character string of the outcome variable.}
#'   \item{covariates}{Returns \code{NULL}; not used in this package.}
#'   \item{INT}{Logical indicating whether the models included an interaction
#'     between treatment and mediator.}
#'   \item{conf.level}{The confidence level.}
#'   \item{control.value}{The value of the treatment variable used as the
#'     control condition.}
#'   \item{treat.value}{The value of the treatment variable used as the
#'     treatment condition.}
#'   \item{nobs}{Integer; the number of observations in the models.}
#'   \item{sims}{Integer; the number of bootstrap replications.}
#'   \item{covar.names}{The pre-treatment covariate names.}
#'
#' @name Mediation
#' @rdname mediation
#' @family Group analysis functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Tingley, D. and Yamamoto, T. and Hirose, K. and Keele, L. and
#'   Imai, K. (2014) mediation: R package for causal mediation analysis.
#'   \emph{Journal of Statistical Software}, \bold{59(5)}, 1--38.
#'   \url{https://dx.doi.org/10.18637/jss.v059.i05}
#' @references Imai, K. and Keele, L. and Yamamoto, T. (2010) Identification
#'   inference, and sensitivity analysis for causal mediation effects.
#'   \emph{Statistical Science}, \bold{25(1)}, 51--71.
#'   \url{https://dx.doi.org/10.1214/10-STS321}
#' @references Imai, K. and Keele, L. and Tingley, D. (2010) A general approach
#'   to causal mediation analysis. \emph{Psychological Methods}, \bold{15(4)},
#'   309--334. \url{https://dx.doi.org/10.1037/a0020761}
#' @references Imai, K. and Keele, L. and Tingley, D. and Yamamoto, T. (2011)
#'   Unpacking the black box of causality: learning about causal mechanisms from
#'   experimental and observational studies. \emph{American Political Science
#'   Review}, \bold{105(4)}, 765--789.
#'   \url{https://dx.doi.org/10.1017/S0003055411000414}
#' @references Imai, K. and Yamamoto, T. (2013) Identification and sensitivity
#'   analysis for multiple causal mechanisms: revisiting evidence from framing
#'   experiments. \emph{Political Analysis}, \bold{21(2)}, 141--171.
#'   \url{https://dx.doi.org/10.1093/pan/mps040}
#' @examples
#' \dontrun{
#' med.EglobWt.FSIQ <- brainGraph_mediate(g[[5]], covars.med, 'E.global.wt',
#'   'Group', 'FSIQ', covar.names=c('age', 'gender'), N=1e4)
#' med.strength.FSIQ <- brainGraph_mediate(g[[5]], covars.med, 'strength',
#'   'Group', 'FSIQ', covar.names=c('age', 'gender'), level='vertex')
#' }

brainGraph_mediate <- function(g.list, covars, mediator, treat,
                               outcome, covar.names, level=c('graph', 'vertex'),
                               control.value=0, treat.value=1, int=FALSE,
                               boot=TRUE, boot.ci.type=c('perc', 'bca'), N=1e3,
                               conf.level=0.95, long=FALSE, ...) {
  region <- NULL
  if (!is.brainGraphList(g.list)) try(g.list <- as_brainGraphList(g.list))
  g.list <- g.list[]

  stopifnot(all(hasName(covars, c(treat, outcome, covar.names))))
  sID <- getOption('bg.subject_id')
  if (!hasName(covars, sID)) covars[, eval(sID) := seq_len(dim(covars)[1L])]
  covars[, eval(sID) := check_sID(get(sID))]
  covars <- droplevels(covars[, c(sID, treat, covar.names, outcome), with=FALSE])
  incomp <- covars[!complete.cases(covars), which=TRUE]
  names(incomp) <- covars[incomp, get(sID)]

  level <- match.arg(level)
  dt.graph <- glm_data_table(g.list, level, mediator)
  DT <- covars[dt.graph, on=sID]
  if (length(incomp) > 0L) DT <- DT[-incomp]

  DT[, eval(treat) := as.factor(get(treat))]
  t.levels <- DT[, levels(get(treat))]
  if (all(c(treat.value, control.value) %in% t.levels)) {
    cat.0 <- control.value
    cat.1 <- treat.value
  } else {
    stopifnot(is.numeric(c(control.value, treat.value)))
    cat.0 <- t.levels[control.value + 1L]
    cat.1 <- t.levels[treat.value + 1L]
  }

  cols <- c(sID, treat, covar.names)
  X.m <- brainGraph_GLM_design(DT[, cols, with=FALSE], ...)
  y.y <- DT[, get(outcome)]
  treatstr <- paste0(treat, cat.1)

  # Different across regions
  DT.m <- melt(DT, id.vars=names(covars), variable.name='region', value.name=mediator)
  regions <- names(dt.graph)[-1L]
  y.m <- as.matrix(DT[, c(sID, regions), with=FALSE], rownames=sID)
  des_args <- list(...)
  if (isTRUE(int)) des_args <- c(des_args, list(int=c(treat, mediator)))
  cols <- append(cols, mediator, after=1L)
  X.y <- lapply(regions, function(r)
                do.call(brainGraph_GLM_design,
                        c(list(covars=DT.m[region == r, cols, with=FALSE]), des_args)))
  names(X.y) <- regions
  attrs <- attributes(X.y[[1L]])[-c(1L, 2L)]
  X.y <- abind::abind(X.y, along=3L)
  attributes(X.y) <- c(attributes(X.y), attrs)

  # Calculate the resampled statistics for each region
  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  res_boot <- boot_mediate(X.m, y.m, X.y, y.y, mediator, treat, treatstr, int, N)

  res_p <- as.data.table(pvalArray(res_boot, N, dim(y.m)[2L]), keep.rownames='region')
  res_boot <- rbindlist(apply(res_boot, 3L, as.data.table), idcol='region')
  res_obs <- res_boot[, .SD[.N], by=region]
  res_boot <- res_boot[, .SD[-.N], by=region]

  conf.limits <- (1 + c(-1, 1) * conf.level) / 2
  if (isTRUE(boot)) {
    boot.ci.type <- match.arg(boot.ci.type)
    bootFun <- switch(boot.ci.type, perc=fastquant, bca=BC.CI2)
    if (boot.ci.type == 'bca') conf.limits <- qnorm(conf.limits)
    res_ci <- res_boot[, lapply(.SD, bootFun, conf.limits, N), by=region]
  }

  if (isFALSE(long)) res_boot <- NULL
  out <- list(level=level, removed.subs=incomp, X.m=X.m, X.y=X.y, y.m=y.m, y.y=y.y,
              res.obs=res_obs, res.ci=res_ci, res.p=res_p,
              boot=boot, boot.ci.type=boot.ci.type, res.boot=res_boot,
              treat=treat, mediator=mediator, outcome=outcome, covariates=NULL, INT=int,
              conf.level=conf.level, control.value=cat.0, treat.value=cat.1,
              nobs=dim(X.m)[1L], sims=N, covar.names=covar.names)
  out$atlas <- guess_atlas(g.list[[1L]])
  class(out) <- c('bg_mediate', class(out))
  return(out)
}

#' Calculate bootstrapped estimates of direct, indirect, and total effects
#'
#' @param X.m,X.y Design matrices for the mediation variable as outcome
#'   (\code{X.m}) and the mediation as a covariate predicting the specified
#'   outcome variable (\code{X.y})
#' @param y.m,y.y Numeric matrix and vector, respectively, containing the
#'   outcome variable data of the above designs
#' @param treatstr Character string denoting the column of the design matrix for
#'   the treatment variable
#' @inheritParams brainGraph_mediate
#' @noRd

boot_mediate <- function(X.m, y.m, X.y, y.y, mediator, treat, treatstr, int, N) {
  b <- NULL

  regions <- dimnames(y.m)[[2L]]
  dimXm <- dim(X.m)
  n <- dimXm[1L]
  pm <- dimXm[2L]
  dfRm <- n - pm
  py <- dim(X.y)[2L]
  ny <- dim(y.m)[2L]

  # Randomization/resampling matrix
  index <- t(replicate(N, sample.int(n, replace=TRUE)))
  index <- rbind(index, seq_len(n))

  # Don't use "drop" when it isn't needed (i.e., vertex-level). It is much slower.
  if (ny == 1L) {
    Xyperm <- function(X, porder) X[porder, , , drop=FALSE]
    Xyfun <- f_beta_3d_g
  } else {
    Xyperm <- function(X, porder) X[porder, , ]
    Xyfun <- f_beta_3d
  }
  # For faster calculation of the "Q" matrix for both "fits.m" and "beta.y"
  diagIndsY <- diag(1, n, py)
  diagIndsM <- diag(1, n, pm)

  # These are used in the loop calculating "effects.tmp"
  # The last 2 columns are (1 - tt[3]) and (1 - tt[4])
  ttMat <- matrix(c(1, 1, 1, 0, 0, 1,
                    0, 0, 1, 0, 0, 1,
                    1, 0, 1, 1, 0, 0,
                    1, 0, 0, 0, 1, 1),
                   nrow=6L)
  effects.tmp <- array(NA, dim=c(n, 4L, ny))

  vnames <- dimnames(X.y)[[2L]]
  treatstrOther <- setdiff(grep(treat, dimnames(X.m)[[2L]], value=TRUE), treatstr)
  Xcols <- c(mediator, treatstr)
  if (isTRUE(int)) {
    treatintstr <- paste0(treatstr, ':', mediator)
    Xcols <- c(Xcols, treatintstr)
    ecols <- 1L:4L
  } else {
    ecols <- c(1L, 3L)  # "d1" == "d0", and "z1" == "z0", so don't recalculate
  }
  bcols <- which(vnames %in% Xcols)

  # When "int=FALSE", reduces to the "Baron & Kenny" approach (for "d1", at least)
  fun_effects <- function(Xy.diff, beta.y, Xcols, bcols, regions) {
    X <- Xy.diff[, Xcols, ]
    b <- beta.y[bcols, ]
    vapply(regions, function(r) X[, , r] %*% b[, r], numeric(n))
  }

  # Loop through the resamples
  #-----------------------------------------------------------------------------
  res <- foreach(b=seq_len(N + 1L), .combine=rbind) %dopar% {
    neworder <- index[b, ]

    # Mediator predictions
    fits.m <- f_beta_m(X.m[neworder, ], y.m[neworder, ], diagIndsM, n, pm, ny, dfRm)
    error <- vapply(fits.m$sigma, function(r) rnorm(n, mean=0, sd=r), numeric(n))

    X.m.t <- X.m.c <- X.m
    X.m.t[, treatstr] <- 1
    X.m.c[, treatstr] <- 0
    X.m.t[, treatstrOther] <- X.m.c[, treatstrOther] <- 0
    PredictM1 <- X.m.t %*% fits.m$coefficients + error
    PredictM0 <- X.m.c %*% fits.m$coefficients + error

    # Outcome predictions
    beta.y <- Xyfun(Xyperm(X.y, neworder), y.y[neworder], regions, diagIndsY, n, py, ny)

    # e1: Mediation(1); e2: Mediation(0); e3: Direct(1); e4: Direct(0)
    # AKA: d1, d0, z1, z0 (respectively) for each "sim")
    for (e in ecols) {
      X.y.t <- X.y.c <- X.y
      X.y.t[, treatstr, ] <- ttMat[1L, e]
      X.y.c[, treatstr, ] <- ttMat[2L, e]
      X.y.t[, mediator, ] <- PredictM1 * ttMat[3L, e] + PredictM0 * ttMat[5L, e]  #PredictMt
      X.y.c[, mediator, ] <- PredictM1 * ttMat[4L, e] + PredictM0 * ttMat[6L, e]  #PredictMc
      if (isTRUE(int)) {
        #X.y.t[, treatstrOther, ] <- X.y.c[, treatstrOther, ] <- 0
        X.y.t[, treatintstr, ] <- X.y.t[, treatstr, ] * X.y.t[, mediator, ]
        X.y.c[, treatintstr, ] <- X.y.c[, treatstr, ] * X.y.c[, mediator, ]
      }
      Xy.diff <- X.y.t - X.y.c
      effects.tmp[, e, ] <- fun_effects(Xy.diff, beta.y, Xcols, bcols, regions)
    }

    return(colMeans(effects.tmp))
  }
  res <- array(res, dim=c(4L, N + 1L, ny))
  if (isFALSE(int)) res[c(2L, 4L), , ] <- res[ecols, , ]
  res2 <- array(dim=c(10L, N + 1L, ny))
  res2[1L:4L, , ] <- res
  res2 <- aperm(res2, c(2L, 1L, 3L))
  dimnames(res2)[2L:3L] <- list(c('d1', 'd0', 'z1', 'z0', 'tau', 'n0', 'n1', 'd.avg', 'z.avg', 'n.avg'), regions)
  res2[, 'd.avg', ] <- (res2[, 'd1', ] + res2[, 'd0', ]) / 2
  res2[, 'z.avg', ] <- (res2[, 'z1', ] + res2[, 'z0', ]) / 2
  res2[, 'tau', ] <- res2[, 'd.avg', ] + res2[, 'z.avg', ]    # i.e., (d1 + d0 + z1 + z0) / 2
  res2[, 'n0', ] <- res2[, 'd0', ] / res2[, 'tau', ]
  res2[, 'n1', ] <- res2[, 'd1', ] / res2[, 'tau', ]
  res2[, 'n.avg', ] <- (res2[, 'n1', ] + res2[, 'n0', ]) / 2
  return(res2)
}

# Only calculate coefficients and sigma for "y.m ~ X.m"
f_beta_m <- function(X, Y, diagMat, n, p, ny, dfR) {
  QR <- qr.default(X, LAPACK=TRUE)
  Q <- qr_Q2(QR, diagMat, n, p)
  R <- qr_R2(QR, p)
  beta <- backsolve(R, crossprod(Q, Y), p)
  beta[QR$pivot, ] <- beta
  ehat <- Y - X %*% beta
  s <- if (ny == 1L) sum(ehat^2) else .colSums(ehat^2, n, ny)
  list(coefficients=beta, sigma=sqrt(s / dfR))
}

# Functions to calculate coefficients for "y.y ~ X.y"; i.e., multiple designs, 1 outcome
f_beta_3d <- function(X, Y, regions, diagMat, n, p, ny) {
  QR <- qr(X, LAPACK=TRUE)
  Q <- lapply(QR, qr_Q2, diagMat, n, p)
  R <- lapply(QR, qr_R2, p)
  beta <- matrix(NaN, p, ny, dimnames=list(NULL, regions))
  for (r in regions) {
    beta[QR[[r]]$pivot, r] <- backsolve(R[[r]], crossprod(Q[[r]], Y), p)
  }
  beta
}

# If "level='graph'", this avoids using "drop=FALSE" above
f_beta_3d_g <- function(X, Y, regions, diagMat, n, p, ny=1L) {
  QR <- qr.default(X[, , 1L], LAPACK=TRUE)
  Q <- qr_Q2(QR, diagMat, n, p)
  R <- qr_R2(QR, p)
  beta <- backsolve(R, crossprod(Q, Y), p)
  beta[QR$pivot, ] <- beta
  dimnames(beta) <- list(NULL, 'graph')
  beta
}

#' Calculate P-values in an array for mediation analysis
#'
#' @param res_boot Numeric array with \eqn{(N+1) \times 10 \times N_y} dimensions,
#'   where \eqn{N} is the number of resamples and \eqn{N_y} is the number of regions
#' @noRd

pvalArray <- function(res_boot, N=dim(res_boot)[1L] - 1L, ny=dim(res_boot)[3L]) {
  seqN <- seq_len(N)
  gt0 <- colSums(res_boot[seqN, , ] > 0)
  lt0 <- N - gt0  # Only differs if there are "x == 0" exactly
  #lt0 <- colSums(res_boot[seqN, , ] < 0)
  pMat <- 2 * pmin.int(gt0, lt0) / N
  zeros <- which(res_boot[N + 1L, , ] == 0)
  if (length(zeros) > 0L) pMat[zeros] <- 1
  dim(pMat) <- c(10L, ny)
  dimnames(pMat) <- dimnames(res_boot)[2L:3L]
  t(pMat)
}

#' Bias-corrected and accelerated confidence intervals
#'
#' @param theta Numeric vector with the bootstrap-resampled statistics
#' @param quants Numeric vector with 2 elements: the confidence limits; i.e.,
#'   \code{qnorm(c((1 - alpha) / 2, (1 + alpha) / 2))}
#' @param N Integer; the number of resamples
#' @noRd

BC.CI2 <- function(theta, quants, N) {
  avg <- sum(theta) / N
  z.inv <- sum(theta < avg) / N
  z <- qnorm(z.inv)
  U <- (N - 1L) * (avg - theta)
  U2 <- U^2
  top <- sum(U * U2)
  under <- 6 * (sum(U2))^1.5
  a <- top / under
  lower.upper <-  pnorm(z + (z + quants) / (1 - a * (z + quants)))
  fastquant(theta, lower.upper, N)
}

#' A pared down "quantile" function with no argument checking
#'
#' This function uses the default \code{type=7} from
#' \code{\link[stats]{quantile}}.
#'
#' @param x Numeric vector of bootstrap-resampled statistics
#' @param probs Numeric vector of probabilities
#' @param N Integer; the number of resamples
#' @noRd

fastquant <- function(x, probs, N) {
  index <- 1 + (N - 1) * probs
  lo <- floor(index)
  hi <- ceiling(index)
  x <- sort(x, partial=unique(c(lo, hi)))
  qs <- x[lo]
  i <- which(index > lo & x[hi] != qs)
  h <- (index - lo)[i]
  qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
  qs
}

################################################################################
# S3 METHODS FOR "bg_mediate"
################################################################################

#' Print a summary from a brainGraph mediation analysis
#'
#' @param object A \code{bg_mediate} object
#' @param mediate Logical indicating whether or not to use the \code{summary}
#'   method from \code{\link[mediation]{mediate}} (default: \code{FALSE}). If
#'   \code{TRUE}, only a single region can be printed.
#' @param region Character string specifying which region's results to
#'   summarize; only relevant if \code{level='vertex'} (default: \code{NULL})
#' @export
#' @rdname mediation

summary.bg_mediate <- function(object, mediate=FALSE, region=NULL, digits=max(3L, getOption('digits') - 2L), ...) {
  stopifnot(inherits(object, 'bg_mediate'))
  Mediator <- treat <- Outcome <- NULL

  DT.obs <- copy(object$res.obs)
  DT.ci <- copy(object$res.ci)
  DT.p <- copy(object$res.p)
  setnames(DT.ci, c('region', paste0(names(DT.ci)[-1L], '.ci')))
  setnames(DT.p, c('region', paste0(names(DT.p)[-1L], '.p')))
  DT.all <- merge(merge(DT.obs, DT.p, by='region'), DT.ci, by='region')
  DT.all[, c('Mediator', 'treat', 'Outcome') := with(object, mediator, treat, outcome)]

  change <- matrix(c('d0', 'd0.p', 'z0', 'z0.p', 'tau', 'tau.p', 'n0', 'n0.p',
                     'b0.acme', 'p0.acme', 'b0.ade', 'p0.ade', 'b.tot', 'p.tot', 'b0.prop', 'p0.prop'),
                   ncol=2L)
  setnames(DT.all, change[, 1L], change[, 2L])
  change_ci <- change[seq.int(1L, 7L, 2L), ]
  change_ci <- cbind(paste0(change_ci[, 1L], '.ci'),
                     sub('[bdnz]([01]?)\\.', 'ci.low\\1.', change_ci[, 2L]))
  change_ci <- cbind(change_ci, sub('low', 'high', change_ci[, 2L]))
  DT.all[, eval(change_ci[, 2L]) := lapply(.SD, function(x) x[1L]), by=region, .SDcols=change_ci[, 1L]]
  DT.all[, eval(change_ci[, 3L]) := lapply(.SD, function(x) x[2L]), by=region, .SDcols=change_ci[, 1L]]

  mainnames <- c('Mediator', 'treat', 'Outcome', 'region')
  acme <- c('b0.acme', 'ci.low0.acme', 'ci.high0.acme', 'p0.acme')
  total <- sub('0.acme', '.tot', acme)
  # Different behavior if mediator-treatment interaction
  if (isTRUE(object$INT)) {
    change1 <- sub('0', '1', change)[-c(5L, 6L), ]
    change1 <- rbind(change1, sub('1', '.avg', change1))
    setnames(DT.all, change1[, 1L], change1[, 2L])
    change_ci1 <- sub('0', '1', change_ci)[-3L, ]
    change_ci1 <- rbind(change_ci1, sub('1', '.avg', change_ci1))
    DT.all[, eval(change_ci1[, 2L]) := lapply(.SD, function(x) x[1L]), by=region, .SDcols=change_ci1[, 1L]]
    DT.all[, eval(change_ci1[, 3L]) := lapply(.SD, function(x) x[2L]), by=region, .SDcols=change_ci1[, 1L]]
    acme <- c(acme, sub('0', '1', acme), sub('0', '.avg', acme))

  } else {
    DT.all[, grep('1|avg', names(DT.all)) := NULL]
  }
  DT.all[, grep('.*.ci', names(DT.all)) := NULL]
  setcolorder(DT.all, c(mainnames, acme, sub('acme', 'ade', acme), total, sub('acme', 'prop', acme)))
  DT.all <- DT.all[, .SD[1L], keyby=region]
  DT.all[, region := as.factor(region)]

  object <- c(object, list(DT.sum=DT.all, region=region, digits=digits, mediate=mediate))
  class(object) <- c('summary.bg_mediate', class(object))
  return(object)
}

#' @aliases summary.bg_mediate
#' @method print summary.bg_mediate
#' @export

print.summary.bg_mediate <- function(x, ...) {
  region <- NULL
  width <- getOption('width') / 4
  dashes <- rep.int('-', width)
  print_title_summary(simpleCap(x$level), '-level mediation results')
  cat('# of observations: ', x$nobs, '\n')

  # Print a table of the model variables
  message('\n', 'Variables', '\n', dashes)
  df <- data.frame(A=c('  Mediator:', ' Treatment:', '    Control condition:',
                       '  Treatment condition:', '   Outcome:'),
                   B=c(x$mediator, x$treat, x$control.value, x$treat.value, x$outcome))
  cov.df <- data.frame(A=c('', 'Covariates:', rep.int('', length(x$covar.names) - 1L)),
                       B=c('', x$covar.names))
  df <- rbind(df, cov.df)
  dimnames(df)[[2L]] <- rep.int('', 2L)
  print(df, right=FALSE, row.names=FALSE)
  cat('\nTreatment-mediator interaction? ', x$INT, '\n\n')

  print_subs_summary(x)

  if (isTRUE(x$boot)) {
    low <- (1 - x$conf.level) / 2
    high <- 1 - low
    message('\n', 'Bootstrapping', '\n', dashes)
    ci <- switch(x$boot.ci.type,
                 perc='Percentile bootstrap',
                 bca='Bias-corrected accelerated')
    cat('Bootstrap CI type: ', ci, '\n')
    cat('# of bootstrap replicates: ', prettyNum(x$sims, ','), '\n')
    cat('Bootstrap CI level: ',
        sprintf('[%s]', paste(paste0(100 * c(low, high), '%'), collapse=' ')), '\n\n')
  }

  if (isTRUE(x$mediate)) {
    if (!requireNamespace('mediation', quietly=TRUE)) stop('Must install the "mediation" package.')
    region <- if (is.null(x$region)) x$DT.sum[, levels(region)[1L]] else x$region
    message('Mediation summary for: ', region, '\n', dashes)
    print(summary(bg_to_mediate(x, region)))
  } else {
    if (is.null(x$region)) {
      regions <- x$DT.sum[, levels(region)]
    } else {
      regions <- x$region
    }
    message('Mediation statistics', '\n', dashes)
    print(x$DT.sum[region %in% regions])
  }
  invisible(x)
}

#' Convert brainGraph results to mediate object
#'
#' \code{\link{bg_to_mediate}} converts the results into an object of class
#' \code{\link[mediation]{mediate}}. In \code{brainGraph}, it is only used for
#' the \code{\link[mediation]{summary.mediate}} method, but you can similarly
#' use its output for the \code{\link[mediation]{plot.mediate}} method.
#'
#' @param x Object output from \code{\link{brainGraph_mediate}}
#' @export
#'
#' @return \code{bg_to_mediate} returns an object of class \code{mediate}
#' @seealso \code{\link[mediation]{mediate}}
#' @rdname mediation

bg_to_mediate <- function(x, region=NULL) {
  if (!inherits(x, c('bg_mediate', 'summary.bg_mediate'))) {
    stop('Use only with \'bg_mediate\' objects!')
  }
  if (x$level == 'graph') {
    res.obs <- x$res.obs
    res.ci <- x$res.ci
    res.p <- x$res.p
  } else {
    regions <- x$res.obs[, unique(region)]
    if (is.null(region)) {
      i <- 1L
    } else {
      i <- which(region == regions)
    }
    res.obs <- x$res.obs[region == regions[i]]
    res.ci <- x$res.ci[region == regions[i]]
    res.p <- x$res.p[region == regions[i]]
  }
  # My object is missing: call=cl, model.y=model.y, model.m=model.m)
  out <- list(d0=res.obs$d0, d1=res.obs$d1, d0.ci=res.ci$d0, d1.ci=res.ci$d1, d0.p=res.p$d0, d1.p=res.p$d1,
              z0=res.obs$z0, z1=res.obs$z1, z0.ci=res.ci$z0, z1.ci=res.ci$z1, z0.p=res.p$z0, z1.p=res.p$z1,
              n0=res.obs$n0, n1=res.obs$n1, n0.ci=res.ci$n0, n1.ci=res.ci$n1, n0.p=res.p$n0, n1.p=res.p$n1,
              tau.coef=res.obs$tau, tau.ci=res.ci$tau, tau.p=res.p$tau,
              d.avg=res.obs$d.avg, d.avg.ci=res.ci$d.avg, d.avg.p=res.p$d.avg,
              z.avg=res.obs$z.avg, z.avg.ci=res.ci$z.avg, z.avg.p=res.p$z.avg,
              n.avg=res.obs$n.avg, n.avg.ci=res.ci$n.avg, n.avg.p=res.p$n.avg,
              boot=x$boot, boot.ci.type=x$boot.ci.type, treat=x$treat, mediator=x$mediator, covariates=x$covariates,
              INT=x$INT, control.value=x$control.value, treat.value=x$treat.value, nobs=x$nobs, sims=x$sims,
              robustSE=FALSE, cluster=NULL)
  class(out) <- 'mediate'
  return(out)
}
