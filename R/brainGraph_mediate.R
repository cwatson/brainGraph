#' Mediation analysis with brain graph measures as mediator variables
#'
#' \code{brainGraph_mediate} performs simple mediation analyses in which a given
#' graph- or vertex-level measure (e.g., \emph{weighted global efficiency})
#' is the mediator \emph{M}. The outcome (or dependent/response) variable
#' \emph{Y} can be a neuropsychological measure (e.g., \emph{IQ}) or
#' can be a disease-specific metric (e.g., recovery time). The treatment
#' variable should be a \code{factor}.
#'
#' This code was adapted closely from \code{\link[mediation]{mediate}} in the
#' \code{mediation} package, and the procedure is exactly the same as theirs
#' (see the references listed below). So, if you use this function, please cite
#' their work.
#'
#' As of \code{brainGraph v2.0.0}, this function has been tested only for a
#' treatment (independent) variable \emph{X} being a 2-level factor (e.g.,
#' disease group, old vs. young, etc.).
#'
#' Allowing for treatment-mediator interaction (setting \code{int=TRUE})
#' currently will only work properly if the mediator is a continuous variable;
#' since the mediator is always a graph metric, this should always be the case.
#'
#' @param covars A data table containing covariates of interest. It must include
#'   columns for \emph{Study.ID}, the treatment variable, \code{covar.names},
#'   and the outcome variable.
#' @param mediator Character string; the name of the graph measure acting as
#'   the \emph{mediating} variable
#' @param treat Character string; the \emph{treatment} variable (e.g.,
#'   \emph{Group})
#' @param outcome Character string; the name of the outcome variable of interest
#'   (e.g., full-scale IQ, memory, etc.)
#' @param covar.names Character vector of the column names in \code{covars} to
#'   include in the models as pre-treatment covariates.
#' @param boot Logical indicating whether or not to perform bootstrapping
#'   (default: \code{TRUE})
#' @param boot.ci.type Character string; which type of CI's to calculate
#'   (default: \code{perc})
#' @param N Integer; the number of bootstrap samples to run (default:
#'   \code{1e3})
#' @param conf.level Numeric; the level of the CI's to calculate (default:
#'   \code{0.95} for the 2.5 and 97.5 percentiles)
#' @param control.value Value of \code{treat} to be used as the control
#'   condition (default: \code{0})
#' @param treat.value Value of \code{treat} to be used as the treatment
#'   condition (default: \code{1})
#' @param long Logical indicating whether or not to return all bootstrap samples
#'   (default: \code{TRUE})
#' @param int Logical indicating whether or not to include an interaction of the
#'   mediator and treatment (default: \code{FALSE})
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM_design}}
#'   (e.g., \code{binarize}) (unused in the \code{summary} method)
#' @inheritParams GLM
#' @export
#' @importFrom RcppEigen fastLmPure
#'
#' @return An object of class \code{bg_mediate} with elements:
#'   \item{level}{Either \code{graph} or \code{vertex}.}
#'   \item{removed.subs}{A character vector of Study.ID's removed due to
#'     incomplete data}
#'   \item{X.m, X.y}{Design matrices for the model with the mediator as the
#'     outcome variable (\code{X.m}) and for the model with the mediator as an
#'     additional predictor (\code{X.y})}
#'   \item{y.m, y.y}{Outomce variables for the associated design matrices above.
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
#' @name MediationAnalysis
#' @aliases brainGraph_mediate
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
#' med.EglobWt.FSIQ <- brainGraph_mediate(dt.G[threshold == thresholds[5]],
#'   covars.med, 'E.global.wt', 'Group', 'FSIQ', covar.names=c('age', 'gender'),
#'   boot=TRUE, N=1e4)
#' med.strength.FSIQ <-
#'   brainGraph_mediate(dt.V[threshold == thresholds[5] & region == 'lcACC'],
#'                      covars.med, 'strength', 'Group', 'FSIQ',
#'                      covar.names=c('age', 'gender'), N=1e3)
#' }

brainGraph_mediate <- function(g.list, covars, mediator, treat,
                               outcome, covar.names, level=c('graph', 'vertex'),
                               boot=TRUE, boot.ci.type=c('perc', 'bca'), N=1e3,
                               conf.level=0.95, control.value=0, treat.value=1,
                               long=TRUE, int=FALSE, ...) {

  Study.ID <- region <- treatintstr <- NULL
  stopifnot(all(c(treat, outcome, covar.names) %in% names(covars)))
  if (!'Study.ID' %in% names(covars)) covars$Study.ID <- as.character(seq_len(nrow(covars)))
  covars <- droplevels(covars[, c('Study.ID', treat, covar.names, outcome), with=F])
  incomp <- covars[!complete.cases(covars), Study.ID]
  covars <- covars[!Study.ID %in% incomp]
  setkey(covars, Study.ID)

  level <- match.arg(level)
  if (level == 'graph') {
    dt.graph <- data.table(Study.ID=vapply(g.list, graph_attr, character(1), 'name'),
                           graph=vapply(g.list, graph_attr, numeric(1), mediator))
  } else if (level == 'vertex') {
    y <- t(vapply(g.list, vertex_attr, numeric(vcount(g.list[[1]])), mediator))
    dt.graph <- data.table(Study.ID=vapply(g.list, graph_attr, character(1), 'name'))
    dt.graph <- cbind(dt.graph, y)
    setnames(dt.graph, 2:ncol(dt.graph), V(g.list[[1]])$name)
  }

  dt.graph <- dt.graph[!Study.ID %in% incomp]
  DT <- merge(covars, dt.graph, on=key(covars))
  DT[, eval(treat) := as.factor(get(treat))]
  DT.m <- melt(DT, id.vars=names(covars), variable.name='region', value.name=mediator)

  t.levels <- DT[, levels(get(treat))]
  if (treat.value %in% t.levels & control.value %in% t.levels) {
    cat.0 <- control.value
    cat.1 <- treat.value
  } else {
    cat.0 <- t.levels[1]
    cat.1 <- t.levels[2]
  }

  X.m <- brainGraph_GLM_design(DT[, c(treat, covar.names), with=F], ...)
  n <- nrow(X.m)
  y.y <- DT[, get(outcome)]
  treatstr <- paste0(treat, cat.1)
  if (isTRUE(int)) treatintstr <- paste0(treatstr, ':', mediator)

  # Different across regions
  regions <- DT.m[, levels(region)]
  X.y <- res_boot <- sapply(regions, function(x) NULL)
  y.m <- matrix(0, nrow=n, ncol=length(regions), dimnames=list(DT.m[region == regions[1], Study.ID], regions))
  for (i in regions) {
    y.m[, i] <- DT.m[region == i, get(mediator)]
    if (isTRUE(int)) {
      X.y[[i]] <- brainGraph_GLM_design(DT.m[region == i, c(mediator, treat, covar.names), with=F],
                                        int=c(treat, mediator), ...)
    } else {
      X.y[[i]] <- brainGraph_GLM_design(DT.m[region == i, c(mediator, treat, covar.names), with=F], ...)
    }
    res_boot[[i]] <- boot_mediate(N, n, X.m, y.m[, i], treat, cat.1, X.y[[i]],
                                  y.y, mediator, treatstr, int, treatintstr)
  }
  res_boot <- rbindlist(res_boot, idcol=TRUE)
  setnames(res_boot, '.id', 'region')
  res_obs <- res_boot[, .SD[.N], by=region]
  res_p <- res_boot[, lapply(.SD, function(x) pval(x[1:N], x[N + 1])), by=region]
  res_boot <- res_boot[, .SD[-.N], by=region]

  low <- (1 - conf.level) / 2
  high <- 1 - low
  boot.ci.type <- match.arg(boot.ci.type)
  if (isTRUE(boot) & boot.ci.type == 'perc') {
    res_ci <- res_boot[, lapply(.SD, quantile, c(low, high), na.rm=TRUE), by=region]
  } else {
    res_ci <- res_boot[, lapply(.SD, BC.CI, low, high), by=region]
  }

  if (!isTRUE(long)) res_boot <- NULL

  out <- list(level=level, removed.subs=incomp, X.m=X.m, X.y=X.y, y.m=y.m, y.y=y.y,
              res.obs=res_obs, res.ci=res_ci, res.p=res_p,
              boot=boot, boot.ci.type=boot.ci.type, res.boot=res_boot,
              treat=treat, mediator=mediator, outcome=outcome, covariates=NULL, INT=int,
              conf.level=conf.level, control.value=cat.0, treat.value=cat.1,
              nobs=n, sims=N, covar.names=covar.names)
  out$atlas <- if (is.null(g.list[1]$atlas)) guess_atlas(g.list[1]) else g.list[1]$atlas
  class(out) <- c('bg_mediate', class(out))
  return(out)
}

boot_mediate <- function(N, n, X.m, y.m, treat, cat.1, X.y,
                         y.y, mediator, treatstr, int, treatintstr) {
  b <- tau <- d1 <- d0 <- z1 <- z0 <- n0 <- n1 <- d.avg <- z.avg <- n.avg <- NULL

  # Randomization/resampling matrix
  A <- matrix(rep(1:n, N), byrow=TRUE, nrow=N)
  index <- t(apply(A, 1, sample, replace=TRUE))
  index <- rbind(index, 1:n)

  # Loop through the resamples
  res <- foreach(b=seq_len(N + 1), .combine='rbind') %dopar% {

    # Mediator predictions
    est.m <- fastLmPure(X.m[index[b, ], ], y.m[index[b, ]], method=2)
    error <- rnorm(n, mean=0, sd=est.m$s)

    X.m.t <- X.m.c <- X.m
    X.m.t[, treatstr] <- 1
    X.m.c[, treatstr] <- 0
    PredictM1 <- X.m.t %*% est.m$coefficients + error
    PredictM0 <- X.m.c %*% est.m$coefficients + error

    # Outcome predictions
    est.y <- fastLmPure(X.y[index[b, ], ], y.y[index[b, ]], method=2)
    effects.tmp <- matrix(NA, nrow=n, ncol=4)
    for (e in 1:4) { # These calculate d1, d0, z1, z0 (respectively) for each "sim")
      tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0), c(1, 0, 1, 1), c(1, 0, 0, 0))
      X.y.t <- X.y.c <- X.y
      #TODO: my replacement of the following 2 lines prob. won't work for #groups > 2
#      cat.t <- ifelse(tt[1], cat.1, cat.0)
#      cat.c <- ifelse(tt[2], cat.1, cat.0)
#      pred.data.t[, treat] <- factor(cat.t, levels=t.levels)
#      pred.data.c[, treat] <- factor(cat.c, levels=t.levels)
      X.y.t[, treatstr] <- tt[1]
      X.y.c[, treatstr] <- tt[2]
      X.y.t[, mediator] <- PredictM1 * tt[3] + PredictM0 * (1 - tt[3])  #PredictMt
      X.y.c[, mediator] <- PredictM1 * tt[4] + PredictM0 * (1 - tt[4])  #PredictMc
      if (isTRUE(int)) {
        X.y.t[, treatintstr] <- X.y.t[, treatstr] * X.y.t[, mediator]
        X.y.c[, treatintstr] <- X.y.c[, treatstr] * X.y.c[, mediator]
      }
      pr.1 <- X.y.t %*% est.y$coefficients
      pr.0 <- X.y.c %*% est.y$coefficients
      pr.mat <- as.matrix(cbind(pr.1, pr.0))
      effects.tmp[, e] <- pr.mat[, 1] - pr.mat[, 2]
    }

    return(t(colMeans(effects.tmp)))
  }
  res <- as.data.table(res)
  setnames(res, c('d1', 'd0', 'z1', 'z0'))
  res[, tau := (d1 + d0 + z1 + z0) / 2]
  res[, n0 := d0 / tau]
  res[, n1 := d1 / tau]
  res[, d.avg := (d1 + d0) / 2]
  res[, z.avg := (z1 + z0) / 2]
  res[, n.avg := (n1 + n0) / 2]
  return(res)
}

pval <- function(x, xhat) {
  out <- ifelse(xhat == 0, 1, 2 * min(sum(x > 0), sum(x < 0)) / length(x))
  return(min(out, 1))
}

BC.CI <- function(theta, low, high) {
  z.inv <- length(theta[theta < mean(theta)]) / length(theta)
  z <- qnorm(z.inv)
  U <- (length(theta) - 1) * (mean(theta) - theta)
  top <- sum(U^3)
  under <- 6 * (sum(U^2))^{3/2}
  a <- top / under
  lower.inv <-  pnorm(z + (z + qnorm(low)) / (1 - a * (z + qnorm(low))))
  lower2 <- lower <- quantile(theta, lower.inv)
  upper.inv <-  pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper2 <- upper <- quantile(theta, upper.inv)
  return(c(lower, upper))
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
#' @method summary bg_mediate
#' @rdname mediation

summary.bg_mediate <- function(object, mediate=FALSE, region=NULL, digits=max(3L, getOption('digits') - 2L), ...) {
  stopifnot(inherits(object, 'bg_mediate'))
  ci.low0.acme <- d0.ci <- ci.high0.acme <- ci.low0.ade <- z0.ci <- ci.high0.ade <- ci.low.tot <-
    tau.ci <- ci.high.tot <- ci.low0.prop <- n0.ci <- ci.high0.prop <- ci.low1.acme <- d1.ci <-
    ci.high1.acme <- ci.low1.ade <- z1.ci <- ci.high1.ade <- ci.low1.prop <- n1.ci <- ci.high1.prop <-
    ci.low.avg.acme <- d.avg.ci <- ci.high.avg.acme <- ci.low.avg.ade <- z.avg.ci <- ci.high.avg.ade <-
    ci.low.avg.prop <- n.avg.ci <- ci.high.avg.prop <- Mediator <- treat <- Outcome <- NULL

  DT.obs <- copy(object$res.obs)
  DT.ci <- copy(object$res.ci)
  DT.p <- copy(object$res.p)
  setnames(DT.ci, c('region', paste0(names(DT.ci)[-1], '.ci')))
  setnames(DT.p, c('region', paste0(names(DT.p)[-1], '.p')))
  DT.all <- merge(merge(DT.obs, DT.p, by='region'), DT.ci, by='region')
  setnames(DT.all,
           c('d0', 'd0.p', 'z0', 'z0.p', 'tau', 'tau.p', 'n0', 'n0.p'),
           c('b0.acme', 'p0.acme', 'b0.ade', 'p0.ade', 'b.tot', 'p.tot', 'b0.prop', 'p0.prop'))
  DT.all[, ci.low0.acme := .SD[1, d0.ci], by=region]
  DT.all[, ci.high0.acme := .SD[2, d0.ci], by=region]
  DT.all[, ci.low0.ade := .SD[1, z0.ci], by=region]
  DT.all[, ci.high0.ade := .SD[2, z0.ci], by=region]
  DT.all[, ci.low0.prop := .SD[1, n0.ci], by=region]
  DT.all[, ci.high0.prop := .SD[2, n0.ci], by=region]
  DT.all[, ci.low.tot := .SD[1, tau.ci], by=region]
  DT.all[, ci.high.tot := .SD[2, tau.ci], by=region]
  DT.all[, Mediator := object$mediator]
  DT.all[, treat := object$treat]
  DT.all[, Outcome := object$outcome]

  # Different behavior if mediator-treatment interaction
  if (isTRUE(object$INT)) {
    setnames(DT.all,
             c('d1', 'd1.p', 'd.avg', 'd.avg.p', 'z1', 'z1.p', 'z.avg', 'z.avg.p', 'n1', 'n1.p', 'n.avg', 'n.avg.p'),
             c('b1.acme', 'p1.acme', 'b.avg.acme', 'p.avg.acme', 'b1.ade', 'p1.ade', 'b.avg.ade', 'p.avg.ade', 'b1.prop', 'p1.prop', 'b.avg.prop', 'p.avg.prop'))
    DT.all[, ci.low1.acme := .SD[1, d1.ci], by=region]
    DT.all[, ci.high1.acme := .SD[2, d1.ci], by=region]
    DT.all[, ci.low1.ade := .SD[1, z1.ci], by=region]
    DT.all[, ci.high1.ade := .SD[2, z1.ci], by=region]
    DT.all[, ci.low1.prop := .SD[1, n1.ci], by=region]
    DT.all[, ci.high1.prop := .SD[2, n1.ci], by=region]
    DT.all[, ci.low.avg.acme := .SD[1, d1.ci], by=region]
    DT.all[, ci.high.avg.acme := .SD[2, d1.ci], by=region]
    DT.all[, ci.low.avg.ade := .SD[1, z1.ci], by=region]
    DT.all[, ci.high.avg.ade := .SD[2, z1.ci], by=region]
    DT.all[, ci.low.avg.prop := .SD[1, n1.ci], by=region]
    DT.all[, ci.high.avg.prop := .SD[2, n1.ci], by=region]
    DT.all[, grep('.*.ci', names(DT.all)) := NULL]
    setcolorder(DT.all, c('Mediator', 'treat', 'Outcome', 'region',
                          'b0.acme', 'ci.low0.acme', 'ci.high0.acme', 'p0.acme',
                          'b1.acme', 'ci.low1.acme', 'ci.high1.acme', 'p1.acme',
                          'b.avg.acme', 'ci.low.avg.acme', 'ci.high.avg.acme', 'p.avg.acme',
                          'b0.ade', 'ci.low0.ade', 'ci.high0.ade', 'p0.ade',
                          'b1.ade', 'ci.low1.ade', 'ci.high1.ade', 'p1.ade',
                          'b.avg.ade', 'ci.low.avg.ade', 'ci.high.avg.ade', 'p.avg.ade',
                          'b.tot', 'ci.low.tot', 'ci.high.tot', 'p.tot',
                          'b0.prop', 'ci.low0.prop', 'ci.high0.prop', 'p0.prop',
                          'b1.prop', 'ci.low1.prop', 'ci.high1.prop', 'p1.prop',
                          'b.avg.prop', 'ci.low.avg.prop', 'ci.high.avg.prop', 'p.avg.prop'))

  } else {
    DT.all[, grep('1', names(DT.all)) := NULL]
    DT.all[, grep('avg', names(DT.all)) := NULL]
    DT.all[, grep('.*.ci', names(DT.all)) := NULL]
    setcolorder(DT.all, c('Mediator', 'treat', 'Outcome', 'region',
                          'b0.acme', 'ci.low0.acme', 'ci.high0.acme', 'p0.acme',
                          'b0.ade', 'ci.low0.ade', 'ci.high0.ade', 'p0.ade',
                          'b.tot', 'ci.low.tot', 'ci.high.tot', 'p.tot',
                          'b0.prop', 'ci.low0.prop', 'ci.high0.prop', 'p0.prop'))
  }
  DT.all <- DT.all[, .SD[1], keyby=region]
  DT.all[, region := as.factor(region)]

  object$DT.sum <- DT.all
  object$region <- region
  object$digits <- digits
  object$mediate <- mediate
  class(object) <- c('summary.bg_mediate', class(object))
  return(object)
}

#' @aliases summary.bg_mediate
#' @method print summary.bg_mediate
#' @keywords internal

print.summary.bg_mediate <- function(x, ...) {
  region <- NULL
  print_title_summary('brainGraph mediation results')
  cat('# of observations: ', x$nobs, '\n')
  cat('Level: ', x$level, '\n')
  cat('Mediator: ', x$mediator, '\n')
  cat('Treatment: ', x$treat, '\n')
  cat('\tControl value: ', x$control.value, '\n')
  cat('\tTreatment value: ', x$treat.value, '\n')
  cat('Outcome: ', x$outcome, '\n')
  cov.df <- data.frame(x$covar.names)
  dimnames(cov.df)[[2]] <- ''
  cat('Pre-treatment covariates:\n ')
  print(cov.df, right=FALSE)

  cat('\nTreatment-mediator interaction? ', x$INT, '\n\n')

  if (isTRUE(x$boot)) {
    low <- (1 - x$conf.level) / 2
    high <- 1 - low
    message('\n', 'Bootstrapping', '\n', rep('-', 30))
    ci <- switch(x$boot.ci.type,
                 perc='Percentile bootstrap',
                 bca='Bias-corrected accelerated')
    cat('Bootstrap CI type: ', ci, '\n')
    cat('# of bootstrap replicates: ', prettyNum(x$sims, ','), '\n')
    cat('Bootstrap CI level: ',
        sprintf('[%s]', paste(paste0(100 * c(low, high), '%'), collapse=' ')), '\n\n')
  }

  if (isTRUE(x$mediate)) {
    message('Mediation summary for: ', x$region, '\n', rep('-', 30))
    print(summary(bg_to_mediate(x, region=x$region)))
  } else {
    if (is.null(x$region)) {
      regions <- x$DT.sum[, levels(region)]
    } else {
      regions <- x$region
    }
    message('Mediation statistics', '\n', rep('-', 20))
    print(x$DT.sum[region %in% regions])
  }
  invisible(x)
}

#' Convert brainGraph results to "mediate" object
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
  if (!inherits(x, 'bg_mediate')) {
    stop('Use only with \'bg_mediate\' objects!')
  }
  if (x$level == 'graph') {
    res.obs <- x$res.obs
    res.ci <- x$res.ci
    res.p <- x$res.p
  } else {
    regions <- x$res.obs[, unique(region)]
    if (is.null(region)) {
      i <- 1
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
