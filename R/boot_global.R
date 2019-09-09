#' Bootstrapping for global graph measures
#'
#' Perform bootstrapping to obtain groupwise standard error estimates of a
#' global graph measure (e.g. \emph{modularity}).
#'
#' The confidence intervals are calculated using the \emph{normal approximation}
#' at the \eqn{100 \times conf}\% level (by default, 95\%).
#'
#' For getting estimates of \emph{weighted global efficiency}, a method for
#' transforming edge weights must be provided. The default is to invert them.
#' See \code{\link{xfm.weights}}.
#'
#' @param densities Numeric vector of graph densities to loop through
#' @param resids An object of class \code{brainGraph_resids} (the output from
#'   \code{\link{get.resid}})
#' @param R Integer; the number of bootstrap replicates. Default: \code{1e3}
#' @param measure Character string of the measure to test. Default: \code{mod}
#' @param conf Numeric; the confidence level for calculating confidence
#'   intervals. Default: \code{0.95}
#' @param .progress Logical indicating whether or not to show a progress bar.
#'   Default: \code{TRUE}
#' @inheritParams xfm.weights
#' @export
#' @importFrom boot boot
#'
#' @return \code{brainGraph_boot} -- an object of class \code{brainGraph_boot}
#'   containing some input variables, in addition to a list of
#'   \code{\link[boot]{boot}} objects (one for each group).
#'
#' @name Bootstrapping
#' @family Group analysis functions
#' @family Structural covariance network functions
#' @seealso \code{\link[boot]{boot}}, \code{\link[boot]{boot.ci}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' boot.E.global <- brainGraph_boot(densities, resids.all, 1e3, 'E.global')
#' }

brainGraph_boot <- function(densities, resids, R=1e3,
                            measure=c('mod', 'E.global', 'Cp', 'Lp', 'assortativity',
                                      'strength', 'mod.wt', 'E.global.wt'),
                            conf=0.95, .progress=TRUE,
                            xfm.type=c('1/w', '-log(w)', '1-w', '-log10(w/max(w))', '-log10(w/max(w)+1)')) {
  stopifnot(inherits(resids, 'brainGraph_resids'))

  # 'statistic' function for the bootstrapping process
  statfun <- function(x, i, measure, res.obj, xfm.type) {
    corrs <- corr.matrix(res.obj[i], densities=densities, rand=TRUE)
    if (measure %in% c('strength', 'mod.wt', 'E.global.wt')) {
      g.boot <- apply(corrs[[1]]$r.thresh, 3, function(y)
                      graph_from_adjacency_matrix(corrs[[1]]$R * y, mode='undirected',
                                                  diag=FALSE, weighted=TRUE))
    } else {
      g.boot <- apply(corrs[[1]]$r.thresh, 3, graph_from_adjacency_matrix,
                      mode='undirected', diag=FALSE)
    }
    res <- switch(measure,
        mod=,mod.wt=vapply(g.boot, function(y) max(cluster_louvain(y)$modularity), numeric(1)),
        E.global=vapply(g.boot, efficiency, numeric(1), 'global'),
        E.global.wt=vapply(g.boot, function(g) efficiency(xfm.weights(g, xfm.type), 'global'), numeric(1)),
        Cp=vapply(g.boot, transitivity, numeric(1), type='localaverage'),
        Lp=vapply(g.boot, mean_distance, numeric(1)),
        assortativity=vapply(g.boot, assortativity_degree, numeric(1)),
        strength=vapply(g.boot, function(y) mean(strength(y)), numeric(1)))
    return(res)
  }

  # Show a progress bar so you aren't left in the dark
  if (isTRUE(.progress)) {
    intfun <- function(data, indices, measure, res.obj, xfm.type) {
      curVal <- get('counter', envir=env) + ncpus
      assign('counter', curVal, envir=env)
      setTxtProgressBar(get('progbar', envir=env), curVal)
      flush.console()
      statfun(data, indices, measure, res.obj, xfm.type)
    }
  } else {
    intfun <- statfun
  }

  if (.Platform$OS.type == 'windows') {
    my.parallel <- 'snow'
    ncpus <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
    cl <- makeCluster(ncpus, type='SOCK')
    clusterEvalQ(cl, library(brainGraph))
  } else {
    my.parallel <- 'multicore'
    ncpus <- detectCores()
    cl <- NULL
  }

  measure <- match.arg(measure)
  xfm.type <- match.arg(xfm.type)
  env <- environment()
  groups <- resids$groups
  my.boot <- setNames(vector('list', length(groups)), groups)
  for (g in groups) {
    counter <- 0
    res.dt <- resids$resids.all[g]
    if (isTRUE(.progress)) progbar <- txtProgressBar(min=0, max=R, style=3)
    my.boot[[g]] <- boot(res.dt, intfun, measure=measure, res.obj=resids[g], xfm.type=xfm.type, R=R,
                         parallel=my.parallel, ncpus=ncpus, cl=cl)
    if (isTRUE(.progress)) close(progbar)
  }

  out <- list(measure=measure, densities=densities, groups=groups, conf=conf, boot=my.boot)
  class(out) <- c('brainGraph_boot', class(out))
  return(out)
}

#' Print a summary from a bootstrap analysis
#'
#' @param object,x A \code{brainGraph_boot} object
#' @importFrom boot boot.ci
#' @export
#' @method summary brainGraph_boot
#' @rdname Bootstrapping

summary.brainGraph_boot <- function(object, ...) {
  kNumDensities <- length(object$densities)
  # Get everything into a data.table
  boot.dt <- with(object,
                  data.table(Group=rep(groups, each=kNumDensities),
                             density=rep(densities, length(groups)),
                             Observed=c(vapply(boot, with, numeric(kNumDensities), t0)),
                             se=c(vapply(boot, function(x) apply(x$t, 2, sd), numeric(kNumDensities)))))
  ci <- with(object,
             vapply(seq_along(densities), function(x)
                    vapply(boot, function(y)
                             boot.ci(y, type='norm', index=x, conf=conf)$normal[2:3],
                           numeric(2)),
                    numeric(2 * length(groups))))
  boot.dt$ci.low <- c(t(ci[2 * (seq_along(object$groups) - 1) + 1, ]))
  boot.dt$ci.high <- c(t(ci[2 * (seq_along(object$groups) - 1) + 2, ]))

  meas.full <- switch(object$measure,
                      mod='Modularity',
                      E.global='Global efficiency',
                      E.global.wt='Global efficiency (weighted)',
                      Cp='Clustering coefficient',
                      Lp='Average shortest path length',
                      assortativity='Degree assortativity',
                      strength='Average strength',
                      mod.wt='Modularity (weighted)')
  boot.sum <- list(meas.full=meas.full, DT.sum=boot.dt, conf=object$conf, R=object$boot[[1]]$R)
  class(boot.sum) <- c('summary.brainGraph_boot', class(boot.sum))
  boot.sum
}

#' @aliases summary.brainGraph_boot
#' @method print summary.brainGraph_boot

print.summary.brainGraph_boot <- function(x, ...) {
  print_title_summary('Bootstrap analysis')
  cat('Graph metric: ', x$meas.full, '\n')
  cat('Number of bootstrap samples generated: ', x$R, '\n')
  conf.pct <- 100 * x$conf
  cat(conf.pct, '% confidence intervals\n\n')
  print(x$DT.sum)
  invisible(x)
}

#' Plot bootstrap output of global graph measures across densities
#'
#' The \code{plot} method returns two \code{ggplot} objects: one with shaded
#' regions based on the standard error, and the other based on confidence
#' intervals (calculated using the normal approximation.
#'
#' @param ... Unused
#' @param alpha A numeric indicating the opacity for
#'   \code{\link[ggplot2]{geom_ribbon}}
#' @export
#' @method plot brainGraph_boot
#' @rdname Bootstrapping
#'
#' @return \code{plot} -- \emph{list} with the following elements:
#'   \item{se}{A ggplot object with ribbon representing standard error}
#'   \item{ci}{A ggplot object with ribbon representing confidence intervals}

plot.brainGraph_boot <- function(x, ..., alpha=0.4) {
  Observed <- Group <- se <- ci.low <- ci.high <- NULL

  boot.sum <- summary(x)
  boot.dt <- boot.sum$DT.sum
  b <- ggplot(boot.dt, aes(x=density, y=Observed, col=Group)) +
    geom_line() +
    labs(y=boot.sum$meas.full)

  se <- b + geom_ribbon(aes(ymin=Observed-se, ymax=Observed+se, fill=Group), alpha=alpha)
  ci <- b + geom_ribbon(aes(ymin=ci.low, ymax=ci.high, fill=Group), alpha=alpha)
  return(list(se=se, ci=ci))
}
