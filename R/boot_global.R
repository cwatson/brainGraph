#' Bootstrapping for global graph measures
#'
#' Performs bootstrapping to get group standard error estimates of a global
#' graph measure (e.g. modularity).
#'
#' The confidence intervals are calculated using the normal approximation and at
#' the 95\% level.
#'
#' @param densities Numeric vector of graph densities to loop through
#' @param resids A data.table of the residuals (from \code{\link{get.resid}})
#' @param R Integer; the number of bootstrap replicates (default: \code{1e3})
#' @param measure Character string of the measure to test (default: \code{mod})
#' @param conf Numeric; the confidence level for calculating confidence
#'   intervals (default: 0.95)
#' @param .progress Logical indicating whether or not to show a progress bar
#'   (default: \code{TRUE})
#' @export
#' @importFrom boot boot
#'
#' @return An object of class \code{brainGraph_boot} containing some input
#'   variables, in addition to a list of \code{\link[boot]{boot}} objects (one
#'   for each group).
#'
#' @family Group analysis functions
#' @family Structural covariance network functions
#' @seealso \code{\link[boot]{boot}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' boot.E.global <- brainGraph_boot(densities, resids.all, 1e3, 'E.global')
#' }

brainGraph_boot <- function(densities, resids, R=1e3,
                            measure=c('mod', 'E.global', 'Cp', 'Lp', 'assortativity',
                                      'strength', 'mod.wt', 'E.global.wt'),
                            conf=0.95, .progress=TRUE) {
  Group <- t0 <- NULL
  if (!is.factor(resids$Group)) resids[, Group := as.factor(Group)]
  groups <- resids[, levels(Group)]
  measure <- match.arg(measure)

  # 'statistic' function for the bootstrapping process
  statfun <- function(x, i, measure) {
    corrs <- corr.matrix(x[i], densities=densities)
    if (measure %in% c('strength', 'mod.wt', 'E.global.wt')) {
      g.boot <- apply(corrs$r.thresh, 3, function(y)
                      graph_from_adjacency_matrix(corrs$R * y, mode='undirected',
                                                  diag=FALSE, weighted=TRUE))
    } else {
      g.boot <- apply(corrs$r.thresh, 3, graph_from_adjacency_matrix,
                      mode='undirected', diag=FALSE)
    }
    res <- switch(measure,
        mod=,mod.wt=vapply(g.boot, function(x) max(modularity(cluster_louvain(x))), numeric(1)),
        E.global=,E.global.wt=vapply(g.boot, efficiency, numeric(1), 'global'),
        Cp=vapply(g.boot, transitivity, numeric(1), type='localaverage'),
        Lp=vapply(g.boot, mean_distance, numeric(1)),
        assortativity=vapply(g.boot, assortativity_degree, numeric(1)),
        strength=vapply(g.boot, function(x) mean(graph.strength(x)), numeric(1)))
    return(res)
  }

  # Show a progress bar so you aren't left in the dark
  if (isTRUE(.progress)) {
    intfun <- function(data, indices, measure) {
      curVal <- get('counter', envir=env) + ncpus
      assign('counter', curVal, envir=env)
      setTxtProgressBar(get('progbar', envir=env), curVal)
      flush.console()
      statfun(data, indices, measure)
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
    my.parallel='multicore'
    ncpus <- detectCores()
    cl <- NULL
  }
  env <- environment()
  my.boot <- sapply(groups, function(x) NULL)
  for (x in names(my.boot)) {
    counter <- 0
    if (isTRUE(.progress)) progbar <- txtProgressBar(min=0, max=R, style=3)
    my.boot[[x]] <- boot(resids[x], intfun, measure=measure, R=R,
                         parallel=my.parallel, ncpus=ncpus, cl=cl)
    if (isTRUE(.progress)) close(progbar)
  }

  out <- list(measure=measure, densities=densities, groups=groups, conf=conf, boot=my.boot)
  class(out) <- c('brainGraph_boot', class(out))
  return(out)
}

#' @inheritParams brainGraph_boot
#' @export
#' @rdname brainGraph_boot
boot_global <- function(densities, resids, R=1e3,
                        measure=c('mod', 'E.global', 'Cp', 'Lp', 'assortativity',
                                  'strength', 'mod.wt', 'E.global.wt'),
                        conf=0.95, .progress=TRUE) {
  .Deprecated('brainGraph_boot')
  brainGraph_boot(densities, resids, R=R, measure=measure, conf=conf, .progress=.progress)
}

#' Print a summary from a bootstrap analysis
#'
#' @param object A \code{brainGraph_boot} object (from
#'   \code{\link{brainGraph_boot}})
#' @param ... Unused
#' @importFrom boot boot.ci
#' @export
#' @method summary brainGraph_boot

summary.brainGraph_boot <- function(object, ...) {
  # Get everything into a data.table
  boot.dt <- with(object,
                  data.table(Group=rep(groups, each=length(densities)),
                             density=rep(densities, length(groups)),
                             Observed=c(sapply(boot, with, t0)),
                             se=c(sapply(boot, function(x) apply(x$t, 2, sd)))))
  ci <- with(object,
             vapply(seq_along(densities), function(x)
                    sapply(boot, function(y)
                           boot.ci(y, type='norm', index=x, conf=conf)$normal[2:3]),
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
  boot.sum <- list(meas.full=meas.full, dt=boot.dt, conf=object$conf, R=object$boot[[1]]$R)
  class(boot.sum) <- c('summary.brainGraph_boot', class(boot.sum))
  boot.sum
}

#' @aliases summary.brainGraph_boot
#' @method print summary.brainGraph_boot

print.summary.brainGraph_boot <- function(x, ...) {
  message('Bootstrap analysis\n', rep('-', getOption('width') / 2))
  cat('Graph metric: ', x$meas.full, '\n')
  cat('Number of bootstrap samples generated: ', x$R, '\n')
  conf.pct <- 100 * x$conf
  cat(conf.pct, '% confidence intervals\n\n')
  print(x$dt)
  invisible(x)
}

#' Plot bootstrap output of global graph measures across densities
#'
#' This function plots the output of \code{\link{brainGraph_boot}}, returning two
#' \code{ggplot} objects: one with shaded regions based on the standard error,
#' and the other based on confidence intervals (calculated using the normal
#' approximation.
#'
#' @param x The object output from \code{\link{brainGraph_boot}}
#' @param ... Not used
#' @param alpha A numeric indicating the opacity for
#'   \code{\link[ggplot2]{geom_ribbon}}
#' @export
#' @method plot brainGraph_boot
#'
#' @return A \emph{list} with the following elements:
#'   \item{se}{A ggplot object with ribbon representing standard error}
#'   \item{ci}{A ggplot object with ribbon representing confidence intervals}
#' @seealso \code{\link{brainGraph_boot}, \link[boot]{boot.ci}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot.brainGraph_boot <- function(x, ..., alpha=0.4) {
  Observed <- Group <- se <- ci.low <- ci.high <- NULL

  boot.sum <- summary(x)
  boot.dt <- boot.sum$dt
  b <- ggplot(boot.dt, aes(x=density, y=Observed, col=Group)) +
    geom_line() +
    labs(y=boot.sum$meas.full)

  se <- b + geom_ribbon(aes(ymin=Observed-se, ymax=Observed+se, fill=Group), alpha=alpha)
  ci <- b + geom_ribbon(aes(ymin=ci.low, ymax=ci.high, fill=Group), alpha=alpha)
  return(list(se=se, ci=ci))
}
