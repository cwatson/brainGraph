#' Bootstrapping for global graph measures
#'
#' This function performs bootstrapping to get group standard error estimates of
#' a global graph measure (e.g. modularity). It will output a list containing a
#' data table with standard errors and 95\% confidence intervals at each density
#' for each group, and 2 ggplot objects for plotting. This function is intended
#' for cortical thickness networks (in which there is only one graph per group),
#' but will obviously work in other scenarios.
#'
#' The 95\% confidence intervals are calculated using the normal approximation.
#'
#' @param densities A vector of graph densities to loop through
#' @param resids A data.table of the residuals (from \code{\link{get.resid}})
#' @param groups A character vector indicating group names
#' @param R The number of bootstrap replicates (default: 1e3)
#' @param measure Character string of the measure to test (default: 'mod')
#' @export
#'
#' @return A list with the following elements:
#' \item{g}{A list of \code{\link[boot]{boot}} objects}
#' \item{dt}{A data table with length \emph{# densities * # groups}}
#' \item{p1}{A ggplot object with ribbon representing standard error}
#' \item{p2}{A ggplot object with ribbon representing 95\% confidence interval}
#' @seealso \code{\link[boot]{boot}, \link[boot]{boot.ci},
#' \link{permute.group}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' boot.res <- boot_global(densities, m$resids, groups, 1e3, 'E.global')
#' }

boot_global <- function(densities, resids, groups, R=1e3, measure='mod') {

  meas <- Group <- se <- ci.low <- ci.high <- NULL
  kNumGroups <- length(groups)
  kNumDensities <- length(densities)
  # 'statistic' function for the bootstrapping process
  statfun <- function(x, i, measure) {
    group <- x[i]
    corrs <- lapply(densities, function(x) corr.matrix(group, density=x))
    g.boot <- lapply(corrs, function(x)
                     graph_from_adjacency_matrix(x$r.thresh, mode='undirected',
                                                 diag=F))
    if (measure == 'mod') {
      res <- vapply(g.boot, function(x) modularity(cluster_louvain(x)), numeric(1))
    } else if (measure == 'E.global') {
      res <- vapply(g.boot, graph.efficiency, numeric(1), 'global')
    } else if (measure == 'Cp') {
      res <- vapply(g.boot, transitivity, numeric(1), type='localaverage')
    } else if (measure == 'Lp') {
      res <- vapply(g.boot, average.path.length, numeric(1))
    } else if (measure == 'assortativity') {
      res <- vapply(g.boot, assortativity.degree, numeric(1))
    }
    return(res)
  }

  # Show a progress bar so you aren't left in the dark
  intfun <- function(data, indices, measure) {
    curVal <- get('counter', envir=env) + ncpus
    assign('counter', curVal, envir=env)
    setTxtProgressBar(get('progbar', envir=env), curVal)
    flush.console()
    statfun(data, indices, measure)
  }

  ncpus <- detectCores()
  env <- environment()
  my.boot <- vector('list', length=kNumGroups)
  for (i in seq_len(kNumGroups)) {
    counter <- 0
    progbar <- txtProgressBar(min=0, max=R, style=3)

    my.boot[[i]] <- boot(resids[groups[i], !'Group', with=F], intfun,
                       measure=measure, R=R, parallel='multicore', ncpus=ncpus)
    close(progbar)
  }

  # Get everything into a data.table
  boot.dt <- data.table(density=rep(densities, kNumGroups),
                        meas=c(sapply(my.boot, with, t0)),
                        se=c(sapply(my.boot, function(x) apply(x$t, 2, sd))),
                        Group=rep(groups, each=kNumDensities))
  bootplot <- ggplot(boot.dt, aes(x=density, y=meas, col=Group)) +
    geom_line() +
    geom_ribbon(aes(ymin=meas-se, ymax=meas+se, fill=Group), alpha=0.3) +
    ylab(measure)

  # Use the estimated normal 95% CI instead of se
  ci <- vapply(seq_along(densities), function(x)
               sapply(my.boot, function(y)
                      boot.ci(y, type='norm', index=x)$normal[2:3]),
               numeric(kNumGroups * 2))
  boot.dt$ci.low <- c(t(ci[2 * (seq_along(groups) - 1) + 1, ]))
  boot.dt$ci.high <- c(t(ci[2 * (seq_along(groups) - 1) + 2, ]))
  bootplot.ci <- ggplot(boot.dt, aes(x=density, y=meas, col=Group)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=ci.low, ymax=ci.high, fill=Group), alpha=0.3) +
    ylab(measure)

  return(list(g=my.boot, dt=boot.dt, p1=bootplot, p2=bootplot.ci))
}
