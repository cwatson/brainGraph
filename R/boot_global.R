#' Bootstrapping for global graph measures
#'
#' This function performs bootstrapping to get group standard error estimates of
#' a global graph measure (e.g. modularity). It will output a list containing
#' the \code{\link[boot]{boot}} objects and a \code{\link{data.table}} with
#' standard errors and 95\% confidence intervals at each density for each group.
#'
#' The 95\% confidence intervals are calculated using the normal approximation.
#'
#' @param densities Numeric vector of graph densities to loop through
#' @param resids A data.table of the residuals (from \code{\link{get.resid}})
#' @param R Integer; the number of bootstrap replicates (default: 1e3)
#' @param measure Character string of the measure to test (default: 'mod')
#' @param .progress Logical indicating whether or not to show a progress bar
#'   (default: \code{TRUE})
#' @param ... Arguments passed to \code{\link[boot]{boot.ci}} (specifically,
#'   \code{conf})
#' @export
#' @importFrom boot boot boot.ci
#'
#' @return A list with two elements:
#'   \item{boot}{A list of \code{\link[boot]{boot}} objects (one for each group)}
#'   \item{dt}{A data table with length \emph{# densities * # groups}}
#'
#' @family Group analysis functions
#' @seealso \code{\link[boot]{boot}, \link[boot]{boot.ci}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' boot.E.global <- boot_global(densities, resids.all, 1e3, 'E.global')
#' }

boot_global <- function(densities, resids, R=1e3,
                        measure=c('mod', 'E.global', 'Cp', 'Lp', 'assortativity',
                                  'strength', 'mod.wt', 'E.global.wt'),
                        .progress=TRUE, ...) {
  Group <- t0 <- NULL
  if (!is.factor(resids$Group)) resids[, Group := as.factor(Group)]
  groups <- resids[, levels(Group)]
  measure <- match.arg(measure)

  # 'statistic' function for the bootstrapping process
  statfun <- function(x, i, measure) {
    group <- as.matrix(x[i, !c('Study.ID', 'Group'), with=F])
    corrs <- lapply(densities, function(x) corr.matrix(group, density=x))
    if (measure %in% c('strength', 'mod.wt', 'E.global.wt')) {
      g.boot <- lapply(corrs, function(x)
                       graph_from_adjacency_matrix(x$R * x$r.thresh, mode='undirected',
                                                   diag=FALSE, weighted=TRUE))
    } else {
      g.boot <- lapply(corrs, function(x)
                       graph_from_adjacency_matrix(x$r.thresh, mode='undirected',
                                                   diag=FALSE))
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
  my.boot <- vector('list', length=length(groups))
  for (i in seq_along(groups)) {
    counter <- 0
    if (isTRUE(.progress)) progbar <- txtProgressBar(min=0, max=R, style=3)
    my.boot[[i]] <- boot(resids[groups[i]], intfun, measure=measure, R=R,
                         parallel=my.parallel, ncpus=ncpus, cl=cl)
    if (isTRUE(.progress)) close(progbar)
  }

  # Get everything into a data.table
  boot.dt <- data.table(density=rep(densities, length(groups)),
                        meas=c(sapply(my.boot, with, t0)),
                        se=c(sapply(my.boot, function(x) apply(x$t, 2, sd))),
                        Group=rep(groups, each=length(densities)))
  ci <- vapply(seq_along(densities), function(x)
               sapply(my.boot, function(y)
                      boot.ci(y, type='norm', index=x, ...)$normal[2:3]),
               numeric(2 * length(groups)))
  boot.dt$ci.low <- c(t(ci[2 * (seq_along(groups) - 1) + 1, ]))
  boot.dt$ci.high <- c(t(ci[2 * (seq_along(groups) - 1) + 2, ]))

  return(list(boot=my.boot, dt=boot.dt))
}
