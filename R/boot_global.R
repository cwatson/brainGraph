#' Bootstrapping for global graph measures
#'
#' This function performs bootstrapping to get group standard error estimates of
#' a global graph measure (e.g. modularity). It will output a list containing
#' the \code{\link[boot]{boot}} objects and a \code{\link{data.table}} with
#' standard errors and 95\% confidence intervals at each density for each group.
#'
#' The 95\% confidence intervals are calculated using the normal approximation.
#'
#' @param densities A vector of graph densities to loop through
#' @param resids A data.table of the residuals (from \code{\link{get.resid}})
#' @param R The number of bootstrap replicates (default: 1e3)
#' @param measure Character string of the measure to test (default: 'mod')
#' @export
#'
#' @return A list with two elements:
#' \item{g}{A list of \code{\link[boot]{boot}} objects (one for each group)}
#' \item{dt}{A data table with length \emph{# densities * # groups}}
#' @seealso \code{\link[boot]{boot}, \link[boot]{boot.ci},
#' \link{permute.group}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' boot.E.global <- boot_global(densities, m$resids, 1e3, 'E.global')
#' }

boot_global <- function(densities, resids, R=1e3,
                        measure=c('mod', 'E.global', 'Cp', 'Lp',
                                  'assortativity')) {

  Group <- t0 <- NULL
  groups <- resids[, levels(Group)]
  measure <- match.arg(measure)

  # 'statistic' function for the bootstrapping process
  statfun <- function(x, i, measure) {
    group <- as.matrix(x[i, !c('Study.ID', 'Group'), with=F])
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
    progbar <- txtProgressBar(min=0, max=R, style=3)

    my.boot[[i]] <- boot(resids[groups[i]], intfun,
                       measure=measure, R=R, parallel=my.parallel, ncpus=ncpus,
                       cl=cl)
    close(progbar)
  }

  # Get everything into a data.table
  boot.dt <- data.table(density=rep(densities, length(groups)),
                        meas=c(sapply(my.boot, with, t0)),
                        se=c(sapply(my.boot, function(x) apply(x$t, 2, sd))),
                        Group=rep(groups, each=length(densities)))
  ci <- vapply(seq_along(densities), function(x)
               sapply(my.boot, function(y)
                      boot.ci(y, type='norm', index=x)$normal[2:3]),
               numeric(2 * length(groups)))
  boot.dt$ci.low <- c(t(ci[2 * (seq_along(groups) - 1) + 1, ]))
  boot.dt$ci.high <- c(t(ci[2 * (seq_along(groups) - 1) + 2, ]))

  return(list(g=my.boot, dt=boot.dt))
}
