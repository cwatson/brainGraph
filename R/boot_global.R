#' Bootstrapping for global graph measures
#'
#' This function performs bootstrapping to get group standard error estimates of
#' a global graph measure (e.g. modularity). It will output a list containing a
#' data table with standard errors and 95\% confidence intervals at each density
#' for each group, and 2 ggplot objects for plotting. This function is intended
#' for cortical thickness networks (in which there is only one graph per group).
#'
#' @param densities A vector of graph densities to loop through
#' @param resids A data.table of the residuals (from \code{\link{get.resid}})
#' @param groups A character vector indicating group names
#' @param num.reps The number of bootstrap replicates (default: 1e3)
#' @param measure Character string of the measure to test (default: 'mod')
#' @export
#'
#' @return A list with the following elements:
#' \item{dt}{A data table with length \emph{\# densities * \# groups}}
#' \item{p1}{A ggplot object with ribbon representing standard error}
#' \item{p2}{A ggplot object with ribbon representing 95\% confidence interval}
#' @seealso \code{\link[boot]{boot}, \link[boot]{boot.ci},
#' \link{permute.global}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' boot.res <- boot_global(densities, m$resids, groups, 1e3, 'E.global')
#' }

boot_global <- function(densities, resids, groups, num.reps=1e3, measure='mod') {

  # 'statistic' function for the bootstrapping process
  statfun <- function(x, i, measure='mod') {
    group <- x[i]
    corrs <- lapply(densities, function(x) corr.matrix(group, density=x))
    g.boot <- lapply(corrs, function(x)
                     simplify(graph.adjacency(x$r.thresh, mode='undirected')))
    if (measure == 'mod') {
      res <- sapply(g.boot, function(x) modularity(cluster_louvain(x)))
    } else if (measure == 'E.global') {
      res <- sapply(g.boot, graph.efficiency, 'global')
    } else if (measure == 'Cp') {
      res <- sapply(g.boot, transitivity, type='localaverage')
    } else if (measure == 'Lp') {
      res <- sapply(g.boot, average.path.length)
    } else if (measure == 'assortativity') {
      res <- sapply(g.boot, assortativity.degree)
    }
    return(res)
  }

  kNumDensities <- length(densities)
  boot1 <- boot(resids[groups[1]][, !'Group', with=F], statfun, R=num.reps,
                   parallel='multicore', ncpus=detectCores())
  boot2 <- boot(resids[groups[2]][, !'Group', with=F], statfun, R=num.reps,
                   parallel='multicore', ncpus=detectCores())
  boot.dt <- data.table(density=rep(densities, 2),
                        meas=c(boot1$t0, boot2$t0),
                        se=c(apply(boot1$t, 2, sd), apply(boot2$t, 2, sd)),
                        Group=rep(groups, each=kNumDensities))
  bootplot <- ggplot(boot.dt, aes(x=density, y=meas, col=Group)) +
    geom_line() +
    geom_ribbon(aes(ymin=meas-se, ymax=meas+se, fill=Group), alpha=0.3)

  # Use the estimated normal 95% CI instead of se
  ci1 <- sapply(seq_len(kNumDensities), function(x) boot.ci(boot1, type='norm',
                                                           index=x)$normal[2:3])
  ci2 <- sapply(seq_len(kNumDensities), function(x) boot.ci(boot2, type='norm',
                                                           index=x)$normal[2:3])
  boot.dt$ci.low <- c(ci1[1, ], ci2[1, ])
  boot.dt$ci.high <- c(ci1[2, ], ci2[2, ])
  bootplot.ci <- ggplot(boot.dt, aes(x=density, y=meas, col=Group)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin=ci.low, ymax=ci.high, fill=Group), alpha=0.3)

  return(list(dt=boot.dt, p1=bootplot, p2=bootplot.ci))
}
