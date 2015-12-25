#' Calculate permutation p-values and plot group differences
#'
#' For a given (global- or vertex-level) graph measure, determine the
#' permutation p-value and create a plot showing group differences, either
#' across densities or regions. You may specify the \eqn{\alpha}-level; a red
#' asterisk is added if \eqn{p < \alpha} and a blue asterisk is added if
#' \eqn{\alpha < p < 0.1} (i.e. a "trend"). You may also choose whether you want
#' a one- or two-sided test.
#'
#' @param g1 List of igraph graph objects for group 1
#' @param g2 List of igraph graph objects for group 2
#' @param perm.dt Data table with the permutation results
#' @param measure Character string for the graph measure of interest
#' @param level Character string, either 'graph' or 'vertex'
#' @param alternative Character string, whether to do a two- or one-sided test
#' (default: 'two.sided')
#' @param alpha Significance level (default: 0.05)
#' @param groups Character vector of group names (default: NULL)
#' @param ylabel Character string for y-axis label (default: NULL)
#' @export
#'
#' @return A list with three elements:
#' \item{dt}{A data table with p-values for each density/region}
#' \item{p1}{A \code{\link[ggplot2]{ggplot}} plotting object}
#' \item{p2}{A \code{\link[ggplot2]{ggplot}} plotting object}
#'
#' @seealso \code{\link{permute.group}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' perms.mod.sig <- perms.sig(g[[1]], g[[2]], perms.all, 'mod', level='graph',
#'   'less', groups, ylabel='Modularity')
#' perms.mod.btwn <- perms.sig(g[[1]], g[[2]], perms.btwn, 'btwn.cent',
#'   level='vertex')
#' }

plot_perm_diffs <- function(g1, g2, perm.dt, measure,
                            level=c('graph', 'vertex'),
                            alternative=c('two.sided', 'less', 'greater'),
                            alpha=0.05, groups=NULL, ylabel=NULL) {

  p <- perm.diff <- obs.diff <- sig <- trend <- yloc <- obs <- Group <- mean.diff <- NULL
  ci.low <- ci.high <- region <- reg.num <- NULL
  densities.perm <- perm.dt[, unique(density)]
  densities.g <- which(round(sapply(g1, graph_attr, 'density'), 2) %in% round(densities.perm, 2))
  g1 <- g1[densities.g]
  g2 <- g2[densities.g]
  N <- perm.dt[, .N, by=density]$N  # Handles diff num. of perm's across densities

  if (is.null(groups)) groups <- c('Group 1', 'Group 2')
  alt <- match.arg(alternative)
  level <- match.arg(level)
  # Graph-level permutations
  #-------------------------------------
  if (level == 'graph') {
    if (!measure %in% names(perm.dt)) {
      stop(sprintf('Measure %s is not valid!', deparse(substitute(measure))))
    }
    if (is.igraph(g1)) {
      meas.obs1 <- graph_attr(g1, measure)
      meas.obs2 <- graph_attr(g2, measure)
    } else {
      meas.obs1 <- sapply(g1, graph_attr, measure)
      meas.obs2 <- sapply(g2, graph_attr, measure)
    }
    perm.dt <- perm.dt[, list(density, perm.diff=get(measure))]
    perm.dt$obs.diff <- rep(meas.obs1 - meas.obs2, times=N)

    if (alt == 'two.sided') {
      perm.dt[, p := (sum(abs(perm.diff) >= abs(unique(obs.diff))) + 1) / (.N + 1),
              by=density]
      ci <- c(alpha / 2, 1 - (alpha / 2))

    } else if (alt == 'less') {
      perm.dt[, p := (sum(perm.diff <= unique(obs.diff)) + 1) / (.N + 1), by=density]
      ci <- c(alpha, NULL)

    } else if (alt == 'greater') {
      perm.dt[, p := (sum(perm.diff >= unique(obs.diff)) + 1) / (.N + 1), by=density]
      ci <- c(NULL, 1 - alpha)
    }

    result.dt <- data.table(Group=rep(groups, each=length(densities.perm)),
                            density=densities.perm,
                            N=rep(N, length(groups)),
                            obs=c(meas.obs1, meas.obs2),
                            perm.diff=rep(perm.dt[, mean(perm.diff), by=density]$V1, length(groups)),
                            p=perm.dt[, unique(p), by=density]$V1,
                            sig='', trend='')
    result.dt[, p.fdr := p.adjust(p, 'fdr')]
    result.dt[p < alpha, sig := '*']
    result.dt[p >= alpha & p < 0.1, trend := '*']
    result.dt[, yloc := round(min(obs) - 0.05 * diff(range(obs)), 3)]
    sigplot <- ggplot(data=result.dt, aes(x=density, y=obs, col=Group)) +
      geom_line() +
      geom_text(aes(y=yloc, label=sig), col='red', size=8) +
      geom_text(aes(y=yloc, label=trend), col='blue', size=8) +
      theme(legend.position='none')
    if (!is.null(ylabel)) {
      sigplot <- sigplot + ylab(ylabel)
    } else {
      sigplot <- sigplot + ylab(measure)
    }

    perm.dt[, mean.diff := mean(perm.diff), by=density]
    perm.dt[, c('ci.low', 'ci.high') := as.list(sort(perm.diff)[.N * ci]), by=density]
    sigplot2 <- ggplot(result.dt[, list(obs.diff=-diff(obs)), by=density],
                       aes(x=as.factor(density))) +
      geom_line(data=perm.dt[, list(ci.low=unique(ci.low)),  by=density],
                aes(x=density, y=ci.low), lty=2) +
      geom_line(data=perm.dt[, list(ci.high=unique(ci.high)), by=density],
                aes(x=density, y=ci.high), lty=2) +
      geom_line(aes(x=density, y=obs.diff), col='red') +
      geom_point(aes(x=density, y=obs.diff), col='red', size=3) +
      geom_line(data=perm.dt, aes(x=density, y=mean.diff), lty=2) +
      #geom_boxplot(data=perm.dt, aes(x=as.factor(density), y=perm.diff),
      #             fill='cyan3', outlier.size=0) +
      #scale_x_discrete(breaks=seq(densities.perm[1],
      #                            densities.perm[length(densities.perm)], by=0.05)) +
      xlab('Density') +
      ylab(sprintf('Permutation difference (%s - %s)', groups[1], groups[2])) +
      ggtitle(ylabel)

  # Vertex-level permutations
  #-------------------------------------
  } else if (level == 'vertex') {
    if (is.igraph(g1)) {
      meas.obs1 <- vertex_attr(g1, measure)
      meas.obs2 <- vertex_attr(g2, measure)
    } else {
      meas.obs1 <- sapply(g1, vertex_attr, measure)
      meas.obs2 <- sapply(g2, vertex_attr, measure)
    }
    perm.dt[, N := .N, by=density]
    perm.m <- melt(perm.dt, id.vars=c('density', 'N'), variable.name='region',
                   value.name='perm.diff')
    setkey(perm.m, density, region)
    perm.m$obs.diff <- rep(as.vector(meas.obs1 - meas.obs2),
                           times=rep(N, each=ncol(perm.dt)-2))

    if (alt == 'two.sided') {
      perm.m[, p := (sum(abs(perm.diff) >= abs(unique(obs.diff))) + 1) / (.N + 1),
             by=list(density, region)]
    } else if (alt == 'less') {
      perm.m[, p := (sum(perm.diff <= unique(obs.diff)) + 1) / (.N + 1),
             by=list(density, region)]
    } else if (alt == 'greater') {
      perm.m[, p := (sum(perm.diff >= unique(obs.diff)) + 1) / (.N + 1),
             by=list(density, region)]
    }
    p.fdr <- perm.m[, list(p=unique(p)), by=list(region, density)][, p.adjust(p, 'fdr'), by=density]$V1
    perm.m$p.fdr <- rep(p.fdr, times=rep(N, each=ncol(perm.dt)-2))

    perm.m.sig <- perm.m[p < alpha]
    perm.m.sig[, reg.num := rep(seq_len(length(unique(region))), each=unique(N)),
               by=density]

    sigplot <- ggplot(perm.m.sig, aes(x=region, y=perm.diff)) +
      geom_boxplot(fill='cyan3', outlier.size=0) +
      geom_segment(data=perm.m.sig[, list(obs.diff=unique(obs.diff)), by=list(density, reg.num)],
                   aes(x=reg.num-0.25, xend=reg.num+0.25, y=obs.diff, yend=obs.diff), col='red') +
      facet_wrap(~ density, scales='free')

      sigplot <- sigplot +
        ylab(sprintf('Permutation difference (%s - %s)', groups[1], groups[2]))
      result.dt <- perm.m
      sigplot2 <- NULL
  }

  return(list(dt=result.dt, p1=sigplot, p2=sigplot2))
}
