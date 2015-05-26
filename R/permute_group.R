#' Permutation test for group difference of graph measures
#'
#' This function draws permutations from linear model residuals to determine the
#' significance of between-group differences of a global graph measure. This
#' function is intended for cortical thickness networks (in which there is only
#' one graph per group), but is not restricted to that type of data.
#'
#' The \emph{graph} "level" will calculate modularity (Louvain algorithm),
#' clustering coefficient, average path length, degree assortativity, global
#' efficiency, and lobe assortativity.
#'
#' The \emph{lobe} "level" is intended to test for group differences in
#' inter-lobar connections, e.g. from the temporal lobe to the rest of the
#' brain.
#'
#' The \emph{asymmetry} "level" is intended to test for group differences in an
#' "asymmetry index", which is a measure of difference in intra-hemispheric
#' connection density.
#'
#' @param densities A vector of graph densities to loop through
#' @param resids A data.table of the residuals (from \code{\link{get.resid}})
#' @param num.subjs A vector of length 2 indicating group sizes
#' @param num.perms The number of permutations to perform (default: 1e3)
#' @param level A character string for the attribute level to calculate
#' differences; either 'graph', 'vertex', 'lobe', or 'asymmetry'
#' @param atlas Character string of the atlas name
#' @param atlas.dt Data table containing the specific atlas data
#' @export
#'
#' @return A data table with values for group differences in modularity, global
#' efficiency, clustering, average path length, and assortativity
#'
#' @seealso \code{\link[igraph]{centr_betw}, \link{count_interlobar},
#' \link{edge_asymmetry}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' out <- permute.group(densities, m$resids, summary(covars$Group), 1e3, 'graph',
#'   atlas='dk', atlas.dt)
#' out <- permute.group(densities, m$resids, summary(covars$Group), 1e3, 'vertex')
#' }

permute.group <- function(densities, resids, num.subjs, num.perms=1e3,
                          level=c('graph', 'vertex', 'lobe', 'asymmetry'),
                          atlas, atlas.dt=NULL) {
  n1 <- as.numeric(num.subjs[1])
  n.all <- sum(num.subjs)

  level <- match.arg(level)
  out <- foreach(i=seq_len(num.perms), .combine='rbind',
                 .packages='plyr', .export='assign_lobes') %dopar% {
    shuffled <- sample(n.all)
    corrs1 <- lapply(densities, function(x)
                     corr.matrix(resids[shuffled[1:n1]][, !'Group', with=F],
                                 density=x))
    corrs2 <- lapply(densities, function(x)
                     corr.matrix(resids[shuffled[(n1 +1):n.all]][, !'Group',
                                 with=F], density=x))
    g1 <- lapply(corrs1, function(x)
                 simplify(graph.adjacency(x$r.thresh, mode='undirected')))
    g2 <- lapply(corrs2, function(x)
                 simplify(graph.adjacency(x$r.thresh, mode='undirected')))

    if (level == 'vertex') {
      btwn.diff <- mapply(function(x, y) centr_betw(x)$res - centr_betw(y)$res,
                          g1, g2, SIMPLIFY=T)
      tmp <- as.data.table(cbind(densities, t(btwn.diff)))

    } else {
      g1 <- lapply(g1, set_graph_attr, 'atlas', atlas)
      g2 <- lapply(g2, set_graph_attr, 'atlas', atlas)
      g1 <- lapply(g1, assign_lobes, atlas.dt)
      g2 <- lapply(g2, assign_lobes, atlas.dt)

      if (level == 'graph') {
        mod1 <- sapply(g1, function(x) modularity(cluster_louvain(x)))
        mod2 <- sapply(g2, function(x) modularity(cluster_louvain(x)))
        mod.diff <- mod1 - mod2
        Cp1 <- sapply(g1, transitivity, type='localaverage')
        Cp2 <- sapply(g2, transitivity, type='localaverage')
        Cp.diff <- Cp1 - Cp2
        Lp1 <- sapply(g1, average.path.length)
        Lp2 <- sapply(g2, average.path.length)
        Lp.diff <- Lp1 - Lp2
        assortativity1 <- sapply(g1, assortativity.degree)
        assortativity2 <- sapply(g2, assortativity.degree)
        assort.diff <- assortativity1 - assortativity2
        E.global1 <- sapply(g1, graph.efficiency, 'global')
        E.global2 <- sapply(g2, graph.efficiency, 'global')
        E.global.diff <- E.global1 - E.global2
  
        assort.lobe1 <- sapply(g1, function(x) assortativity_nominal(x, V(x)$lobe))
        assort.lobe2 <- sapply(g2, function(x) assortativity_nominal(x, V(x)$lobe))
        assort.lobe.diff <- assort.lobe1 - assort.lobe2
        tmp <- data.table(density=densities, mod=mod.diff, E.global=E.global.diff,
                   Cp=Cp.diff, Lp=Lp.diff, assortativity=assort.diff,
                   assortativity.lobe=assort.lobe.diff)
  
      } else if (level == 'lobe') {
        t1 <- as.data.table(ldply(g1, count_interlobar, 'Temporal', atlas.dt))
        t2 <- as.data.table(ldply(g2, count_interlobar, 'Temporal', atlas.dt))
        tdiff <- t1 - t2
        tmp <- data.table(density=densities, diff=tdiff)
  
      } else if (level == 'asymmetry') {
        asymm1 <- ldply(g1, edge_asymmetry)$asymm
        asymm2 <- ldply(g2, edge_asymmetry)$asymm
        adiff <- asymm1 - asymm2
        tmp <- data.table(density=densities, asymm=adiff)
      }
    }

    tmp
  }
  return(out)
}
