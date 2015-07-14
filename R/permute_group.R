#' Permutation test for group difference of graph measures
#'
#' This function draws permutations from linear model residuals to determine the
#' significance of between-group differences of a global graph measure. This
#' function is intended for cortical thickness networks (in which there is only
#' one graph per group), but is not restricted to that type of data.
#'
#' The \emph{graph} "level" will calculate modularity (Louvain algorithm),
#' clustering coefficient, average path length, degree assortativity, global
#' efficiency, lobe assortativity, and asymmetry.
#'
#' The \emph{lobe} "level" is intended to test for group differences in
#' inter-lobar connections, e.g. from the temporal lobe to the rest of the
#' brain.
#'
#' @param density Numeric; the density of the resultant graphs
#' @param resids A data table of the residuals (from \code{\link{get.resid}})
#' @param num.subjs A vector of length 2 indicating group sizes
#' @param num.perms The number of permutations to perform (default: 1e3)
#' @param level A character string for the attribute level to calculate
#' differences; either 'graph', 'vertex', or 'lobe'
#' @param atlas Character string of the atlas name
#' @param atlas.dt A data table containing the specific atlas data
#' @export
#'
#' @return A data table with values for group differences in modularity, global
#' efficiency, clustering, average path length, and assortativity (etc.)
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

permute.group <- function(density, resids, num.subjs, num.perms=1e3,
                          level=c('graph', 'vertex', 'lobe'),
                          atlas, atlas.dt=NULL) {
  n1 <- as.numeric(num.subjs[1])
  n.all <- sum(num.subjs)

  level <- match.arg(level)
  out <- foreach(i=seq_len(num.perms), .combine='rbind',
                 .export='assign_lobes') %dopar% {
    shuffled <- sample(n.all)
    corrs1 <- corr.matrix(resids[shuffled[1:n1]][, !'Group', with=F],
                          density=density)
    corrs2 <- corr.matrix(resids[shuffled[(n1+1):n.all]][, !'Group', with=F],
                          density=density)
    g1 <- graph_from_adjacency_matrix(corrs1$r.thresh, mode='undirected', diag=F)
    g2 <- graph_from_adjacency_matrix(corrs2$r.thresh, mode='undirected', diag=F)

    if (level == 'vertex') {
      btwn.diff <- centr_betw(g1)$res - centr_betw(g2)$res
      tmp <- as.data.table(cbind(density, t(btwn.diff)))

    } else {
      g1$atlas <- g2$atlas <- atlas
      g1 <- assign_lobes(g1, atlas.dt, rand=T)
      g2 <- assign_lobes(g2, atlas.dt, rand=T)

      if (level == 'graph') {
        mod1 <- modularity(cluster_louvain(g1))
        mod2 <- modularity(cluster_louvain(g2))
        mod.diff <- mod1 - mod2
        Cp1 <- transitivity(g1, type='localaverage')
        Cp2 <- transitivity(g2, type='localaverage')
        Cp.diff <- Cp1 - Cp2
        Lp1 <- average.path.length(g1)
        Lp2 <- average.path.length(g2)
        Lp.diff <- Lp1 - Lp2
        assortativity1 <- assortativity.degree(g1)
        assortativity2 <- assortativity.degree(g2)
        assort.diff <- assortativity1 - assortativity2
        E.global1 <- graph.efficiency(g1, 'global')
        E.global2 <- graph.efficiency(g2, 'global')
        E.global.diff <- E.global1 - E.global2

        assort.lobe1 <- assortativity_nominal(g1, V(g1)$lobe)
        assort.lobe2 <- assortativity_nominal(g2, V(g2)$lobe)
        assort.lobe.diff <- assort.lobe1 - assort.lobe2
        asymm1 <- edge_asymmetry(g1)$asymm
        asymm2 <- edge_asymmetry(g2)$asymm
        adiff <- asymm1 - asymm2
        tmp <- data.table(density=density, mod=mod.diff, E.global=E.global.diff,
                   Cp=Cp.diff, Lp=Lp.diff, assortativity=assort.diff,
                   assortativity.lobe=assort.lobe.diff, asymm=adiff)

      } else if (level == 'lobe') {
        t1 <- as.data.table(count_interlobar(g1, 'Temporal', atlas.dt))
        t2 <- as.data.table(count_interlobar(g2, 'Temporal', atlas.dt))
        tdiff <- t1 - t2
        tmp <- data.table(density=density, diff=tdiff)

      }
    }

    tmp
  }
  return(out)
}
