#' Permutation test for group difference of graph measures
#'
#' This function draws permutations from linear model residuals to determine the
#' significance of between-group differences of a global or vertex-wise graph
#' measure. This function is intended for cortical thickness networks (in which
#' there is only one graph per group), but is not restricted to that type of
#' data.
#'
#' The \emph{graph} "level" will calculate modularity (Louvain algorithm),
#' clustering coefficient, average path length, degree assortativity, global
#' efficiency, lobe assortativity, and edge asymmetry.
#'
#' The \emph{vertex} "level" will calculate either betweenness centrality or
#' vulnerability (of which only the maximum is taken, thus making this
#' technically a \emph{graph} level measure).
#'
#' The \emph{lobe} "level" is intended to test for group differences in
#' inter-lobar connections, e.g. from the temporal lobe to the rest of the
#' brain.
#'
#' @param permSet The set of permutations to loop through (obtained from
#' \code{\link[permute]{shuffleSet}})
#' @param density Numeric; the density of the resultant graphs
#' @param resids A data table of the residuals (from \code{\link{get.resid}})
#' @param level A character string for the attribute level to calculate
#' differences; either 'graph', 'vertex', or 'lobe'
#' @param atlas Character string of the atlas name
#' @param atlas.dt A data table containing the specific atlas data
#' @param measure A character string, either 'btwn.cent' or 'vulnerability'
#' (specific to the vertex \emph{level})
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
#' m <- get.resid(all.thick, covars)
#' myPerms <- shuffleSet(n=nrow(m$resids), nset=1e3)
#' out <- permute.group(myPerms, densities, m$resids, 1e3, 'graph', atlas='dk',
#'    atlas.dt)
#' out <- permute.group(myPerms, densities, m$resids, 1e3, 'vertex')
#' }

permute.group <- function(permSet, density, resids,
                          level=c('graph', 'vertex', 'lobe'),
                          atlas, atlas.dt=NULL, measure=c('btwn.cent', 'vulnerability')) {
  i <- NULL
  level <- match.arg(level)
  measure <- match.arg(measure)

  groups <- as.numeric(resids$Group)
  out <- foreach(i=seq_len(nrow(permSet)), .combine='rbind',
                 .export='assign_lobes') %dopar% {
    corrs1 <- corr.matrix(resids[which(groups[permSet[i, ]] == 1), !'Group', with=F],
                          density=density)
    corrs2 <- corr.matrix(resids[which(groups[permSet[i, ]] == 2), !'Group', with=F],
                          density=density)
    g1 <- graph_from_adjacency_matrix(corrs1$r.thresh, mode='undirected', diag=F)
    g2 <- graph_from_adjacency_matrix(corrs2$r.thresh, mode='undirected', diag=F)

    if (level == 'vertex') {
      if (measure == 'vulnerability') {
        g1$E.global <- graph.efficiency(g1, 'global')
        g2$E.global <- graph.efficiency(g2, 'global')
        V(g1)$degree <- degree(g1)
        V(g2)$degree <- degree(g2)
        vuln.diff <- max(vulnerability(g1), .parallel=F) - max(vulnerability(g2), .parallel=F)
        tmp <- data.table(density=density, vulnerability=vuln.diff)
      } else {
        btwn.diff <- centr_betw(g1)$res - centr_betw(g2)$res
        tmp <- as.data.table(cbind(density, t(btwn.diff)))
      }

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
