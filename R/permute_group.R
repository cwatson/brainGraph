#' Permutation test for group difference of graph measures
#'
#' This function draws permutations from linear model residuals to determine the
#' significance of between-group differences of a global graph measure. This
#' function is intended for cortical thickness networks (in which there is only
#' one graph per group), but is not restricted to that type of data.
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
#' @param atlas Character string for the specific atlas to use
#' @param atlas.list List containing the specific atlas data
#' differences; either 'graph', 'vertex', 'lobe', or 'asymmetry'
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
#' out <- permute.group(densities, m$resids, summary(covars$Group), 1e3, 'graph')
#' out <- permute.group(densities, m$resids, summary(covars$Group), 1e3, 'lobe',
#'   atlas, atlas.list)
#' }

permute.group <- function(densities, resids, num.subjs, num.perms=1e3,
                          level=c('graph', 'vertex', 'lobe', 'asymmetry'),
                          atlas=NULL, atlas.list=NULL) {
  n1 <- num.subjs[1]
  n.all <- sum(num.subjs)

  level <- match.arg(level)
  out <- foreach(i=seq_len(num.perms), .combine='rbind') %dopar% {
    shuffled <- sample(n.all)
    corrs1 <- lapply(densities, function(x)
                     corr.matrix(resids[shuffled[1:n1]][, !'Group', with=F],
                                 density=x))
    corrs2 <- lapply(densities, function(x)
                     corr.matrix(resids[shuffled[(n1 +1):n.all]][, !'Group', with=F],
                                 density=x))
    g1 <- lapply(corrs1, function(x)
                 simplify(graph.adjacency(x$r.thresh, mode='undirected')))
    g2 <- lapply(corrs2, function(x)
                 simplify(graph.adjacency(x$r.thresh, mode='undirected')))

    if (level == 'vertex') {
      btwn.diff <- mapply(function(x, y) centr_betw(x)$res - centr_betw(y)$res,
                          g1, g2, SIMPLIFY=T)
      tmp <- as.data.table(cbind(densities, t(btwn.diff)))

    } else if (level == 'graph') {
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

      g1 <- lapply(g1, assign_lobes, atlas, atlas.list)
      g2 <- lapply(g2, assign_lobes, atlas, atlas.list)
      assort.lobe1 <- sapply(g1, function(x) assortativity_nominal(x, V(x)$lobe))
      assort.lobe2 <- sapply(g2, function(x) assortativity_nominal(x, V(x)$lobe))
      assort.lobe.diff <- assort.lobe1 - assort.lobe2
      tmp <- data.table(density=densities, mod=mod.diff, E.global=E.global.diff,
                 Cp=Cp.diff, Lp=Lp.diff, assortativity=assort.diff,
                 assortativity.lobe=assort.lobe.diff)
      #tmp <- data.table(density=densities, assortativity.lobe=assort.lobe.diff)

    } else if (level == 'lobe') {
      g1 <- lapply(g1, assign_lobes, atlas, atlas.list)
      g2 <- lapply(g2, assign_lobes, atlas, atlas.list)

      t1 <- as.data.table(ldply(g1, count_interlobar, 'Temporal', atlas.list))
      t2 <- as.data.table(ldply(g2, count_interlobar, 'Temporal', atlas.list))
      tdiff <- t1 - t2
      tmp <- data.table(density=densities, diff=tdiff)

    } else if (level == 'asymmetry') {
      g1 <- lapply(g1, assign_lobes, atlas, atlas.list)
      g2 <- lapply(g2, assign_lobes, atlas, atlas.list)

      asymm1 <- ldply(g1, edge_asymmetry)$asymm
      asymm2 <- ldply(g2, edge_asymmetry)$asymm
      adiff <- asymm1 - asymm2
      tmp <- data.table(density=densities, diff=adiff)
    }

    tmp
  }
  return(out)
}
