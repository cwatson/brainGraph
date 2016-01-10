#' Permutation test for group difference of graph measures
#'
#' This function draws permutations from linear model residuals to determine the
#' significance of between-group differences of a global or vertex-wise graph
#' measure. This function is intended for cortical thickness networks (in which
#' there is only one graph per group), but can be extended to other types of
#' data.
#'
#' The \emph{graph} "level" will calculate modularity (Louvain algorithm),
#' clustering coefficient, average path length, degree assortativity, global
#' efficiency, lobe assortativity, and edge asymmetry.
#'
#' The \emph{vertex} "level" will calculate a vertex-wise measure. Currently,
#' you can choose betweenness centrality, degree, nodal efficiency, k-nearest
#' neighbor degree, transitivity, or vulnerability.
#'
#' The \emph{lobe} "level" is intended to test for group differences in number
#' of inter-lobar connections, e.g. from the temporal lobe to the rest of the
#' brain.
#'
#' The \emph{other} "level" allows you to pass your own function to do
#' permutations with. This is useful if you want to calculate something that I
#' haven't hard-coded (e.g. number of hubs between groups). It must take as its
#' own arguments: "g1", "g2", and "density".
#'
#' @param permSet A matrix of the set of permutations to loop through; the
#' number of rows equals the desired number of permutations and the number of
#' columns equals the total number of subjects across groups
#' @param density Numeric; the density of the resultant graphs
#' @param resids A data table of the residuals (from \code{\link{get.resid}})
#' @param level A character string for the attribute level to calculate
#' differences; either 'graph', 'vertex', 'lobe', or 'other'
#' @param atlas Character string of the atlas name
#' @param measure A character string, either 'btwn.cent', 'degree', 'E.nodal',
#' 'knn', or 'transitivity', 'vulnerability' (specific to the vertex \emph{level})
#' @param .function A custom function you can pass (if \emph{level} is 'other')
#' @export
#'
#' @return A data table with values for group differences in modularity, global
#' efficiency, clustering, average path length, and assortativity (etc.)
#'
#' @seealso \code{\link[igraph]{centr_betw}, \link{vulnerability},
#' \link{count_interlobar}, \link{edge_asymmetry}, \link{graph.efficiency}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' m <- get.resid(all.thick, covars)
#' myPerms <- shuffleSet(n=nrow(m$resids), nset=1e3)
#' out <- permute.group(myPerms, densities[N], m$resids, 'graph', atlas='dk')
#' out <- permute.group(myPerms, densities[N], m$resids, 'vertex')
#' out <- permute.group(myPerms, densities[N], m$resids, 'other',
#'   .function=myFun)
#' }

permute.group <- function(permSet, density, resids,
                          level=c('graph', 'vertex', 'lobe', 'other'), atlas,
                          measure=c('btwn.cent', 'degree', 'E.nodal',
                                    'knn', 'transitivity', 'vulnerability'),
                          .function=NULL) {
  i <- NULL
  level <- match.arg(level)
  if (level == 'other') {
    if (is.null(.function) | !is.function(.function)) {
      stop(paste('Argument ".function" must be a function!'))
    }
  }
  measure <- match.arg(measure)

  groups <- as.numeric(resids$Group)
  out <- foreach(i=seq_len(nrow(permSet)), .combine='rbind',
                 .export='assign_lobes') %dopar% {
    corrs1 <- corr.matrix(as.matrix(resids[which(groups[permSet[i, ]] == 1),
                          !c('Study.ID', 'Group'), with=F]),
                          density=density)
    corrs2 <- corr.matrix(as.matrix(resids[which(groups[permSet[i, ]] == 2),
                          !c('Study.ID', 'Group'), with=F]),
                          density=density)
    g1 <- graph_from_adjacency_matrix(corrs1$r.thresh, mode='undirected', diag=F)
    g2 <- graph_from_adjacency_matrix(corrs2$r.thresh, mode='undirected', diag=F)

    # Vertex-level
    #-----------------------------------
    if (level == 'vertex') {
      if (measure == 'vulnerability') {
        g1$E.global <- graph.efficiency(g1, 'global')
        g2$E.global <- graph.efficiency(g2, 'global')
        V(g1)$degree <- degree(g1)
        V(g2)$degree <- degree(g2)
        vuln.diff <- vulnerability(g1, .parallel=F) -
          vulnerability(g2, .parallel=F)
        tmp <- as.data.table(cbind(density, t(vuln.diff)))

      } else if (measure == 'degree') {
        deg.diff <- degree(g1) - degree(g2)
        tmp <- as.data.table(cbind(density, t(deg.diff)))

      } else if (measure == 'E.nodal') {
        E.nodal.diff <- graph.efficiency(g1, 'nodal') - graph.efficiency(g2, 'nodal')
        tmp <- as.data.table(cbind(density, t(E.nodal.diff)))

      } else if (measure == 'knn') {
        knn.diff <- graph.knn(g1)$knn - graph.knn(g2)$knn
        tmp <- as.data.table(cbind(density, t(knn.diff)))

      } else if (measure == 'transitivity') {
        transitivity.diff <- transitivity(g1, type='local', isolates='zero') -
          transitivity(g2, type='local', isolates='zero')
        tmp <- as.data.table(cbind(density, t(transitivity.diff)))

      } else {
        btwn.diff <- centr_betw(g1)$res - centr_betw(g2)$res
        tmp <- as.data.table(cbind(density, t(btwn.diff)))
      }

    # Custom function
    #-----------------------------------
    } else if (level == 'other') {
      tmp <- .function(g1, g2, density)

    } else {
      g1$atlas <- g2$atlas <- atlas
      g1 <- assign_lobes(g1, rand=T)
      g2 <- assign_lobes(g2, rand=T)

      # Graph-level
      #-----------------------------------
      if (level == 'graph') {
        mod1 <- modularity(cluster_louvain(g1))
        mod2 <- modularity(cluster_louvain(g2))
        mod.diff <- mod1 - mod2
        Cp1 <- transitivity(g1, type='localaverage')
        Cp2 <- transitivity(g2, type='localaverage')
        Cp.diff <- Cp1 - Cp2
        Lp1 <- mean_distance(g1)
        Lp2 <- mean_distance(g2)
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

      # "Lobe" level
      #-----------------------------------
      } else if (level == 'lobe') {
        t1 <- as.data.table(count_interlobar(g1, 'Temporal'))
        t2 <- as.data.table(count_interlobar(g2, 'Temporal'))
        tdiff <- t1 - t2
        tmp <- data.table(density=density, diff=tdiff)
      }
    }

    tmp
  }
  return(out)
}
