#' Permutation test for group difference of graph measures
#'
#' This function draws permutations from linear model residuals to determine the
#' significance of between-group differences of a global or vertex-wise graph
#' measure. This function is intended for cortical thickness networks (in which
#' there is only one graph per group), but can be extended to other types of
#' data. If you would like to calculate differences in the area-under-the-curve
#' (AUC) across densities, then supply a vector to the \code{densities}
#' argument.
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
#' @param permSet Numeric matrix of the set of permutations to loop through; the
#'   number of rows equals the desired number of permutations and the number of
#'   columns equals the total number of subjects across groups
#' @param densities Numeric vector of graph densities
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
#' @family Group analysis functions
#' @seealso \code{\link[igraph]{centr_betw}, \link{vulnerability},
#' \link{count_interlobar}, \link{edge_asymmetry}, \link{efficiency}}
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

permute.group <- function(permSet, densities, resids,
                          level=c('graph', 'vertex', 'lobe', 'other'), atlas,
                          measure=c('btwn.cent', 'degree', 'E.nodal', 'ev.cent',
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
    corrs <- lapply(unique(groups), function(x)
                    lapply(densities, function(y)
                           corr.matrix(as.matrix(resids[which(groups[permSet[i, ]] == x),
                          !c('Study.ID', 'Group'), with=F]),
                          density=y)))
    g <- lapply(corrs, lapply, function(x)
                graph_from_adjacency_matrix(x$r.thresh, mode='undirected', diag=F))

    # Vertex-level
    #-----------------------------------
    if (level == 'vertex') {
      if (measure == 'vulnerability') {
        for (jj in seq_along(g)) {
          for (kk in seq_along(g[jj])) {
            g[[jj]][[kk]]$E.global <- efficiency(g[[jj]][[kk]], 'global')
          }
        }
        meas.list <- lapply(g, function(x) t(sapply(x, vulnerability)))

      } else if (measure == 'degree') {
        meas.list <- lapply(g, function(x) t(sapply(x, degree)))
      } else if (measure == 'E.nodal') {
        meas.list <- lapply(g, function(x) t(sapply(x, efficiency, 'nodal')))
        #TODO: add more
      } else if (measure == 'ev.cent') {
        meas.list <- lapply(g, function(x) t(sapply(x, function(y) centr_eigen(y)$vector)))
      } else if (measure == 'knn') {
        meas.list <- lapply(g, function(x) t(sapply(x, function(y) graph.knn(y)$knn)))
      } else if (measure == 'transitivity') {
        meas.list <- lapply(g, function(x)
                            t(sapply(x, transitivity,type='local', isolates='zero')))
      } else {
        meas.list <- lapply(g, function(x) t(sapply(x, function(y) centr_betw(y)$res)))
      }
      my.diff <- sapply(seq_along(V(g[[1]][[1]])), function(x)
                         auc_diff(densities, cbind(meas.list[[1]][, x], meas.list[[2]][, x])))
      tmp <- as.data.table(t(my.diff))
      setnames(tmp, 1:ncol(tmp), V(g[[1]][[1]])$name)

    # Custom function
    #-----------------------------------
    } else if (level == 'other') {
      tmp <- .function(g, densities)

    } else {
      for (jj in seq_along(g)) {
        for (kk in seq_along(g[[jj]])) {
          g[[jj]][[kk]]$atlas <- atlas
          g[[jj]][[kk]] <- assign_lobes(g[[jj]][[kk]], rand=TRUE)
        }
      }

      # Graph-level
      #-----------------------------------
      if (level == 'graph') {
        mod <- sapply(g, sapply, function(x) modularity(cluster_louvain(x)))
        mod.diff <- auc_diff(densities, mod)
        Cp <- sapply(g, sapply, function(x) transitivity(x, type='localaverage'))
        Cp.diff <- auc_diff(densities, Cp)
        Lp <- sapply(g, sapply, mean_distance)
        Lp.diff <- auc_diff(densities, Lp)
        assort <- sapply(g, sapply, assortativity_degree)
        assort.diff <- auc_diff(densities, assort)
        E.global <- sapply(g, sapply, efficiency, 'global')
        E.global.diff <- auc_diff(densities, E.global)
#        E.local <- sapply(g, sapply, efficiency, 'local', use.parallel=F)
#        E.local.diff <- auc_diff(densities, E.local)

        assort.lobe <- sapply(g, sapply, function(x)
                              assortativity_nominal(x, as.integer(factor(V(x)$lobe))))
        assort.lobe.diff <- auc_diff(densities, assort.lobe)
        asymm <- sapply(g, sapply, function(x) edge_asymmetry(x)$asymm)
        asymm.diff <- auc_diff(densities, asymm)
        tmp <- data.table(mod=mod.diff, E.global=E.global.diff,# E.local=E.local.diff,
                   Cp=Cp.diff, Lp=Lp.diff, assortativity=assort.diff,
                   asymm=asymm.diff, assortativity.lobe=assort.lobe.diff)

      # "Lobe" level
      #-----------------------------------
#      } else if (level == 'lobe') {
#        t1 <- as.data.table(count_interlobar(g1, 'Temporal'))
#        t2 <- as.data.table(count_interlobar(g2, 'Temporal'))
#        tdiff <- t1 - t2
#        tmp <- data.table(density=density, diff=tdiff)
      }
    }

    tmp
  }
  if (length(densities) == 1) out <- cbind(density=densities, out)
  return(out)
}
