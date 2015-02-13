#' Permutation test for group difference of global graph measures
#'
#' This function draws permutations from linear model residuals to determine the
#' significance of between-group differences of a global graph measure. This
#' function is intended for cortical thickness networks (in which there is only
#' one graph per group).
#'
#' @param densities A vector of graph densities to loop through
#' @param resids A data.table of the residuals (from \code{\link{get.resid}})
#' @param num.subjs A vector of length 2 indicating group sizes
#' @param num.perms The number of permutations to perform (default: 1e3)
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' out <- permute.global(densities, m$resids, summary(covars$Group), 1e3)
#' }

permute.global <- function(densities, resids, num.subjs, num.perms=1e3) {
  n1 <- num.subjs[1]
  n.all <- sum(num.subjs)

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
    data.table(density=densities, mod=mod.diff, E.global=E.global.diff, Cp=Cp.diff,
               Lp=Lp.diff, assortativity=assort.diff)
  }
  return(out)
}
