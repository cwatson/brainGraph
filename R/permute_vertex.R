#' Permutation test for group difference of graph vertex measures
#'
#' This function performs a permutation test looking at the max t-statistic for
#' a given vertex-level graph attribute (e.g. degree). This is the same
#' principle as that of Nichols & Holmes (2001) used in voxelwise MRI analysis.
#'
#' @param g1 A list of igraph graph objects for group 1
#' @param g2 A list of igraph graph objects for group 2
#' @param alpha The significance level (default: 0.05)
#' @param num.perms The number of permutations (default: 1e3)
#' @param measure A character string of the measure to test
#' @export
#'
#' @return A list with the following elements:
#' \item{g}{A graph representing group differences, with \emph{p.perm} included
#' as a vertex-level attribute (it is actually 1 - \emph{p.perm})}
#' \item{p.max}{The proportion of the number of times the maximum t-statistic of
#' the permuted groups was greater than the observed maximum t-statistic}
#' \item{thresh}{The (1 - alpha) \%ile maximum t-statistic of all permuted values}
#'
#' @seealso \code{\link{sample}, \link{group.graph.diffs}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Nichols TE & Holmes AP (2001). \emph{Nonparametric permutation
#' tests for functional neuroimaging: A primer with examples.} Human Brain
#' Mapping, 15(1):1-25.
#' @examples
#' \dontrun{
#' g.diff <- permute.groups(g.wt[[1]][[3]], g.wt[[2]][[3]], measure='degree')
#' }

permute.vertex <- function(g1, g2, alpha=0.05, num.perms=1e3, measure) {
  combined <- c(g1, g2)
  g.diffs <- group.graph.diffs(g1, g2, measure)
  max.observed <- max(V(g.diffs)$size2)

  n1 <- length(g1)
  n.all <- length(combined)
  max.rand <- foreach(i=seq_len(num.perms), .combine='c') %dopar% {
    shuffled <- sample(n.all)
    g1.rand <- combined[shuffled[1:n1]]
    g2.rand <- combined[shuffled[(n1 + 1):n.all]]

    max(V(group.graph.diffs(g1.rand, g2.rand, measure))$size2)
  }

  p.max <- (sum(abs(max.rand) >= abs(max.observed)) + 1) / (num.perms + 1)
  thresh <- sort(max.rand)[(1 - alpha) * num.perms]
  Nv <- vcount(g.diffs)
  V(g.diffs)$p.perm <- 1 - sapply(seq_len(Nv), function(x)
                                  sum(max.rand >= V(g.diffs)$size2[x]) / num.perms)
  return(list(g=g.diffs, p.max=p.max, thresh=thresh))
}
