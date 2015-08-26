#' Permutation test for group difference of graph vertex measures
#'
#' This function performs a permutation test looking at the max t-statistic for
#' a given vertex-level graph attribute (e.g. degree). This is the same
#' principle as that of Nichols & Holmes (2001) used in voxelwise MRI analysis.
#'
#' By default, a t-test is used when calling \code{\link{group.graph.diffs}}; if
#' you would like to do e.g. a linear model, then see that function's help
#' section and pass the appropriate arguments.
#'
#' @param g1 A list of igraph graph objects for group 1
#' @param g2 A list of igraph graph objects for group 2
#' @param measure A character string of the measure to test
#' @param test Character string for the test to use (passed to
#' \code{\link{group.graph.diffs}}); one of 't.test', 'wilcox.test', or 'lm'
#' @param alpha The significance level (default: 0.05)
#' @param N The number of permutations (default: 1e3)
#' @param ... Other arguments passed to \code{\link{group.graph.diffs}}
#' @export
#'
#' @return A list with the following elements:
#' \item{g}{A graph representing group differences, with \emph{p.perm} included
#' as a vertex-level attribute (it is actually 1 - \emph{p.perm})}
#' \item{p.max}{The proportion of the number of times the maximum t-statistic of
#' the permuted groups was greater than the observed maximum t-statistic}
#' \item{thresh}{The \eqn{1 - \alpha} \%ile maximum t-statistic of all permuted values}
#'
#' @seealso \code{\link{sample}, \link{group.graph.diffs}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Nichols TE & Holmes AP (2001). \emph{Nonparametric permutation
#' tests for functional neuroimaging: A primer with examples.} Human Brain
#' Mapping, 15(1):1-25.
#' @examples
#' \dontrun{
#' g.diff <- permute.vertex(g.wt[[1]][[3]], g.wt[[2]][[3]], measure='degree')
#' g.diff <- permute.vertex(g[[1]][[5]], g[[2]][[5]], measure='degree',
#'     test='lm', covars=covars.dti)
#' }

permute.vertex <- function(g1, g2, measure, test=c('t.test', 'wilcox.test', 'lm'),
                           alpha=0.05, N=1e3, ...) {
  combined <- c(g1, g2)
  g.diffs <- group.graph.diffs(g1, g2, measure)
  max.observed <- max(V(g.diffs)$size2)

  n1 <- length(g1)
  n.all <- length(combined)
  max.rand <- foreach(i=seq_len(N), .combine='c') %dopar% {
    shuffled <- sample(n.all)
    g1.rand <- combined[shuffled[1:n1]]
    g2.rand <- combined[shuffled[(n1 + 1):n.all]]

    if (test == 'lm') {
      max(V(group.graph.diffs(g1.rand, g2.rand, measure=measure, test=test,
                              permute=TRUE, perm.order=shuffled, ...))$size2,
          na.rm=T)
    } else {
      max(V(group.graph.diffs(g1.rand, g2.rand, measure, test, ...))$size2,
          na.rm=T)
    }
  }

  p.max <- (sum(abs(max.rand) >= abs(max.observed)) + 1) / (N + 1)
  thresh <- sort(max.rand)[(1 - alpha) * N]
  Nv <- vcount(g.diffs)
  V(g.diffs)$p.perm <- 1 - vapply(seq_len(Nv), function(x)
                                  sum(max.rand >= V(g.diffs)$size2[x]) / N,
                                  numeric(1))
  return(list(g=g.diffs, p.max=p.max, thresh=thresh))
}
