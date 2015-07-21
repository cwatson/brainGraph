#' Do between-group tests at each vertex for a given graph measure
#'
#' This function takes two lists of graphs (the length of each equaling the
#' number of subjects per group) and performs either a 2-sample t-test or
#' 2-sample Wilcoxon test at each vertex for a given network measure (e.g.
#' vertex degree).
#'
#' @param g1 A list of igraph graph objects for group 1
#' @param g2 A list of igraph graph objects for group 2
#' @param measure A character string of the measure to test
#' @param test A character string for the test to use, either 't.test' (default)
#' or 'wilcox.test'
#' @export
#'
#' @return A graph with vertex attributes:
#' \item{size2}{equals the t-statistic}
#' \item{size}{\emph{size2} transformed to be positive values (for
#' visualization purposes)}
#' \item{p}{Equal to \eqn{1-p}}
#' \item{p.adj}{Equal to \eqn{1-p_{FDR}} (the FDR-adjusted p-value)}
#'
#' @seealso \code{\link[stats]{t.test}, \link[stats]{wilcox.test},
#' \link[stats]{p.adjust}, \link{vec.transform}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' g.diffs.btwn <- group.graph.diffs(g1, g2, 'btwn.cent', test='wilcox.test')
#' }

group.graph.diffs <- function(g1, g2, measure, test=c('t.test', 'wilcox.test')) { 
  statistic <- p.value <- NULL  # To appease R CHECK

  Nv <- vcount(g1[[1]])
  meas1 <- vapply(g1, vertex_attr, numeric(Nv), measure)
  meas2 <- vapply(g2, vertex_attr, numeric(Nv), measure)

  test <- match.arg(test)
  if (test == 't.test') {
    stats <- mapply(function(x, y) t.test(meas1[x, ], meas2[y, ]),
                      seq_len(Nv), seq_len(Nv), SIMPLIFY=F)
  } else if (test == 'wilcox.test') {
    stats <- mapply(function(x, y) wilcox.test(meas1[x, ], meas2[y, ]),
                      seq_len(Nv), seq_len(Nv), SIMPLIFY=F)
  }

  g.diffs <- g1[[1]]
  V(g.diffs)$size2 <- vapply(stats, with, numeric(1), statistic)
  V(g.diffs)$size <- vec.transform(V(g.diffs)$size2, 0, 20)
  V(g.diffs)$p <- 1 - vapply(stats, with, numeric(1), p.value)
  V(g.diffs)$p.adj <- 1 - p.adjust(1 - V(g.diffs)$p, 'fdr')
  return(g.diffs)
}
