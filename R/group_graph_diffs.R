#' Do between-group t-tests at each vertex for a given graph measure
#'
#' This function takes two lists of graphs (the length of each equaling the
#' number of subjects per group) and performs a 2-sample t-test at each vertex
#' for a given network measure (e.g. vertex degree).
#'
#' @param g1 A list of igraph graph objects for group 1
#' @param g2 A list of igraph graph objects for group 2
#' @param measure A character string of the measure to test
#'
#' @export
#'
#' @return A graph with vertex attributes: \emph{size2} equals the t-statistic;
#' \emph{size} is transformed to be positive values for visualization; \emph{p}
#' is equal to \eqn{1 - p}; and \emph{p.adj} is 1 - the FDR-adjusted p-value
#'
#' @seealso \code{\link[stats]{t.test}, \link[stats]{p.adjust},
#' \link{vec.transform}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

group.graph.diffs <- function(g1, g2, measure) { 
  Nv <- vcount(g1[[1]])
  meas1 <- sapply(g1, vertex_attr, measure)
  meas2 <- sapply(g2, vertex_attr, measure)
  t.stats <- mapply(function(x, y) t.test(meas1[x, ], meas2[y, ]),
                    seq_len(Nv), seq_len(Nv), SIMPLIFY=F)

  g.diffs <- g1[[1]]
  V(g.diffs)$size2 <- sapply(t.stats, function(x) x$statistic)
  V(g.diffs)$size <- vec.transform(V(g.diffs)$size2, 0, 20)
  V(g.diffs)$p <- sapply(t.stats, function(x) x$p.value)
  V(g.diffs)$p.adj <- p.adjust(V(g.diffs)$p, 'fdr')
  V(g.diffs)$p <- 1 - V(g.diffs)$p
  V(g.diffs)$p.adj <- 1 - V(g.diffs)$p.adj
  return(g.diffs)
}
