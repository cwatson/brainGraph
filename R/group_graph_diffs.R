#' Do between-group tests at each vertex for a given graph measure
#'
#' This function takes two lists of graphs (the length of each equaling the
#' number of subjects per group) and performs either a linear model, 2-sample
#' t-test, or 2-sample Wilcoxon test at each vertex for a given network measure
#' (e.g. vertex degree).
#'
#' The linear model choice is currently pretty inflexible. You will need to
#' provide a data frame of covariates, of which \emph{Study.ID} and \emph{Group}
#' need to be column names. Additionally, all graphs must have a \emph{name}
#' attribute (at the graph level) which matches the \emph{Study.ID} for a given
#' subject. This function will then return the p-value, t-statistic, and
#' parameter estimate related to the \emph{Group} covariate.
#'
#' @param g1 A list of igraph graph objects for group 1
#' @param g2 A list of igraph graph objects for group 2
#' @param measure A character string of the measure to test
#' @param test A character string for the test to use, either 't.test' (default)
#' or 'wilcox.test'
#' @param covars A data frame of covariates; only needed if using \emph{lm}
#' (default: NULL)
#' @param permute Logical; should be \emph{TRUE} if being called from
#' \code{\link{permute.vertex}} (default: FALSE)
#' @param perm.order A vector indicating the order that permuted subjects are
#' in; only necessary if being called from \code{\link{permute.vertex}}
#' @export
#'
#' @return A graph with vertex attributes:
#' \item{size2}{equals the t-statistic}
#' \item{size}{\emph{size2} transformed to be positive values (for
#' visualization purposes)}
#' \item{p}{Equal to \eqn{1-p}}
#' \item{p.fdr}{Equal to \eqn{1-p_{FDR}} (the FDR-adjusted p-value)}
#'
#' @seealso \code{\link[stats]{t.test}, \link[stats]{wilcox.test},
#' \link[stats]{p.adjust}, \link{vec.transform}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' g.diffs.btwn <- group.graph.diffs(g1, g2, 'btwn.cent', test='wilcox.test')
#' g.diffs.btwn <- group.graph.diffs(g.way[[1]][[5]], g.way[[2]][[5]], 'btwn.cent',
#'                                   test='lm', covars=covars.dti)
#' }

group.graph.diffs <- function(g1, g2, measure,
                              test=c('t.test', 'wilcox.test', 'lm'),
                              covars=NULL, permute=FALSE, perm.order=NULL) {
  statistic <- p.value <- Study.ID <- Group <- region <- NULL  # To appease R CHECK

  Nv <- vcount(g1[[1]])
  meas1 <- vapply(g1, vertex_attr, numeric(Nv), measure)
  meas2 <- vapply(g2, vertex_attr, numeric(Nv), measure)

  g.diffs <- g1[[1]]
  test <- match.arg(test)
  if (test == 'lm') {
    # Perform a linear model at each vertex
    #-------------------------------------------------------
    if (is.null(covars)) {
      stop('Please provide a data frame of covariates!')
    }
    id1 <- vapply(g1, graph_attr, character(1), 'name')
    id2 <- vapply(g2, graph_attr, character(1), 'name')
    id.all <- c(id1, id2)
    meas.id.dt <- data.table(Study.ID=id.all)
    meas.id.dt <- cbind(meas.id.dt, t(cbind(meas1, meas2)))
    setnames(meas.id.dt, 2:(nrow(meas1) + 1), V(g.diffs)$name)
    setkey(meas.id.dt, Study.ID)

    if (isTRUE(permute)) {
      if (is.null(perm.order)) {
        stop('You must provide the permuted order of subjects!')
      }
      covars.shuff <- copy(covars)
      covars.shuff <- covars.shuff[perm.order]
      groups <- covars[, levels(Group)]
      covars.shuff[, Group := rep(groups, times=unname(table(covars$Group)))]
      setkey(covars.shuff, Study.ID)
      DT <- merge(covars.shuff, meas.id.dt)
    } else {
      DT <- merge(covars, meas.id.dt)
    }
      DT.tidy <- melt(DT, id.vars=names(covars),
                                   variable.name='region', value.name='measure')
      myLM <- paste('measure ~',
                          paste(c('Group',
                                  names(covars[, !c('Group', 'Study.ID'), with=F])),
                                collapse='+'))
      DT.tidy[, beta := summary(lm(as.formula(myLM), .SD))$coef[2, 1], by=region]
      DT.tidy[, t := summary(lm(as.formula(myLM), .SD))$coef[2, 3], by=region]
      DT.tidy[, p := summary(lm(as.formula(myLM), .SD))$coef[2, 4], by=region]
      p <- DT.tidy[, unique(p), by=region]$V1
      p.fdr <- p.adjust(p, 'fdr')
      DT.tidy[, p.fdr := rep(p.fdr, each=nrow(covars))]
      V(g.diffs)$size2 <- DT.tidy[, unique(t), by=region]$V1
      V(g.diffs)$size <- vec.transform(V(g.diffs)$size2, 0, 20)
      V(g.diffs)$beta <- DT.tidy[, unique(beta), by=region]$V1
      V(g.diffs)$p <- 1 - p
      V(g.diffs)$p.fdr <- 1 - p.fdr
      V(g.diffs)$meas.diff <- rowMeans(meas1) - rowMeans(meas2)

  } else {
    # Do a simple statistical test
    #-------------------------------------------------------
    if (test == 't.test') {
      stats <- Map(function(x, y) t.test(meas1[x, ], meas2[y, ]),
                   seq_len(Nv), seq_len(Nv))
    } else if (test == 'wilcox.test') {
      stats <- Map(function(x, y) wilcox.test(meas1[x, ], meas2[y, ]),
                   seq_len(Nv), seq_len(Nv))
    }

    V(g.diffs)$size2 <- vapply(stats, with, numeric(1), statistic)
    V(g.diffs)$size <- vec.transform(V(g.diffs)$size2, 0, 20)
    V(g.diffs)$p <- 1 - vapply(stats, with, numeric(1), p.value)
    V(g.diffs)$p.fdr <- 1 - p.adjust(1 - V(g.diffs)$p, 'fdr')
    V(g.diffs)$meas.diff <- rowMeans(meas1) - rowMeans(meas2)
  }
  return(g.diffs)
}
