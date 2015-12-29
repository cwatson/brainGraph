#' Perform between-group tests at each vertex for a given vertex measure
#'
#' This function takes two lists of graphs (the length of each equaling the
#' number of subjects per group) and performs either a linear model, 2-sample
#' t-test, or 2-sample Wilcoxon test at each vertex for a given vertex measure
#' (e.g. \emph{degree}).
#'
#' You will need to provide a \code{data.table} of covariates, of which
#' \emph{Study.ID} and \emph{Group} need to be column names. Additionally, all
#' graphs must have a \emph{name} attribute (at the graph level) which matches
#' the \emph{Study.ID} for a given subject. This function will then return the
#' p-value, t-statistic, and parameter estimate associated with the \emph{Group}
#' covariate.
#'
#' @param g1 A list of \code{igraph} graph objects for group 1
#' @param g2 A list of \code{igraph} graph objects for group 2
#' @param measure A character string of the vertex measure of interest
#' @param test A character string for the test to use, either 't.test',
#' 'wilcox.test', or 'lm' (default: 't.test')
#' @param covars A \code{data.table} of covariates; needed if using \emph{lm}
#' (default: NULL)
#' @param permute Logical; should be \emph{TRUE} if being called from
#' \code{\link{permute.vertex}} (default: FALSE)
#' @param perm.order A vector indicating the order that permuted subjects are
#' in; only necessary if being called from \code{\link{permute.vertex}}
#' @export
#' @importFrom RcppEigen fastLmPure
#'
#' @return A graph with vertex attributes:
#' \item{size2}{equals the t-statistic}
#' \item{size}{\emph{size2} transformed to be positive values (for
#' visualization purposes)}
#' \item{p}{Equal to \eqn{1-p}}
#' \item{p.fdr}{Equal to \eqn{1-p_{FDR}} (the FDR-adjusted p-value)}
#' \item{beta}{Equal to the parameter estimate (if 'lm' is used)}
#'
#' @seealso \code{\link{permute.vertex}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' g.diffs.btwn <- group.graph.diffs(g1, g2, 'btwn.cent', test='wilcox.test')
#' g.diffs.btwn <- group.graph.diffs(g[[1]][[5]], g[[2]][[5]], 'btwn.cent',
#'                                   test='lm', covars=covars.dti)
#' }

group.graph.diffs <- function(g1, g2, measure,
                              test=c('t.test', 'wilcox.test', 'lm'),
                              covars=NULL, permute=FALSE, perm.order=NULL) {
  statistic <- p.value <- Study.ID <- Group <- region <- se <- NULL

  Nv <- vcount(g1[[1]])
  meas1 <- vapply(g1, vertex_attr, numeric(Nv), measure)
  meas2 <- vapply(g2, vertex_attr, numeric(Nv), measure)

  g.diffs <- make_empty_graph(Nv, directed=FALSE)
  g.diffs$atlas <- g1[[1]]$atlas
  V(g.diffs)$name <- V(g1[[1]])$name
  g.diffs <- assign_lobes(g.diffs)

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
    setnames(meas.id.dt, 2:(Nv + 1), V(g.diffs)$name)
    setkey(meas.id.dt, Study.ID)

    if (isTRUE(permute)) {
      if (is.null(perm.order)) {
        stop('You must provide the permuted order of subjects!')
      }
      covars.shuff <- copy(covars)
      covars.shuff <- covars.shuff[perm.order]
      covars.shuff$Group <- covars$Group
      setkey(covars.shuff, Group, Study.ID)
      covars.mat <- as.matrix(covars.shuff[, lapply(.SD, as.numeric), .SDcols=2:ncol(covars)])
      DT <- merge(covars.shuff, meas.id.dt)
    } else {
      DT <- merge(covars, meas.id.dt)
      covars.mat <- as.matrix(covars[, lapply(.SD, as.numeric), .SDcols=2:ncol(covars)])
    }
      DT.tidy <- melt(DT, id.vars=names(covars),
                                   variable.name='region', value.name='measure')
      setkey(DT.tidy, Group, Study.ID)
      z <- which(names(covars) == 'Group')
      DT.tidy[, beta := fastLmPure(cbind(1, covars.mat), measure)$coef[z], by=region]
      DT.tidy[, se := fastLmPure(cbind(1, covars.mat), measure)$se[z], by=region]
      DT.tidy[, df := fastLmPure(cbind(1, covars.mat), measure)$df.resid, by=region]
      DT.tidy[, t := beta / se, by=region]
      DT.tidy[, p := 2 * (1 - pt(abs(t), df=df)), by=region]
      p <- DT.tidy[, unique(p), by=region]$V1
      p.fdr <- p.adjust(p, 'fdr')
      V(g.diffs)$size2 <- DT.tidy[, unique(t), by=region]$V1
      V(g.diffs)$size <- vec.transform(V(g.diffs)$size2, 0, 20)
      V(g.diffs)$beta <- DT.tidy[, unique(beta), by=region]$V1
      V(g.diffs)$p <- 1 - p
      V(g.diffs)$p.fdr <- 1 - p.fdr

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
  }
  V(g.diffs)$meas.diff <- rowMeans(meas1) - rowMeans(meas2)
  return(g.diffs)
}
