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
#' the \emph{Study.ID} for a given subject. If you do not provide covariates,
#' the code will pull group membership from the graphs' \emph{Group} graph
#' attributes and do a test of group differences. This function will then return
#' the p-value, t-statistic, and parameter estimate associated with the
#' \emph{Group} covariate.
#'
#' If you would like to test whether a vertex attribute is associated with a
#' different outcome variable (for example, \emph{betweenness centrality} and
#' \emph{full-scale IQ}), then specify the relevant outcome variable in the
#' function call, and provide the data in the covariates table. This currently
#' only works for single-group data.
#'
#' You may optionally do permutation testing by permuting the labels for subject
#' group. This is the same principle as that of Nichols & Holmes (2001) used in
#' voxelwise MRI analyses and implemented in FSL's \emph{randomise}.
#'
#' @param g A list of \code{igraph} graph objects for all subjects (if you have
#'   multiple groups, you must concatenate the separate group lists)
#' @param measure A character string of the vertex measure of interest
#' @param outcome A character string of the name of the outcome variable; by
#'   default, it is ignored
#' @param test A character string for the test to use, either 'lm', 't.test', or
#'   'wilcox.test' (default: 'lm')
#' @param alternative Character string, whether to do a two- or one-sided test
#'   (default: 'two.sided')
#' @param covars A \code{data.table} of covariates; needed if using \emph{lm}
#'   (default: NULL)
#' @param permute Logical indicating whether or not to permute group labels
#'   (default: FALSE)
#' @param N Integer; number of permutations to create (default: 5e3)
#' @param alpha Numeric; the significance level (default: 0.05)
#' @export
#' @importFrom RcppEigen fastLmPure
#' @importFrom permute shuffleSet
#'
#' @return A list containing:
#' \item{g}{A graph with vertex attributes: \emph{size2} (t-statistic),
#'   \emph{size} (the t-stat transformed for visualization purposes), \emph{p}
#'   (equal to \eqn{1-p}), \emph{p.fdr} (equal to \eqn{1-p_{FDR}}, the
#'   FDR-adjusted p-value), \emph{beta} (the parameter estimate, if \emph{lm} or
#'   \emph{t.test} is used), \emph{se} (the standard error of \emph{beta}),
#'   \emph{df} (graph-level attribute of the degrees of freedom)}
#' \item{perm}{A list containing: \emph{null.dist} (the null distribution of
#'   maximum t-statistics), \emph{thresh} (the t-statistic value corresponding
#'   to 100\eqn{1 - \alpha}\% of the null distribution)}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Nichols TE & Holmes AP (2001). \emph{Nonparametric permutation
#'   tests for functional neuroimaging: A primer with examples.} Human Brain
#'   Mapping, 15(1):1-25.
#' @examples
#' \dontrun{
#' g.diffs.btwn <- SPM(g=c(g.norm[[1]][[1]], g.norm[[2]][[1]]),
#'   measure='btwn.cent', test='wilcox.test')
#' g.diffs.btwn <- SPM(g=c(g[[1]][[5]], g[[2]][[5]]),
#'   measure='btwn.cent', test='lm', covars=covars.dti)
#' g.corrs.btwn.IQ <- SPM(g=g2, measure='btwn.cent', outcome='IQ', test='lm',
#'   covars=cbind(covars.dti, IQ))
#' }

SPM <- function(g, measure, outcome=measure,
                test=c('lm', 't.test', 'wilcox.test'),
                alternative=c('two.sided', 'less', 'greater'),
                covars=NULL,
                permute=FALSE, N=5e3, alpha=0.05) {
  statistic <- p.value <- Study.ID <- Group <- region <- se <- i <- V1 <- p <- df <- NULL
  tmax.null.dist <- tmax.thresh <- NA

  Nv <- vcount(g[[1]])
  meas <- vapply(g, vertex_attr, numeric(Nv), measure)
  g.diffs <- make_empty_brainGraph(g[[1]])

  id <- vapply(g, graph_attr, character(1), 'name')
  groups <- vapply(g, graph_attr, character(1), 'Group')
  meas.id.dt <- data.table(Study.ID=id, Group=groups)
  setkey(meas.id.dt, Group, Study.ID)
  meas.id.dt <- cbind(meas.id.dt, t(meas))
  setnames(meas.id.dt, 3:(Nv + 2), V(g.diffs)$name)

  alt <- match.arg(alternative)
  test <- match.arg(test)
  # Perform a Wilcoxon test at every vertex
  #---------------------------------------------------------
  if (test == 'wilcox.test') {
    f.wil <- function(x, ...) {
      stats <- wilcox.test(x, ...)
      stats$parameter <- NA
      list(stats$statistic, stats$parameter, stats$p.value)
    }

    DT.tidy <- melt(meas.id.dt, id.vars=c('Study.ID', 'Group'),
                    variable.name='region', value.name=measure)
    DT.tidy[, c('t', 'df', 'p') := f.wil(get(measure) ~ Group, alternative=alt), by=region]
    V(g.diffs)$df <- DT.tidy[, unique(df), by=region]$V1
    V(g.diffs)$meas.diff <- DT.tidy[, mean(get(measure)), by=list(region, Group)][, -diff(V1), by=region]$V1
    V(g.diffs)$p <- 1 - DT.tidy[, unique(p), by=region]$V1

  } else {
  # Specify a linear model at every vertex
  #---------------------------------------------------------
    f <- function(x, y, z) {
      est <- fastLmPure(x, y)
      list(est$coef[z], est$se[z], est$df.resid)
    }
    calc_stats <- function(DT, covars) {
      DT.tidy <- melt(DT, id.vars=names(covars),
                      variable.name='region', value.name=measure)
      setkeyv(DT.tidy, key(DT))
      z <- which(names(DT) == 'Group')
      covars.mat <- cbind(1, as.matrix(DT[, lapply(.SD, as.numeric), .SDcols=names(covars)[-1]]))
      DT.tidy[, c('beta', 'se', 'df') := f(covars.mat, get(measure), z), by=region]
      return(DT.tidy)
    }

    if (test == 'lm' & is.null(covars)) stop('Please provide a data table of covariates!')
    if (test == 't.test') covars <- meas.id.dt[, c('Study.ID', 'Group'), with=F]
    covars[, Group := as.factor(Group)]
    setkeyv(covars, key(meas.id.dt))

    DT <- covars[meas.id.dt]

    # Remove rows with NA if 'outcome' differs from 'measure'
    if (outcome != measure) {
      DT <- DT[complete.cases(DT)]
      covars[, Group := NULL]
      DT[, Group := NULL]
      setkey(DT, Study.ID)
      z <- which(names(DT) == outcome)
      covars.mat <- cbind(1, as.matrix(DT[, lapply(.SD, as.numeric), .SDcols=names(covars)[-c(1, z)]]))
      DT.tidy <- melt(DT, id.vars=names(covars),
                      variable.name='region', value.name='measure')
      setkeyv(DT.tidy, key(DT))
      DT.tidy[, c('beta', 'se', 'df') := f(cbind(covars.mat, measure), get(outcome), z), by=region]
    } else {
      DT.tidy <- calc_stats(DT, covars)
      V(g.diffs)$meas.diff <- DT.tidy[, mean(get(measure)), by=list(region, Group)][, -diff(V1), by=region]$V1
      if (isTRUE(permute)) {
        tmax.observed <- DT.tidy[, max(beta / se)]
        perm.order <- shuffleSet(n=nrow(DT), nset=N)
        tmax.null.dist <- foreach (i=seq_len(N), .combine='c') %dopar% {
          DT.shuff <- DT[perm.order[i, ]]
          DT.shuff$Group <- DT$Group
          setkeyv(DT.shuff, key(DT))
          calc_stats(DT.shuff, DT.shuff[, names(covars), with=F])[, max(beta / se)]
        }
        tmax.thresh <- sort(tmax.null.dist)[floor((1 - alpha) * N) + 1]
      }
    }

    DT.tidy[, t := beta / se, by=region]

    if (alt == 'two.sided') {
      V(g.diffs)$p <- 1 - DT.tidy[, unique(2 * (1 - pt(abs(t), df=df))), by=region]$V1
    } else if (alt == 'less') {
      V(g.diffs)$p <- 1 - DT.tidy[, unique(1 - pt(t, df=df)), by=region]$V1
    } else if (alt == 'greater') {
      V(g.diffs)$p <- 1 - DT.tidy[, unique(pt(t, df=df)), by=region]$V1
    }
    V(g.diffs)$beta <- DT.tidy[, unique(beta), by=region]$V1
    V(g.diffs)$se <- DT.tidy[, unique(se), by=region]$V1
    g.diffs$df <- DT.tidy[, unique(df)]

  }

  V(g.diffs)$p.fdr <- 1 - p.adjust(1 - V(g.diffs)$p, 'fdr')
  V(g.diffs)$size2 <- DT.tidy[, unique(t), by=region]$V1
  V(g.diffs)$size <- vec.transform(V(g.diffs)$size2, 0, 20)
  if (isTRUE(permute)) {
    V(g.diffs)$p.perm <- 1 - vapply(V(g.diffs)$size2, function(x)
                                    sum(tmax.null.dist >= x) / N, numeric(1))
  }

  return(list(g=g.diffs, perm=list(null.dist=tmax.null.dist, thresh=tmax.thresh)))
}
