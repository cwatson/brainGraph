#' Run linear models at each vertex of a graph
#'
#' This function takes a list of \code{igraph} graphs and specifies a linear
#' model at each vertex for a given vertex measure (e.g. \emph{degree}).
#' Similarly, you may choose to test for a graph-level measure by specifying
#' \code{level='graph'}.
#'
#' You will need to provide a \code{data.table} of covariates, of which
#' \emph{Study.ID} needs to be the first column. Additionally, all graphs must
#' have a \emph{name} attribute (at the graph level) which matches the
#' \emph{Study.ID} for a given subject. If you do not provide covariates,
#' the code will pull group membership from the graphs' \emph{Group} attributes
#' and do a test of group differences. This function returns, for each region,
#' the contrast of parameter estimates (i.e., \eqn{\gamma}), standard error of
#' the contrast, t-statistic, \eqn{100 (1 - \alpha)}\% confidence interval,
#' p-value, and FDR-adjusted p-value.
#'
#' To test whether a vertex attribute is associated with a different outcome
#' (e.g., \emph{betweenness centrality} and \emph{full-scale IQ}), then specify
#' the relevant outcome variable in the function call, and provide the data in
#' the covariates table.
#'
#' You may optionally do permutation testing by permuting the labels for subject
#' group. This is the same principle as that of Nichols & Holmes (2001) used in
#' voxelwise MRI analyses and implemented in FSL's \emph{randomise}.
#'
#' @param g.list A list of \code{igraph} graph objects for all subjects (if you
#'   have multiple groups, you must concatenate the separate group lists)
#' @param covars A \code{data.table} of covariates
#' @param measure A character string of the vertex measure of interest
#' @param con.vec A numeric vector specifying the contrast of interest
#' @param outcome A character string of the name of the outcome variable; by
#'   default, it is ignored
#' @param X A numeric matrix, if you wish to supply your own design matrix
#'   (default: \code{NULL})
#' @param con.name Character string of the contrast name (default: \code{NULL})
#' @param alternative Character string, whether to do a two- or one-sided test
#'   (default: \code{'two.sided'})
#' @param alpha Numeric; the significance level (default: 0.05)
#' @param level Character string; either \code{vertex} (default) or
#'   \code{graph}
#' @param permute Logical indicating whether or not to permute group labels
#'   (default: \code{FALSE})
#' @param N Integer; number of permutations to create (default: 5e3)
#' @param perms Matrix of permutations, if you would like to provide your own
#'   (default: \code{NULL})
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM_design}}
#' @export
#' @importFrom permute shuffleSet
#'
#' @return A list containing:
#' \item{g}{A graph with vertex attributes: \emph{size2} (t-statistic),
#'   \emph{size} (the t-stat transformed for visualization purposes), \emph{p}
#'   (equal to \eqn{1-p}), \emph{p.fdr} (equal to \eqn{1-p_{FDR}}, the
#'   FDR-adjusted p-value), \emph{gamma} (the contrast of parameter estimaties,
#'   \emph{se} (the standard error of \emph{gamma}); and graph attributes:
#'   \emph{df} (degrees of freedom), \emph{name} (contrast name), \emph{outcome}
#'   (the outcome variable)}
#' \item{DT}{A data table with an entry for each vertex (region)}
#' \item{X}{A numeric matrix; a copy of the \emph{design matrix}}
#' \item{perm}{A list containing: \emph{null.dist} (the null distribution of
#'   maximum t-statistics), \emph{thresh} (the t-statistic value corresponding
#'   to \eqn{100 \times (1 - \alpha)}\% of the null distribution)}
#'
#' @family GLM functions
#' @family Group analysis functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Nichols TE & Holmes AP (2001). \emph{Nonparametric permutation
#'   tests for functional neuroimaging: A primer with examples.} Human Brain
#'   Mapping, 15(1):1-25.
#' @examples
#' \dontrun{
#' g.lm <- brainGraph_GLM(c(g.norm[[1]][[6]], g.norm[[2]][[6]]),
#'   measure='strength',
#'   covars=covars.all[Group == groups[2],
#'                     c(names(covars.dti)[-2], 'Age.Fontan'),
#'                     with=F]
#'   con.vec=c(0, 0, 0, 1))
#'
#' # Test for the group diff between "Neonate" and "Non-Neonate" in nodal eff.
#' covars <- covars.all[Group == groups[2] & tract == 1,
#'   c(names(covars.dti)[-2], 'Age.op1.cat'), with=F]
#' covars[, Age.op1.cat := factor(Age.op1.cat)]
#' g.lm <- brainGraph_GLM(g.norm[[2]][[6]], measure='E.nodal.wt',
#'   covars=covars, con.vec=c(0, 0, 0, 0, 1))
#' }
brainGraph_GLM <- function(g.list, covars, measure, con.vec, outcome=measure,
        X=NULL, con.name=NULL, alternative=c('two.sided', 'less', 'greater'),
        alpha=0.05, level=c('vertex', 'graph'),
        permute=FALSE, N=5e3, perms=NULL, ...) {
  Study.ID <- region <- Outcome <- Covariate <- p.fdr <- p <- Contrast <- i <-
    t.stat <- p.perm <- NULL
  covars <- droplevels(covars)
  incomp <- covars[!complete.cases(covars), Study.ID]

  level <- match.arg(level)
  if (level == 'vertex') {
    A <- t(vapply(g.list, vertex_attr, numeric(vcount(g.list[[1]])), measure))
    DT <- data.table(Study.ID=vapply(g.list, graph_attr, character(1), 'name'))
    DT <- cbind(DT, A)
    setnames(DT, 2:ncol(DT), V(g.list[[1]])$name)
  } else if (level == 'graph') {
    DT <- data.table(Study.ID=vapply(g.list, graph_attr, character(1), 'name'),
                     graph=vapply(g.list, graph_attr, numeric(1), measure))
  }

  alt <- match.arg(alternative)
  DT.cov <- merge(covars, DT, by='Study.ID')

  if (is.null(X)) X <- brainGraph_GLM_design(DT.cov[!Study.ID %in% incomp, names(covars), with=F], ...)
  if (is.vector(con.vec)) con.vec <- t(con.vec)
  stopifnot(ncol(X) == ncol(con.vec))

  DT.m <- melt(DT.cov, id.vars=names(covars), variable.name='region', value.name=measure)
  DT.m <- DT.m[!Study.ID %in% incomp]
  if (outcome != measure) {
    i <- which(colnames(X) == outcome)
    X <- X[, -i]
    DT.lm <- DT.m[, brainGraph_GLM_fit(cbind(X, get(measure)), get(outcome), con.vec, alpha, alt), by=region]
    DT.lm[, Outcome := outcome]
    DT.lm[, Covariate := measure]
  } else {
    DT.lm <- DT.m[, brainGraph_GLM_fit(X, get(measure), con.vec, alpha, alt), by=region]
    DT.lm[, Outcome := measure]
  }
  setnames(DT.lm, 'p.val', 'p')
  DT.lm[, p.fdr :=  p.adjust(p, 'fdr')]
  if (!is.null(con.name)) DT.lm[, Contrast := con.name]

  # Return a graph w/ vertex attributes of statistics
  g.diffs <- make_empty_brainGraph(g.list[[1]]$atlas)
  if (!is.null(con.name)) g.diffs$name <- con.name
  g.diffs$outcome <- measure
  V(g.diffs)$p <- 1 - DT.lm$p
  V(g.diffs)$p.fdr <- 1 - DT.lm$p.fdr
  V(g.diffs)$gamma <- DT.lm$gamma
  V(g.diffs)$se <- DT.lm$se
  V(g.diffs)$size2 <- DT.lm$t.stat
  V(g.diffs)$size <- vec.transform(V(g.diffs)$size2, 0, 20)

  # Permutation testing (similar to FSL's randomise)
  tmax.null.dist <- tmax.thresh <- NA
  if (isTRUE(permute)) {
    if (is.null(perms)) perms <- shuffleSet(n=nrow(X), nset=N)
    groupcol <- grep('Group', colnames(X))
    tmax.obs <- DT.lm[, max(t.stat, na.rm=TRUE)]
    tmax.null.dist <- foreach (i=seq_len(N), .combine='c') %dopar% {
      X.shuff <- X[perms[i, ], ]
      X.shuff[, groupcol] <- X[, groupcol]
      DT.shuff <- DT.cov[perms[i, ]]
      DT.shuff$Group <- DT.cov$Group
      #setkeyv(DT.shuff, key(DT.cov))
      DT.m.shuff <- melt(DT.shuff, id.vars=names(covars),
                         variable.name='region', value.name=measure)
      DT.m.shuff[, brainGraph_GLM_fit(X.shuff, get(measure), con.vec, alpha, alt), by=region][, max(t.stat, na.rm=TRUE)]
    }
    tmax.thresh <- sort(tmax.null.dist)[floor((1 - alpha) * N) + 1]
    DT.lm[, p.perm := (sum(tmax.null.dist >= t.stat, na.rm=TRUE) + 1) / (N + 1), by=region]
    V(g.diffs)$p.perm <- 1 - DT.lm$p.perm
  }

  return(list(g=g.diffs, DT=DT.lm, X=X,
              perm=list(null.dist=tmax.null.dist, alpha=alpha, thresh=tmax.thresh)))
}

#' Fit linear models and calculate statistics
#'
#' This function is partly a wrapper for \code{\link[RcppEigen]{fastLmPure}},
#' fitting a linear model and return the contrast of parameter estimates,
#' standard error, t-statistic, and p-value (given a contrast of interest).
#'
#' For speed purposes (if it is called from \code{\link{brainGraph_GLM}} and
#' permutation testing is done), this function does not do argument checking.
#'
#' @param X Numeric matrix; the "design matrix"
#' @param y Numeric vector of the outcome variable
#' @param con.vec Numeric vector of the contrast of interest
#' @param alpha Numeric; the significance level (default: 0.05)
#' @param alternative Character string for the alternative hypothesis (default:
#'   \code{'two.sided'})
#' @export
#' @importFrom RcppEigen fastLmPure
#'
#' @return A list containing:
#'   \item{gamma}{The contrast of parameter estimates}
#'   \item{se}{The standard error}
#'   \item{t.stat}{The t-statistic}
#'   \item{p.val}{The p-value}
#'   \item{ci.low}{The lower confidence limit}
#'   \item{ci.high}{The upper confidence limit}
#'
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

brainGraph_GLM_fit <- function(X, y, con.vec, alpha=0.05,
                               alternative=c('two.sided', 'less', 'greater')) {
  N <- nrow(X)
  p <- ncol(X) - 1
  df <- N - p - 1
  est <- RcppEigen::fastLmPure(X, y, method=2)
  b <- est$coefficients
  gamma <- as.numeric(con.vec %*% b)
  sigma.squared <- est$s^2
  var.covar <- sigma.squared * solve(crossprod(X))
  #TODO: if `con.vec` is a col vec instead, can use crossprod here too, it's faster
  #se <- as.numeric(sqrt(crossprod(con.vec, var.covar) %*% con.vec))
  #TODO: can also try `tcrossprod`
  #se <- as.numeric(sqrt(con.vec %*% tcrossprod(var.covar, con.vec)))
  se <- as.numeric(sqrt(con.vec %*% var.covar %*% t(con.vec)))
  t.stat <- as.numeric(gamma / se)

  alt <- match.arg(alternative)
  if (alt == 'two.sided') {
    p.val <- 2 * pt(abs(t.stat), df=df, lower.tail=FALSE)
  } else if (alt == 'less') {
    p.val <- pt(t.stat, df=df)
  } else if (alt == 'greater') {
    p.val <- pt(t.stat, df=df, lower.tail=FALSE)
  }
  ci.low <- gamma - qt(alpha / 2, df, lower.tail=F) * se
  ci.high <- gamma + qt(1 - (alpha / 2), df) * se
  list(gamma=gamma, se=se, t.stat=t.stat, p.val=p.val,
       alpha=alpha, ci.low=ci.low, ci.high=ci.high)
}

#' Create a design matrix for linear model analysis
#'
#' This function takes a \code{data.table} of covariates and returns a
#' \emph{design matrix} to be used in linear model analysis.
#'
#' There are three different ways to code factors: \emph{dummy}, \emph{effects},
#' or \emph{cell-means} (chosen by the argument \code{coding}). To understand
#' the difference, see Chapter 7 of the User Guide.
#'
#' The argument \code{mean.center} allows you to mean-center any non-factor
#' variables (including any dummy/indicator covariates). The argument
#' \code{binarize} will convert the given factor variable(s) into numeric
#' variable(s).
#'
#' The \code{int} argument specifies which variables should interact with the
#' \emph{Group} factor variable. This argument accepts either numeric variables
#' (e.g., \emph{Age}) and other factor variables (e.g., \emph{Sex}) if you are
#' running a two-way ANOVA. See Chapter 7 of the User Guide for examples.
#'
#' @param covars A \code{data.table} of covariates
#' @param coding Character string indicating how factor variables will be coded
#'   (default: \code{'dummy'})
#' @param mean.center Logical indicating whether to mean center non-factor
#'   variables (default: \code{FALSE})
#' @param binarize Character string specifying the column name(s) of the
#'   covariate(s) to be converted from type \code{factor} to \code{numeric}
#'   (default: \code{NULL})
#' @param int Character string specifying the column name(s) of the
#'   covariate(s) to test for an interaction (default: \code{NULL})
#' @export
#'
#' @return A numeric matrix
#'
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

brainGraph_GLM_design <- function(covars, coding=c('dummy', 'effects', 'cell.means'),
                                  mean.center=FALSE, binarize=NULL, int=NULL) {

  X <- matrix(1, nrow=nrow(covars), ncol=1)

  if (!is.null(binarize)) {
    for (b in binarize) covars[, eval(b) := as.numeric(get(b)) - 1]
  }

  nums <- which(sapply(covars, is.numeric))
  if (isTRUE(mean.center)) {
    for (n in nums) covars[[n]] <- covars[[n]] - mean(covars[[n]])
  }
  if (length(nums) > 0) X <- cbind(X, as.matrix(covars[, nums, with=F]))

  factors <- which(sapply(covars, class) == 'factor')
  coding <- match.arg(coding)
  for (f in factors) {
    cov.name <- names(covars)[f]
    cov.levels <- covars[, levels(get(cov.name))]
    cov.vec <- covars[, as.numeric(get(cov.name))]

    if (coding == 'cell.means') {
      for (i in 1:max(cov.vec)) {
        cov.vec.sub <- ifelse(cov.vec == i, 1, 0)
        X <- cbind(X, cov.vec.sub)
      }
      if (all(X[, 1] == 1)) X <- X[, -1]  # Remove intercept term
      p <- ncol(X)
      colnames(X)[(p - max(cov.vec) + 1):p] <- paste0(cov.name, cov.levels)

    } else {
      code <- ifelse(coding == 'dummy', 0, -1)
      for (i in 2:max(cov.vec)) {
        cov.vec.sub <- ifelse(cov.vec == i, 1, 0)
        cov.vec.sub <- ifelse(cov.vec == 1, code, cov.vec.sub)
        X <- cbind(X, cov.vec.sub)
      }
      p <- ncol(X)
      colnames(X)[(p - (max(cov.vec) - 1) + 1):p] <- paste0(cov.name, cov.levels[-1])
      colnames(X)[1] <- 'Intercept'
    }
  }

  if (!is.null(int)) {
    intnames <- colnames(X)[grep(int, colnames(X))]
    groupnames <- colnames(X)[grep('Group', colnames(X))]
    if (int %in% names(factors) && coding == 'cell.means') {
      for (i in seq_along(intnames)) {
        X <- cbind(X, X[, intnames[i]] * X[, groupnames])
        p2 <- ncol(X)
        colnames(X)[(p + 1):p2] <- paste(groupnames, intnames[i], sep=':')
        p <- ncol(X)
      }
      X <- X[, -which(colnames(X) %in% intnames)]
      X <- X[, -which(colnames(X) %in% groupnames)]
    } else {
      X <- cbind(X, X[, intnames] * X[, groupnames])
      p2 <- ncol(X)
      colnames(X)[(p + 1):p2] <- paste(groupnames, intnames, sep=':')
    }
    if (!int %in% names(factors)) X <- X[, -which(colnames(X) == int)]
  }

  return(X)
}
