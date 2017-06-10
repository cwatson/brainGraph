#' Multi-threshold permutation correction
#'
#' Applies the \emph{multi-threshold permutation correction (MTPC)} method to
#' perform inference in graph theory analyses of brain MRI data.
#'
#' This is a multi-step procedure: (steps 3-4 are the time-consuming steps)
#' \enumerate{
#'   \item Apply thresholds \eqn{\tau} to the networks, and compute network
#'     metrics for all networks and thresholds. (already done beforehand)
#'   \item Compute test statistics \eqn{S_{obs}} for each threshold. (done by
#'     \code{\link{brainGraph_GLM}})
#'   \item Permute group assignments and compute test statistics for each
#'     permutation and threshold. (done by \code{\link{brainGraph_GLM}})
#'   \item Build a null distribution of the maximum statistic across thresholds
#'     (and across brain regions) for each permutation. (done by
#'     \code{\link{brainGraph_GLM}})
#'   \item Determine the critical value, \eqn{S_{crit}} from the null
#'     distribution of maximum statistics.
#'   \item Identify clusters where \eqn{S_{obs} > S_{crit}} and compute the AUC
#'     for these clusters (denoted \eqn{A_{MTPC}}).
#'   \item Compute a critical AUC (\eqn{A_{crit}}) from the mean of the
#'     supra-critical AUC's for the permuted tests.
#'   \item Reject \eqn{H_0} if \eqn{A_{MTPC} > A_{crit}}.
#' }
#'
#' @param g.list A list of lists of \code{igraph} graph objects for all
#'   thresholds and subjects
#' @param thresholds Numeric vector of the thresholds applied to the raw
#'   connectivity matrices.
#' @param N Integer; number of permutations (default: 500)
#' @param perms Matrix of permutations if you would like to provide your own
#'   (default: \code{NULL})
#' @param alpha Numeric; the significance level (default: 0.05)
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM}}
#' @export
#'
#' @return A list with the following elements:
#'   \item{res.glm}{List with length equal to the number of thresholds; each
#'     list element is the output from \code{\link{brainGraph_GLM}}}
#'   \item{DT}{A \code{data.table} for all thresholds, combined from the outputs
#'     of \code{\link{brainGraph_GLM}}}
#'   \item{stats}{A list containing \code{S.mtpc} (the max. observed
#'     statistic), \code{tau.mtpc} (the threshold of the max. observed
#'     statistic), \code{S.crit} (the critical statistic value), and
#'     \code{A.crit} (the critical AUC)}
#'   \item{null.dist}{Numeric matrix with \code{N} rows and number of columns
#'     equal to the number of thresholds. Each element is the maximum statistic
#'     for that permutation and threshold.}
#'   \item{perm.order}{Numeric matrix; the permutation set applied for all
#'     thresholds (each row is a separate permutation)}
#'
#' @family Group analysis functions
#' @family GLM functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Drakesmith M, Caeyenberghs K, Dutt A, Lewis G, David AS, Jones
#'   DK (2015). \emph{Overcoming the effects of false positives and threshold
#'   bias in graph theoretical analyses of neuroimaging data.} NeuroImage,
#'   118:313-333.
#' @examples
#' \dontrun{
#' diffs.mtpc <- mtpc(g.list=g.norm, thresholds=thresholds, N=N,
#'      covars=covars.dti, measure='E.nodal.wt', coding='effects',
#'      con.vec=c(0, 0, 0, 0, -2), alt='greater',
#'      binarize=c('Sex', 'Scanner'), con.name='Group 1 > Group 2')
#' }

mtpc <- function(g.list, thresholds, N=500, perms=NULL, alpha=0.05, ...) {
  A.crit <- A.mtpc <- DT <- n.rand <- region <- t.stat <- threshold <- values <- NULL
  stopifnot(all(lengths(g.list) == length(thresholds)))
  if (is.null(perms)) perms <- shuffleSet(n=sum(vapply(g.list, function(x) length(x[[1]]), numeric(1))), nset=N)

  res.glm <- lapply(seq_along(thresholds), function(x)
                    brainGraph_GLM(c(g.list[[1]][[x]], g.list[[2]][[x]]), N=N,
                                   permute=TRUE, perms=perms, alpha=alpha, ...))
  for (i in seq_along(thresholds)) res.glm[[i]]$DT[, threshold := thresholds[i]]
  mtpc.all <- rbindlist(lapply(res.glm, with, DT))
  mtpc.all[, n.rand := N]

  null.dist.all <- vapply(res.glm, function(x) x$perm$null.dist, numeric(N))
  null.dist.max <- apply(null.dist.all, 1, max)
  S.crit <- sort(null.dist.max)[floor(N * (1 - alpha)) + 1]
  mtpc.all[, S.crit := S.crit]

  crit.regions <- mtpc.all[, rle(t.stat > S.crit), by=region][values == TRUE & lengths > 2, unique(region)]
  mtpc.all[, A.mtpc := 0]
  mtpc.all[region %in% crit.regions, A.mtpc := auc_rle(t.stat, S.crit, thresholds), by=region]

  S.mtpc <- mtpc.all[, max(t.stat, na.rm=TRUE)]
  tau.mtpc <- mtpc.all[t.stat == S.mtpc, threshold]

  null.crit <- which(apply(null.dist.all, 1, function(x)
                           any(with(rle(x > S.crit), values == TRUE & lengths > 2))))
  if (length(null.crit) > 0) {
    if (length(null.crit) > 1) {
      A.crit <- mean(apply(null.dist.all[null.crit, ], 1, auc_rle, S.crit, thresholds))
    } else {
      A.crit <- auc_rle(null.dist.all[null.crit, ], S.crit, thresholds)
    }
    mtpc.all[, A.crit := A.crit]
  } else {
    mtpc.all[, A.crit := 0]
  }

  mtpc.stats <- list(S.mtpc=S.mtpc, tau.mtpc=tau.mtpc,
                     S.crit=S.crit, A.crit=mtpc.all[1, A.crit])
  return(list(res.glm=res.glm, DT=mtpc.all, stats=mtpc.stats,
              null.dist=null.dist.all, perm.order=perms))
}

auc_rle <- function(t.stat, S.crit, thresholds) {
  runs <- rle(t.stat > S.crit)
  myruns <- which(runs$values == TRUE & runs$lengths > 1)
  runs.len.cumsum <- cumsum(runs$lengths)
  ends <- runs.len.cumsum[myruns]
  newindex <- ifelse(myruns > 1, myruns - 1, 0)
  starts <- runs.len.cumsum[newindex] + 1
  if (0 %in% newindex) starts <- c(1, starts)

  x <- Map(function(a, b) thresholds[a:b], starts, ends)
  y <- Map(function(a, b) t.stat[a:b], starts, ends)
  auc <- mapply(auc_diff, x, y)
  return(sum(auc))
}
