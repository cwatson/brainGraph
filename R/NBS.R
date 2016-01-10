#' Network-based statistic for brain MRI data
#'
#' Calculates the \emph{network-based statistic (NBS)}, which allows for
#' family-wise error (FWE) control over network data, introduced for brain MRI
#' data by Zalesky et al. Accepts a three-dimensional array of all subjects'
#' connectivity matrices and a \code{data.table} of covariates, and creates a
#' null distribution of the largest connected component size by permuting
#' subjects. If you would like to perform a t-test at each element, then supply
#' a covariates \code{data.table} with only a \emph{Study.ID} and \emph{Group}
#' column. The graph that is returned by this function will have a
#' \code{t.stat} edge attribute which is the t-statistic for that particular
#' connection.
#'
#' @param A Three-dimensional array of all subjects' connectivity matrices
#' @param covars A \code{data.table} of covariates
#' @param alternative Character string, whether to do a two- or one-sided test
#' (default: 'two.sided')
#' @param p.init Numeric; the initial p-value threshold (default: 0.001)
#' @param N Integer; the number of permutations (default: 1e3)
#' @param symmetric Logical indicating if input matrices are symmetric (default:
#'   FALSE)
#' @export
#' @importFrom RcppEigen fastLmPure
#' @importFrom permute shuffleSet
#'
#' @return A list containing:
#' \item{g.nbs}{The \code{igraph} graph object based on the initial threshold}
#' \item{obs}{Integer vector of the observed connected component sizes}
#' \item{perm}{Integer vector of the permutation distribution of largest
#'   connected component sizes}
#' \item{p.perm}{Numeric vector of the permutation p-values for each component}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Zalesky A., Fornito A., Bullmore E.T. (2010) \emph{Network-based
#'   statistic: identifying differences in brain networks}. NeuroImage,
#'   53(4):1197-1207.
#' @examples
#' \dontrun{
#' max.comp.nbs <- NBS(A.norm.sub[[1]], covars.dti, N=5e3)
#' }

NBS <- function(A, covars, alternative=c('two.sided', 'less', 'greater'),
                p.init=0.001, N=1e3, symmetric=FALSE) {
  i <- k <- NULL
  if (!'Group' %in% names(covars)) stop('"covars" must have a "Group" column')
  alt <- match.arg(alternative)
  if (alt == 'two.sided') {
    pfun <- function(T, df) return(2 * (pt(abs(T), df=df, lower.tail=F)))
  } else if (alt == 'less') {
    pfun <- function(T, df) return(pt(T, df=df, lower.tail=F))
  } else if (alt == 'greater') {
    pfun <- function(T, df) return(pt(T, df=df))
  }
  Nv <- nrow(A)
  covars.mat <- as.matrix(covars[, lapply(.SD, as.numeric), .SDcols=2:ncol(covars)])
  df <- nrow(covars.mat) - ncol(covars.mat) - 1
  z <- which(names(covars) == 'Group')

  # Calculate initial p-values; threshold and create a graph
  #---------------------------------------------------------
  if (isTRUE(symmetric)) {
    inds.upper <- which(upper.tri(A[, , 1]), arr.ind=TRUE)
    p <- matrix(NaN, nrow=Nv, ncol=Nv)
    est.vec <- foreach(i=seq_len(nrow(inds.upper)), .combine='rbind') %dopar% {
      est <- fastLmPure(cbind(1, covars.mat),
                        A[inds.upper[i, 1], inds.upper[i, 2], ], method=2)

      data.table(est$coef[z], est$se[z])
    }
    p.vec <- with(est.vec, pfun(V1 / V2, df))
    p[upper.tri(p, diag=FALSE)] <- p.vec
    p <- ifelse(p < p.init, 1, 0)
    g.nbs <- graph_from_adjacency_matrix(p, diag=F, mode='undirected')
    E(g.nbs)$t.stat <- with(est.vec, V1 / V2)[which(!is.na(p.vec) & p.vec < p.init)]
  } else {
    T.tmp <- coefs <- stderrs <- rep(0, Nv)
    T <- foreach(i=seq_len(Nv), .combine='rbind') %dopar% {
      for (j in seq_len(Nv)) {
        est <- fastLmPure(cbind(1, covars.mat), A[i, j, ], method=2)
        coefs[j] <- est$coef[z]
        stderrs[j] <- est$se[z]
      }
      T.tmp <- coefs / stderrs
      T.tmp
    }
    p <- ifelse(pfun(T, df) < p.init, 1, 0)
    p[is.na(p)] <- T[is.nan(T)] <- 0
    T.max <- ifelse(abs(T) > t(abs(T)), T, t(T))
    p <- ifelse(abs(p) > t(abs(p)), p, t(p))
    p[lower.tri(p)] <- 0
    g.nbs <- graph_from_adjacency_matrix(p, diag=F, mode='undirected')
    E(g.nbs)$t.stat <- T.max[p == 1]
  }
  comps <- unique(components(g.nbs)$csize)

  # Create a null distribution of maximum components sizes
  #---------------------------------------------------------
  myPerms.nbs <- shuffleSet(n=nrow(covars), nset=N)
  if (isTRUE(symmetric)) {
    coefs <- stderrs <- rep(0, length=nrow(inds.upper))
    p <- matrix(NaN, nrow=Nv, ncol=Nv)
    comps.perm <- foreach (k=seq_len(N), .combine='c') %dopar% {
      covars.mat.tmp <- as.matrix(covars.mat[myPerms.nbs[k, ], ])
      covars.mat.tmp[, z-1] <- covars.mat[, z-1]
      for (m in seq_len(nrow(inds.upper))) {
        est <- fastLmPure(cbind(1, covars.mat),
                          A[inds.upper[m, 1], inds.upper[m, 2], ][myPerms.nbs[k, ]],
                          method=2)
        coefs[m] <- est$coef[z]
        stderrs[m] <- est$se[z]
      }
      p.vec <- pfun(coefs / stderrs, df)
      p[upper.tri(p, diag=FALSE)] <- p.vec
      p <- ifelse(p < p.init, 1, 0)
      max(components(graph_from_adjacency_matrix(p, diag=F, mode='undirected'))$csize)
    }
  } else {
    coefs <- stderrs <- matrix(0, nrow=Nv, ncol=Nv)
    comps.perm <- foreach (k=seq_len(N), .combine='c') %dopar% {
      covars.mat.tmp <- as.matrix(covars.mat[myPerms.nbs[k, ], ])
      covars.mat.tmp[, z-1] <- covars.mat[, z-1]
      for (m in seq_len(Nv)) {
        for (n in seq_len(Nv)) {
          est <- fastLmPure(cbind(1, covars.mat.tmp),
                            A[m, n, ][myPerms.nbs[k, ]], method=2)
          coefs[m, n] <- est$coef[z]
          stderrs[m, n] <- est$se[z]
        }
      }
      p.tmp <- pfun(coefs / stderrs, df)
      p.tmp <- ifelse(p.tmp < p.init, 1, 0)
      max(components(graph_from_adjacency_matrix(p.tmp, diag=F, mode='undirected'))$csize)
    }
  }
  p.perm <- sapply(comps, function(x) (sum(x <= comps.perm) + 1) / (N + 1))
  return(list(g.nbs=g.nbs, obs=comps, perm=comps.perm, p.perm=p.perm))
}
