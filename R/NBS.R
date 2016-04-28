#' Network-based statistic for brain MRI data
#'
#' Calculates the \emph{network-based statistic (NBS)}, which allows for
#' family-wise error (FWE) control over network data, introduced for brain MRI
#' data by Zalesky et al. Accepts a three-dimensional array of all subjects'
#' connectivity matrices and a \code{data.table} of covariates, and creates a
#' null distribution of the largest connected component size by permuting
#' subjects across groups. The covariates \code{data.table} must have (at least)
#' a \emph{Group} column.
#'
#' The graph that is returned by this function will have a \code{t.stat} edge
#' attribute which is the t-statistic for that particular connection, along with
#' a \code{p} edge attribute, which is the p-value for that connection.
#' Additionally, each vertex will have a \code{p.nbs} attribute representing
#' \eqn{1 - } the p-value associated with that vertex's component.
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
#' \item{p.init}{Numeric; the initial p-value threshold used}
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
  i <- k <- value <- Var1 <- Var2 <- Var3 <- se <- p <- NULL
  if (!'Group' %in% names(covars)) stop('"covars" must have a "Group" column')
  if (dim(A)[3] != nrow(covars)) stop('"covars" length doesn\'t match array length')

  alt <- match.arg(alternative)
  if (alt == 'two.sided') {
    pfun <- function(q, df) return(2 * (pt(abs(q), df=df, lower.tail=F)))
  } else if (alt == 'less') {
    pfun <- function(q, df) return(pt(q, df=df, lower.tail=F))
  } else if (alt == 'greater') {
    pfun <- function(q, df) return(pt(q, df=df))
  }
  Nv <- nrow(A)
  X <- cbind(1, as.matrix(covars[, lapply(.SD, as.numeric), .SDcols=2:ncol(covars)]))
  df <- nrow(X) - ncol(X)
  z <- which(names(covars) == 'Group')

  # Calculate initial p-values; threshold and create a graph
  #---------------------------------------------------------
  A.m <- setDT(melt(A))
  if (isTRUE(symmetric)) {
    inds.upper <- as.data.table(which(upper.tri(A[, , 1]), arr.ind=TRUE))
    setnames(inds.upper, c('Var1', 'Var2'))
    setkey(inds.upper, Var1, Var2)
    setkey(A.m, Var1, Var2)
    A.m <- A.m[inds.upper]
    setkey(A.m, Var3)
  }
  A.m.sub <- A.m[A.m[, sum(value) > 0, by=list(Var1, Var2)]$V1]
  A.m.sub[, t := with(fastLmPure(X, value, method=2), coefficients / se)[z], by=list(Var1, Var2)]
  A.m.sub[abs(t) > 1.5, p := pfun(t, df)]
  T.dt <- A.m.sub[p < p.init, list(t=unique(t), p=unique(p)), by=list(Var1, Var2)]
  T.mat <- p.mat <- matrix(0, Nv, Nv)
  for (i in seq_len(nrow(T.dt))) {
    T.mat[T.dt$Var1[i], T.dt$Var2[i]] <- T.dt$t[i]
    p.mat[T.dt$Var1[i], T.dt$Var2[i]] <- T.dt$p[i]
  }
  inds.transpose <- which(abs(T.mat) > t(abs(T.mat)), arr.ind=TRUE)

  T.max <- ifelse(abs(T.mat) > t(abs(T.mat)), T.mat, t(T.mat))
  for (i in seq_len(nrow(inds.transpose))) {
    p.mat[inds.transpose[i, 2], inds.transpose[i, 1]] <- p.mat[inds.transpose[i, 1], inds.transpose[i, 2]]
  }
  g.nbs <- graph_from_adjacency_matrix(T.max, diag=F, mode='undirected', weighted=TRUE)
  E(g.nbs)$t.stat <- E(g.nbs)$weight
  E(g.nbs)$p <- 1 - E(graph_from_adjacency_matrix(p.mat, diag=F, mode='undirected', weighted=TRUE))$weight
  if (any(E(g.nbs)$weight < 0)) g.nbs <- delete_edge_attr(g.nbs, 'weight')
  clusts <- components(g.nbs)
  comps <- sort(unique(clusts$csize), decreasing=T)

  # Create a null distribution of maximum component sizes
  #---------------------------------------------------------
  myPerms.nbs <- shuffleSet(n=nrow(covars), nset=N)
  comps.perm <- foreach (k=seq_len(N), .combine='c') %dopar% {
    X.tmp <- as.matrix(X[myPerms.nbs[k, ], ])
    X.tmp[, 'Group'] <- X[, 'Group']
    A.m.tmp <- setDT(melt(A[, , myPerms.nbs[k, ]]))
    if (isTRUE(symmetric)) {
      setkey(A.m.tmp, Var1, Var2)
      A.m.tmp <- A.m.tmp[inds.upper]
      setkey(A.m.tmp, Var3)
    }
    A.m.tmp.sub <- A.m.tmp[A.m.tmp[, sum(value) > 0, by=list(Var1, Var2)]$V1]
    A.m.tmp.sub[, t := with(fastLmPure(X, value, method=2), coefficients / se)[z], by=list(Var1, Var2)]
    A.m.tmp.sub[abs(t) > 1.5, p := ifelse(pfun(t, df) < p.init, 1, 0)]
    T.dt.tmp <- A.m.tmp.sub[p == 1, list(t=unique(t), p=unique(p)), by=list(Var1, Var2)]
    T.mat.tmp <- matrix(0, Nv, Nv)
    for (i in seq_len(nrow(T.dt.tmp))) T.mat.tmp[T.dt.tmp$Var1[i], T.dt.tmp$Var2[i]] <- T.dt.tmp$t[i]
    T.max.tmp <- ifelse(abs(T.mat.tmp) > t(abs(T.mat.tmp)), T.mat.tmp, t(T.mat.tmp))
    max(components(graph_from_adjacency_matrix(T.max.tmp, diag=F, mode='undirected', weighted=TRUE))$csize)
  }
  p.perm <- sapply(comps, function(x) (sum(x <= comps.perm) + 1) / (N + 1))
  x <- clusts$membership
  V(g.nbs)$comp <- match(x, order(table(x), decreasing=TRUE))
  V(g.nbs)$p.nbs <- NA
  for (i in seq_along(comps)) V(g.nbs)[V(g.nbs)$comp == i]$p.nbs <- 1 - p.perm[i]

  return(list(g.nbs=g.nbs, obs=comps, perm=comps.perm, p.perm=p.perm, p.init=p.init))
}
