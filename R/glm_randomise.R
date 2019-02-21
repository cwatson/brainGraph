#' Partition a design matrix into columns of interest and nuisance
#'
#' Consider the matrix formulation of the \emph{general linear model}:
#' \deqn{\mathbf{Y} = \mathbf{M} \psi + \in}
#' where \eqn{Y} is the vector of outcomes, \eqn{M} is the full design matrix
#' (including nuisance covariates), \eqn{\psi} is the vector of parameter
#' estimates, and \eqn{\in} is the vector of error terms. In a permutation
#' framework, algorithms are applied differently depending on the
#' presence/absence of nuisance covariates; thus the model is separated
#' depending on the contrast of interest:
#' \deqn{\mathbf{Y} = \mathbf{X}\beta + \mathbf{Z}\gamma + \in}
#' where \eqn{\mathbf{X}} contains covariates of interest, \eqn{\mathbf{Z}}
#' contains nuisance covariates, and \eqn{\beta} and \eqn{\gamma} are the
#' associated parameter estimates.
#'
#' @param M Numeric matrix; the full design matrix
#' @param con.mat Numeric matrix; the contrast matrix
#' @param part.method Character string; the method of partitioning (default:
#'   \code{beckmann})
#' @importFrom MASS Null
#' @keywords internal
#'
#' @return A list containing:
#'   \item{X}{Numeric matrix for the covariates of interest}
#'   \item{Z}{Numeric matrix for the nuisance covariates}
#'   \item{eCm}{The \emph{effective contrast}, equivalent to the original, for
#'     the partitioned model \code{[X, Z]} and considering all covariates}
#'   \item{eCx}{Same as \code{eCx}, but considering only \code{X}}
#' @references Guttman I. Linear Models: An Introduction. Wiley, New York, 1982.
#' @references Smith SM, Jenkinson M, Beckmann C, Miller K, Woolrich M (2007).
#'   \emph{Meaningful design and contrast estimability in fMRI.} NeuroImage,
#'   34(1):127-36. \url{https://dx.doi.org/10.1016/j.neuroimage.2006.09.019}

partition <- function(M, con.mat, part.method=c('beckmann', 'guttman')) {
  part.method <- match.arg(part.method)
  if (part.method == 'guttman') {
    idx <- which(con.mat != 0, arr.ind=TRUE)[, 2]
    X <- M[, idx, drop=FALSE]
    Z <- M[, -idx, drop=FALSE]
    eCm <- cbind(con.mat[, idx, drop=FALSE], con.mat[, -idx, drop=FALSE])
  } else if (part.method == 'beckmann') {
    Q <- solve(crossprod(M))
    cdc <- solve(con.mat %*% Q %*% t(con.mat))
    X <- M %*% Q %*% t(con.mat) %*% cdc

    Cu <- MASS::Null(t(con.mat))
    Cv <- Cu - (t(con.mat) %*% cdc %*% con.mat %*% Q %*% Cu)
    Z <- M %*% Q %*% Cv %*% solve(t(Cv) %*% Q %*% Cv)
    eCm <- cbind(diag(ncol(X)), matrix(0, nrow=ncol(X), ncol=ncol(Z)))
  }
  eCx <- eCm[, 1:ncol(X), drop=FALSE]
  return(list(X=X, Z=Z, eCm=eCm, eCx=eCx))
}

#' Helper function to setup for randomise
#'
#' \code{setup_randomise} is used to setup the data/objects for any function
#' that does permutations for GLM-based analysis.
#'
#' The tasks performed by this function are:
#' \enumerate{
#'   \item Separate the design matrix into the independent variables of interest
#'     and nuisance variables (which is contrast-dependent) using
#'     \code{\link{partition}}.
#'   \item Calculate the new contrast(s) based on this new design matrix (also
#'     through \code{\link{partition}}
#'   \item Calculate the inverse of the cross produce of the full model
#'   \item Calculate the "hat" and residual-forming matrices due to nuisance
#'     alone
#'   \item For F contrasts, return an inverse of the cross product between the
#'     contrast and the inverted design matrix, and the contrast's rank
#'
#' @param nC Integer; the number of contrasts
#' @inheritParams GLM
#' @keywords internal
#' @rdname randomise
#' @return A list containing:
#'   \item{Mp}{The full partitioned model, joined}
#'   \item{Rz}{The residual-forming matrix}
#'   \item{MtM}{The inverse of the cross product of the full model}
#'   \item{eC}{The \emph{effective contrast}, equivalent to the original, for
#'     the partitioned model \code{[X, Z]} and considering all covariates}
#'   \item{dfR}{The residual degrees of freedom of the full partitioned model}
#'   \item{CMtM}{(only for F-contrasts) The effective contrast multiplied by the
#'     inverse of the cross-product of the full model.}
#'   \item{rkC}{(only for F-contrasts) The rank of the effective contrast
#'     matrix.}

setup_randomise <- function(X, con.mat, con.type, nC) {
  n <- nrow(X)
  Mp <- MtM <- Rz <- eC <- vector('list', length=nC)

  for (j in seq_len(nC)) {
    if (con.type == 'f') {
      parts <- partition(X, con.mat, 'beckmann')
    } else {
      parts <- partition(X, con.mat[j, , drop=FALSE], 'beckmann')
    }
    Mp[[j]] <- with(parts, cbind(X, Z))
    MtM[[j]] <- solve(crossprod(Mp[[j]]))
    Hz <- with(parts, Z %*% solve(crossprod(Z)) %*% t(Z))
    Rz[[j]] <- diag(n) - Hz
    eC[[j]] <- parts$eCm
  }
  dfR <- nrow(Mp[[1]]) - qr(Mp[[1]])$rank

  out <- list(Mp=Mp, Rz=Rz, MtM=MtM, eC=eC, dfR=dfR)
  if (con.type == 'f') {
    CMtM <- solve(eC[[1]] %*% MtM[[1]] %*% t(eC[[1]]))
    rkC <- qr(eC[[1]])$rank
    out <- c(out, list(CMtM=CMtM, rkC=rkC))
  }
  return(out)
}

#' Randomize and fit a model and find the maximum statistic
#'
#' @param DT \code{data.table} with outcome variables
#' @inheritParams GLM
#' @keywords internal

randomise <- function(ctype, N, perms, DT, nC, measure, X, con.mat, alternative) {
  i <- region <- numer <- se <- perm <- NULL
  randMats <- setup_randomise(X, con.mat, ctype, nC)
  Mp <- randMats$Mp; Rz <- randMats$Rz; MtM <- randMats$MtM; eC <- randMats$eC; dfR <- randMats$dfR
  if (ctype == 'f') {CMtM <- randMats$CMtM; rkC <- randMats$rkC}

  maxfun <- switch(alternative,
                   two.sided=function(x) max(abs(x), na.rm=TRUE),
                   less=function(x) min(x, na.rm=TRUE),
                   greater=function(x) max(x, na.rm=TRUE))
  null.dist <- vector('list', length=nC)
  perm.order <- rep(seq_len(N), each=DT[, length(unique(region))])

  for (j in seq_len(nC)) {
    # T-contrasts
    if (ctype == 't') {
      null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        DT[, brainGraph_GLM_fit_t(Mp[[j]], Rz[[j]][perms[i, ], ] %*% get(measure), MtM[[j]], eC[[j]]), by=region]
      }
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][, maxfun(gamma / se), by=perm][, !'perm']

    # F-contrasts
    } else if (ctype == 'f') {
      null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        DT[, brainGraph_GLM_fit_f(Mp[[j]], Rz[[j]][perms[i, ], ] %*% get(measure), dfR, eC[[j]], rkC, CMtM), by=region]
      }
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][, max(numer / (se / dfR), na.rm=TRUE), by=perm][, !'perm']
    }
  }

  null.dist <- rbindlist(null.dist, idcol='contrast')
  return(null.dist)
}

randomise_nbs <- function(ctype, N, perms, DT, nC, skip, p.init, alternative, Nv) {
  se <- perm <- Var1 <- Var2 <- i <- value <- stat <- numer <- NULL
  randMats <- setup_randomise(X, con.mat, ctype, nC)
  Mp <- randMats$Mp; Rz <- randMats$Rz; MtM <- randMats$MtM; eC <- randMats$eC; dfR <- randMats$dfR

  if (ctype == 't') {
    statfun <- switch(alternative,
                      two.sided=function(stat, df) abs(stat) > qt(p.init / 2, df, lower.tail=FALSE),
                      less=function(stat, df) stat < qt(p.init, df),
                      greater=function(stat, df) stat > qt(p.init, df, lower.tail=FALSE))
  } else {
    statfun <- function(stat, dfN, dfD) stat > qf(p.init / 2, dfN, dfD, lower.tail=FALSE)
    CMtM <- randMats$CMtM; rkC <- randMats$rkC
  }
  maxfun.mat <- switch(alternative,
                       two.sided=function(mat) ifelse(abs(mat) > t(abs(mat)), mat, t(mat)),
                       less=function(mat) pmin(mat, t(mat)),
                       greater=function(mat) pmax(mat, t(mat)))
  null.dist <- comps.perm <- vector('list', length=nC)
  perm.order <- rep(seq_len(N), each=DT[, length(unique(interaction(Var1, Var2)))])

  for (j in seq_len(nC)) {
    if (j %in% skip) next
    # T-contrasts
    if (ctype == 't') {
      null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        DT[, brainGraph_GLM_fit_t(Mp[[j]], Rz[[j]][perms[i, ], ] %*% value, MtM[[j]], eC[[j]]), by=list(Var1, Var2)]
      }
      null.dist[[j]][, stat := gamma / se]
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][statfun(stat, dfR), list(Var1, Var2, stat, perm)]

    # F-contrasts
    } else {
      null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
        DT[, brainGraph_GLM_fit_f(Mp[[j]], Rz[[j]][perms[i, ], ] %*% value, dfR, eC[[j]], rkC, CMtM), by=list(Var1, Var2)]
      }
      null.dist[[j]][, stat := numer / (se / dfR)]
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][statfun(stat, rkC, dfR), list(Var1, Var2, stat, perm)]
    }

    comps.perm[[j]] <- foreach(i=null.dist[[j]][, unique(perm)], .combine='c') %dopar% {
      T.mat.tmp <- matrix(0, Nv, Nv)
      T.mat.tmp[null.dist[[j]][perm == i, cbind(Var1, Var2)]] <- null.dist[[j]][perm == i, stat]
      T.max.tmp <- maxfun.mat(T.mat.tmp)
      max(components(graph_from_adjacency_matrix(T.max.tmp, diag=F, mode='undirected', weighted=TRUE))$csize)
    }
    if (length(comps.perm[[j]]) < N) comps.perm[[j]] <- c(comps.perm[[j]], rep(0, N - length(comps.perm[[j]])))
    comps.perm[[j]] <- data.table(perm=comps.perm[[j]])
  }
  comps.perm <- rbindlist(comps.perm, idcol='contrast')
  return(comps.perm)
}
