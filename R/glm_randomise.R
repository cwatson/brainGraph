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
#' @inheritParams GLM
#' @importFrom MASS Null ginv
#' @keywords internal
#'
#' @return A list containing:
#'   \item{X}{Numeric matrix for the covariates of interest}
#'   \item{Z}{Numeric matrix for the nuisance covariates}
#'   \item{eCm}{The \emph{effective contrast}, equivalent to the original, for
#'     the partitioned model \code{[X, Z]} and considering all covariates}
#'   \item{eCx}{Same as \code{eCm}, but considering only \code{X}}

partition <- function(M, con.mat, part.method=c('beckmann', 'guttman', 'ridgway')) {
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

  } else if (part.method == 'ridgway') {
    rZ <- qr(M)$rank - qr(con.mat)$rank
    pinvC <- MASS::ginv(con.mat)
    C0 <- diag(ncol(M)) - t(con.mat) %*% MASS::ginv(t(con.mat))
    tmpX <- M %*% pinvC
    tmpZ <- svd(M %*% C0)$u
    Z <- tmpZ[, 1:rZ]
    X <- tmpX - Z %*% MASS::ginv(Z) %*% tmpX
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
#'   \item Calculate the \dQuote{hat} and residual-forming matrices due to nuisance
#'     alone
#'   \item For F contrasts, return an inverse of the cross product between the
#'     contrast and the inverted design matrix, and the contrast's rank
#'
#' @inheritParams GLM
#' @inheritParams randomise
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
#'   \item{rkC}{(only for F-contrasts) Rank of the effective contrast matrix.}
#'   \item{Xp}{(only for Smith method) List of matrices of covariates of interest}
#'   \item{Zp}{(only for Smith method) List of matrices of nuisance covariates}

setup_randomise <- function(perm.method, part.method, X, contrasts, con.type, nC) {
  n <- nrow(X)
  Mp <- MtM <- Rz <- eC <- Xp <- Zp <- CMtM <- vector('list', length=nC)
  rkC <- rep(0L, nC)

  if (con.type == 't') contrasts <- matrix2list(contrasts)
  for (j in seq_len(nC)) {
    parts <- partition(X, contrasts[[j]], part.method)
    Mp[[j]] <- with(parts, cbind(X, Z))
    MtM[[j]] <- solve(crossprod(Mp[[j]]))

    # The hat and residual-forming matrices are different for diff. methods
    if (perm.method == 'freedmanLane') {
      Hz <- with(parts, Z %*% solve(crossprod(Z)) %*% t(Z))

    } else if (perm.method == 'terBraak') {
      # NOTE: this is really "Hm", and "Rz" is really "Rm"
      Hz <- Mp[[j]] %*% MtM[[j]] %*% t(Mp[[j]])

    } else if (perm.method == 'smith') {
      Xp[[j]] <- parts$X; Zp[[j]] <- parts$Z
      Hz <- Zp[[j]] %*% solve(crossprod(Zp[[j]])) %*% t(Zp[[j]])
    }

    Rz[[j]] <- diag(n) - Hz
    eC[[j]] <- parts$eCm
    if (con.type == 'f') {
      CMtM[[j]] <- solve(eC[[j]] %*% MtM[[j]] %*% t(eC[[j]]))
      rkC[j] <- qr(eC[[j]])$rank
    }
  }
  dfR <- nrow(Mp[[1]]) - qr(Mp[[1]])$rank

  out <- list(Mp=Mp, Rz=Rz, MtM=MtM, eC=eC, dfR=dfR)
  if (con.type == 'f') out <- c(out, list(CMtM=CMtM, rkC=rkC))
  if (perm.method == 'smith') out <- c(out, list(Xp=Xp, Zp=Zp))
  return(out)
}

#' Randomize and fit a model and find the maximum statistic
#'
#' @param DT \code{data.table} with outcome variables
#' @param nC Integer; the number of contrasts
#' @inheritParams GLM
#' @return A \code{data.table} with \code{N} rows and columns specifying the
#'   region, permutation number, effect size (either \code{gamma} or the
#'   numerator for the \emph{F-statistic}), and standard error
#' @keywords internal

randomise <- function(perm.method, part.method, ctype, N, perms, DT, nC, outcome, X, contrasts) {
  i <- region <- numer <- se <- perm <- NULL
  randMats <- setup_randomise(perm.method, part.method, X, contrasts, ctype, nC)
  Mp <- randMats$Mp; Rz <- randMats$Rz; MtM <- randMats$MtM; eC <- randMats$eC; dfR <- randMats$dfR
  if (ctype == 'f') {CMtM <- randMats$CMtM; rkC <- randMats$rkC}

  null.dist <- vector('list', length=nC)
  perm.order <- rep(seq_len(N), each=DT[, length(unique(region))])

  if (perm.method == 'smith') Mp <- lapply(seq_len(nC), function(x) Rz[[x]] %*% randMats$Xp[[x]])
  #-----------------------------------------------------------------------------
  # Loop through contrasts
  #-----------------------------------------------------------------------------
  for (j in seq_len(nC)) {

    # T-contrasts
    #-------------------------------------------------------
    if (ctype == 't') {
      if (perm.method == 'smith') {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          M <- cbind(Mp[[j]][perms[i, ], ], randMats$Zp[[j]])
          MtM <- solve(crossprod(M))
          DT[, brainGraph_GLM_fit_t(M, get(outcome), MtM, eC[[j]]), by=region]
        }
      } else {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          DT[, brainGraph_GLM_fit_t(Mp[[j]], Rz[[j]][perms[i, ], ] %*% get(outcome), MtM[[j]], eC[[j]]), by=region]
        }
      }

    # F-contrasts
    #-------------------------------------------------------
    } else if (ctype == 'f') {
      if (perm.method == 'smith') {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          M <- cbind(Mp[[j]][perms[i, ], ], randMats$Zp[[j]])
          MtM <- solve(crossprod(M))
          CMtM <- solve(eC[[j]] %*% MtM %*% t(eC[[j]]))
          DT[, brainGraph_GLM_fit_f(M, get(outcome), dfR, eC[[j]], rkC[j], CMtM), by=region]
        }
      } else {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          DT[, brainGraph_GLM_fit_f(Mp[[j]], Rz[[j]][perms[i, ], ] %*% get(outcome), dfR, eC[[j]], rkC[j], CMtM[[j]]), by=region]
        }
      }
    }
    null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
  }
  null.dist <- rbindlist(null.dist, idcol='contrast')
  return(null.dist)
}

#' Randomize and fit a model for NBS and return the maximum statistic
#'
#' @param skip Integer vector of the contrast(s) to skip
#' @inheritParams NBS
#' @rdname randomise
#' @keywords internal

randomise_nbs <- function(perm.method, part.method, ctype, N, perms, DT, nC, skip, p.init, X, contrasts, alternative, Nv) {
  se <- perm <- Var1 <- Var2 <- i <- value <- stat <- numer <- NULL
  randMats <- setup_randomise(perm.method, part.method, X, contrasts, ctype, nC)
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

  #-----------------------------------------------------------------------------
  # Loop through contrasts
  #-----------------------------------------------------------------------------
  for (j in seq_len(nC)) {
    if (j %in% skip) next

    # T-contrasts
    #-------------------------------------------------------
    if (ctype == 't') {
      if (perm.method == 'smith') {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          M <- cbind(Rz[[j]][perms[i, ], ] %*% randMats$Xp[[j]], randMats$Zp[[j]])
          MtM <- solve(crossprod(M))
          DT[, brainGraph_GLM_fit_t(M, value, MtM, eC[[j]]), by=list(Var1, Var2)]
        }
      } else {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          DT[, brainGraph_GLM_fit_t(Mp[[j]], Rz[[j]][perms[i, ], ] %*% value, MtM[[j]], eC[[j]]), by=list(Var1, Var2)]
        }
      }
      null.dist[[j]][, stat := gamma / se]
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][statfun(stat, dfR), list(Var1, Var2, stat, perm)]

    # F-contrasts
    #-------------------------------------------------------
    } else {
      if (perm.method == 'smith') {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          M <- cbind(Rz[[j]][perms[i, ], ] %*% randMats$Xp[[j]], randMats$Zp[[j]])
          MtM <- solve(crossprod(M))
          CMtM <- solve(eC[[j]] %*% MtM %*% t(eC[[j]]))
          DT[, brainGraph_GLM_fit_f(M, value, dfR, eC[[j]], rkC[j], CMtM), by=list(Var1, Var2)]
        }
      } else {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          DT[, brainGraph_GLM_fit_f(Mp[[j]], Rz[[j]][perms[i, ], ] %*% value, dfR, eC[[j]], rkC[j], CMtM[[j]]), by=list(Var1, Var2)]
        }
      }
      null.dist[[j]][, stat := numer / (se / dfR)]
      null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
      null.dist[[j]] <- null.dist[[j]][statfun(stat, rkC[j], dfR), list(Var1, Var2, stat, perm)]
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
