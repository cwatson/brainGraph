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
    cdc <- solve(con.mat %*% tcrossprod(Q, con.mat))
    X <- M %*% tcrossprod(Q, con.mat) %*% cdc

    Cu <- MASS::Null(t(con.mat))
    Cv <- Cu - (crossprod(con.mat, cdc) %*% con.mat %*% Q %*% Cu)
    Z <- M %*% Q %*% Cv %*% solve(crossprod(Cv, Q) %*% Cv)
    px <- dim(X)[2L]
    eCm <- cbind(diag(px), matrix(0, nrow=px, ncol=dim(Z)[2L]))

  } else if (part.method == 'ridgway') {
    rZ <- qr(M)$rank - qr(con.mat)$rank
    pinvC <- MASS::ginv(con.mat)
    C0 <- diag(dim(M)[2L]) - t(con.mat) %*% MASS::ginv(t(con.mat))
    tmpX <- M %*% pinvC
    tmpZ <- svd(M %*% C0)$u
    Z <- tmpZ[, 1:rZ]
    X <- tmpX - Z %*% MASS::ginv(Z) %*% tmpX
    px <- dim(X)[2L]
    eCm <- cbind(diag(px), matrix(0, nrow=px, ncol=dim(Z)[2L]))
  }
  eCx <- eCm[, 1:dim(X)[2L], drop=FALSE]
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
#' }
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
  n <- dim(X)[1L]
  Mp <- MtM <- Rz <- eC <- Xp <- Zp <- CMtM <- vector('list', length=nC)
  rkC <- rep(0L, nC)

  if (con.type == 't') contrasts <- matrix2list(contrasts)
  for (j in seq_len(nC)) {
    parts <- partition(X, contrasts[[j]], part.method)
    Mp[[j]] <- with(parts, cbind(X, Z))
    MtM[[j]] <- solve(crossprod(Mp[[j]]))

    # The hat and residual-forming matrices are different for diff. methods
    if (perm.method == 'freedmanLane') {
      Hz <- with(parts, Z %*% tcrossprod(solve(crossprod(Z)), Z))

    } else if (perm.method == 'terBraak') {
      # NOTE: this is really "Hm", and "Rz" is really "Rm"
      Hz <- Mp[[j]] %*% tcrossprod(MtM[[j]], Mp[[j]])

    } else if (perm.method == 'smith') {
      Xp[[j]] <- parts$X; Zp[[j]] <- parts$Z
      Hz <- Zp[[j]] %*% tcrossprod(solve(crossprod(Zp[[j]])), Zp[[j]])
    }

    Rz[[j]] <- diag(n) - Hz
    eC[[j]] <- parts$eCm
    if (con.type == 'f') {
      CMtM[[j]] <- solve(eC[[j]] %*% tcrossprod(MtM[[j]], eC[[j]]))
      rkC[j] <- qr(eC[[j]])$rank
    }
  }
  dfR <- dim(Mp[[1]])[1L] - qr(Mp[[1]])$rank

  out <- list(Mp=Mp, Rz=Rz, MtM=MtM, eC=eC, dfR=dfR)
  if (con.type == 'f') out <- c(out, list(CMtM=CMtM, rkC=rkC))
  if (perm.method == 'smith') out <- c(out, list(Xp=Xp, Zp=Zp))
  return(out)
}

#' Randomize and fit a model and find the maximum statistic
#'
#' @param nC Integer; the number of contrasts
#' @param skip Integer vector indicating which (if any) contrasts to skip (if
#'   there were no significant differences anywhere)
#' @param DT \code{data.table} with outcome variables
#' @param mykey The \code{key} to key by (to differentiate NBS and other GLM
#'   analyses). For GLM, it is \code{'region'}; for NBS, it is
#' @inheritParams GLM
#' @return A \code{data.table} with \code{N} rows and columns specifying the
#'   region, permutation number, effect size (either \code{gamma} or the
#'   numerator for the \emph{F-statistic}), and standard error
#' @keywords internal

randomise <- function(perm.method, part.method, N, perms,
                      contrasts, con.type, nC, skip, DT, outcome, X, mykey) {
  i <- region <- numer <- se <- perm <- Var1 <- Var2 <- value <- stat <- NULL
  randMats <- setup_randomise(perm.method, part.method, X, contrasts, con.type, nC)
  Mp <- randMats$Mp; Rz <- randMats$Rz; MtM <- randMats$MtM; eC <- randMats$eC; dfR <- randMats$dfR
  if (con.type == 'f') {CMtM <- randMats$CMtM; rkC <- randMats$rkC}

  null.dist <- vector('list', length=nC)
  if (grepl(',', mykey)) {
    v <- DT[, eval(parse(text=paste0('interaction(', mykey, ')')))]
    perm.order <- rep(seq_len(N), each=length(unique(v)))
    mykey <- parse(text=paste0('list(', paste0(strsplit(mykey, ',')[[1]], collapse=','), ')'))
  } else {
    perm.order <- rep(seq_len(N), each=DT[, length(unique(get(mykey)))])
  }

  if (perm.method == 'smith') Mp <- lapply(seq_len(nC), function(x) Rz[[x]] %*% randMats$Xp[[x]])
  #-----------------------------------------------------------------------------
  # Loop through contrasts
  #-----------------------------------------------------------------------------
  # Change name to avoid using "get", which is slow
  setnames(DT, outcome, 'outcome')
  for (j in setdiff(seq_len(nC), skip)) {

    # T-contrasts
    #-------------------------------------------------------
    if (con.type == 't') {
      if (perm.method == 'smith') {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          M <- cbind(Mp[[j]][perms[i, ], ], randMats$Zp[[j]])
          MtM <- solve(crossprod(M))
          DT[, brainGraph_GLM_fit_t(M, outcome, MtM, eC[[j]]), by=eval(mykey)]
        }
      } else {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          DT[, brainGraph_GLM_fit_t(Mp[[j]], Rz[[j]][perms[i, ], ] %*% outcome, MtM[[j]], eC[[j]]), by=eval(mykey)]
        }
      }

    # F-contrasts
    #-------------------------------------------------------
    } else if (con.type == 'f') {
      if (perm.method == 'smith') {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          M <- cbind(Mp[[j]][perms[i, ], ], randMats$Zp[[j]])
          MtM <- solve(crossprod(M))
          CMtM <- solve(eC[[j]] %*% tcrossprod(MtM, eC[[j]]))
          DT[, brainGraph_GLM_fit_f(M, outcome, dfR, eC[[j]], rkC[j], CMtM), by=eval(mykey)]
        }
      } else {
        null.dist[[j]] <- foreach(i=seq_len(N), .combine='rbind') %dopar% {
          DT[, brainGraph_GLM_fit_f(Mp[[j]], Rz[[j]][perms[i, ], ] %*% outcome, dfR, eC[[j]], rkC[j], CMtM[[j]]), by=eval(mykey)]
        }
      }
    }
    null.dist[[j]] <- cbind(null.dist[[j]], data.table(perm=perm.order))
  }
  null.dist <- rbindlist(null.dist, idcol='contrast')
  return(null.dist)
}
