#' GLM non-parametric permutation testing
#'
#' \code{randomise} and \code{randomise_3d} perform non-parametric permutation
#' testing for analyses in which there is a single or multiple design matrix per
#' region, respectively. In the latter case, \code{X} should be a 3D array.
#'
#' @name randomise
#' @rdname randomise
NULL

#' Partition a design matrix into columns of interest and nuisance
#'
#' \code{partition} partitions a full design matrix into separate matrices of
#' covariates of interest and nuisance covariates based on a given contrast and
#' partition method.
#'
#' @section Model partitioning:
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
#' The manner of partitioning depends on the method. For example, for the
#' \code{guttman} method, \code{X} is formed from the columns of \code{contrast}
#' that have non-zero entries.
#'
#' @param M Numeric matrix or array of the full design matrix(es)
#' @param contrast For \code{partition}, a numeric matrix with 1 or more rows
#'   (for T and F contrasts, respectively) representing a \emph{single
#'   contrast}.
#' @inheritParams GLM
#' @importFrom MASS Null ginv
#' @importFrom abind abind adrop
#' @export
#' @rdname randomise
#'
#' @return \code{partition} returns a list containing:
#'   \item{Mp}{Numeric array; the combined partitioned arrays}
#'   \item{X}{Numeric array for the covariates of interest}
#'   \item{Z}{Numeric array for the nuisance covariates}
#'   \item{eCm}{The \emph{effective contrast}, equivalent to the original, for
#'     the partitioned model \code{[X, Z]} and considering all covariates}
#'   \item{eCx}{Same as \code{eCm}, but considering only \code{X}}
#' @references Beckmann, C.F. and Jenkinson, M. and Smith, S.M. (2001) General
#'   multi-level linear mdoelling for group analysis in FMRI. Tech Rep.
#'   University of Oxford, Oxford.
#' @references Guttman, I. (1982) \emph{Linear Models: An Introduction}. Wiley,
#'   New York.
#' @references Ridgway, G.R. (2009) Statistical analysis for longitudinal MR
#'   imaging of dementia. PhD thesis.

partition <- function(M, contrast, part.method=c('beckmann', 'guttman', 'ridgway')) {
  # Convert single matrices to array, for simplicity. Will drop extra dim at the end
  dimM <- dim(M)
  if (length(dimM) == 2L) {
    M <- array(M, dim=c(dimM, 1L), dimnames=c(dimnames(M), NULL))
    dimM <- c(dimM, 1L)
  }

  part.method <- match.arg(part.method)
  if (part.method == 'guttman') {
    idx <- unique(which(contrast != 0, arr.ind=TRUE)[, 2L])
    X <- M[, idx, , drop=FALSE]
    Z <- M[, -idx, , drop=FALSE]
    k <- length(idx)
    eCm <- cbind(contrast[, idx, drop=FALSE], contrast[, -idx, drop=FALSE])

  } else {
    subjs <- dimnames(M)[[1L]]
    regions <- dimnames(M)[[3L]]
    k <- qr.default(contrast)$rank
    rZ <- qr.default(M[, , 1L])$rank - k
    X <- array(0, dim=c(dimM[1L], k, dimM[3L]), dimnames=list(subjs, dimnames(contrast)[[1L]], regions))
    Z <- array(0, dim=c(dimM[1L], rZ, dimM[3L]), dimnames=list(subjs, NULL, regions))
    if (part.method == 'beckmann') {
      Cu <- MASS::Null(t(contrast))

      for (r in seq_len(dimM[3L])) {
        D <- inv(M[, , r])
        cdc <- inv(contrast %*% D, contrast, transpose=TRUE)
        Pc <- crossprod(contrast, cdc) %*% contrast %*% D
        Cv <- Cu - Pc %*% Cu
        F3 <- inv(Cv, D %*% Cv)
        X[, , r] <- M[, , r] %*% tcrossprod(D, contrast) %*% cdc
        Z[, , r] <- M[, , r] %*% D %*% Cv %*% F3
      }

    } else if (part.method == 'ridgway') {
      pinvC <- MASS::ginv(contrast)
      C0 <- diag(dimM[2L]) - crossprod(contrast, t(pinvC))

      for (r in seq_len(dimM[3L])) {
        tmpX <- M[, , r] %*% pinvC
        tmpZ <- svd(M[, , r] %*% C0)$u
        Z[, , r] <- tmpZ[, 1L:rZ, drop=FALSE]
        X[, , r] <- tmpX - Z[, , r] %*% MASS::ginv(Z[, , r]) %*% tmpX
      }
    }
    eCm <- cbind(diag(k), matrix(0, nrow=k, ncol=rZ))
  }
  eCx <- eCm[, 1L:k, drop=FALSE]
  if (dimM[3L] == 1L) {
    X <- abind::adrop(X, drop=3L)
    Z <- abind::adrop(Z, drop=3L)
  }
  Mp <- abind::abind(X, Z, along=2L)
  return(list(Mp=Mp, X=X, Z=Z, eCm=eCm, eCx=eCx))
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
#'   \item Calculate the inverse of the cross product of the full model
#'   \item Calculate the \dQuote{hat} and residual-forming matrices due to nuisance
#'     alone
#'   \item For F contrasts, return an inverse of the cross product between the
#'     contrast and the inverted design matrix, and the contrast's rank
#' }
#'
#' For F-contrasts, the \code{CMtM} matrix is
#' \deqn{(C (M^T M)^{-1} C^T)^{-1}}
#'
#' @inheritParams GLM
#' @inheritParams randomise
#' @return A list containing, for each of the following, a list with length
#'   equal to the number of contrasts:
#'   \item{Mp}{The full partitioned model, joined. For Still-White, this is
#'     actually just the \code{Xp} matrix.}
#'   \item{Rz}{The residual-forming matrix}
#'   \item{MtM}{The inverse of the cross product of the full model \code{Mp}}
#'   \item{eC}{The \emph{effective contrast}, equivalent to the original, for
#'     the partitioned model \code{[X, Z]} and considering all covariates}
#'   \item{CMtM}{The effective contrast \code{eC} multiplied by the cross
#'     product of the full model. For F-contrasts it is the inverse; for
#'     T-contrasts, the diagonal.}
#'   \item{rkC}{Rank of the effective contrast matrix.}
#'   \item{Xp}{(for Draper-Stoneman and Smith methods) List of matrices of
#'     covariates of interest}
#'   \item{Zp}{(for Draper-Stoneman and Smith methods) List of matrices of
#'     nuisance covariates}
#' @noRd

setup_randomise <- function(perm.method, part.method, X, contrasts, con.type, nC) {
  if (con.type == 't') contrasts <- matrix2list(contrasts)
  stopifnot(length(contrasts) == nC)
  Mp <- MtM <- Rz <- eC <- Xp <- Zp <- CMtM <- vector('list', length=nC)
  rkC <- rep.int(0L, nC)

  dimX <- dim(X)
  for (j in seq_len(nC)) {
    parts <- partition(X, contrasts[[j]], part.method)
    if (perm.method == 'stillWhite') {
      Mp[[j]] <- parts$X
      eC[[j]] <- parts$eCx
    } else {
      Mp[[j]] <- parts$Mp
      eC[[j]] <- parts$eCm
    }
    MtM[[j]] <- inv(Mp[[j]])
    Xp[[j]] <- parts$X
    Zp[[j]] <- parts$Z
    rkC[j] <- qr.default(eC[[j]])$rank

    Xpinv <- function(x) x %*% pinv(x)
    # Handle multiple design matrices
    if (length(dimX) == 3L) {
      In <- array(diag(dimX[1L]), dim=dimX[c(1L, 1L, 3L)])
      Hz <- switch(perm.method,
                   draperStoneman=, manly=array(0, dim=dim(In)),
                   terBraak=array(apply(Mp[[j]], 3L, Xpinv), dim=dim(In)),
                   array(apply(Zp[[j]], 3L, Xpinv), dim=dim(In)))  # All others
      CXtXfun <- if (part.method == 'guttman') cxtxfun_3d(con.type, transpose=FALSE)
        else cxtxfun_3d(con.type)
      CMtM[[j]] <- CXtXfun(eC[[j]], MtM[[j]], rkC[j], dimX[3L])
    } else {
      In <- diag(dimX[1L])
      # For Draper-Stoneman and Manly, Rz will be the identity matrix
      # NOTE: for ter Braak, this is really "Hm", and "Rz" is really "Rm"
      Hz <- switch(perm.method,
                   draperStoneman=, manly=matrix(0, dimX[1L], dimX[1L]),
                   freedmanLane=, smith=, stillWhite=Zp[[j]] %*% pinv(Zp[[j]]),
                   terBraak=Mp[[j]] %*% pinv(Mp[[j]]))

      CMtM[[j]] <- if (con.type == 'f') inv(eC[[j]] %*% MtM[[j]], eC[[j]], transpose=TRUE)
        else diag_sq(tcrossprod(eC[[j]] %*% MtM[[j]], eC[[j]]), rkC[j])
    }
    Rz[[j]] <- In - Hz
    if (perm.method %in% c('draperStoneman', 'smith')) {
      Xp[[j]] <- if (length(dimX) == 2L) Rz[[j]] %*% Xp[[j]]
        else array(vapply(seq_len(dimX[3L]), function(r) Rz[[j]][, , r] %*% Xp[[j]][, , r],
                          numeric(prod(dim(Xp[[j]])[1L:2L]))), dim=dim(Xp[[j]]))
      dimnames(Xp[[j]]) <- dimnames(parts$X)
    }
  }

  out <- list(Mp=Mp, Rz=Rz, MtM=MtM, eC=eC, CMtM=CMtM, rkC=rkC)
  if (perm.method %in% c('draperStoneman', 'smith', 'stillWhite')) out <- c(out, list(Xp=Xp, Zp=Zp))
  return(out)
}

#' @section Permutation methods:
#' The permutation methods can be split into 2 groups, depending on which part
#' of the model they permute. For full details, see \emph{Winkler et al., 2014}.
#' \describe{
#'   \item{Permute Y}{Freedman-Lane, Manly, and ter Braak}
#'   \item{Permute X}{Smith, Draper-Stoneman, and Still-White}
#' }
#'
#' Depending on the size of the data, it may be faster to use a method that
#' permutes \code{Y} instead of \code{X}. For example, in \code{\link{NBS}} with
#' dense matrices (more than 400-500 edges), it will be somewhat faster to use
#' the \dQuote{Smith} method compared to \dQuote{Freedman-Lane}. If using
#' \code{\link{brainGraph_GLM}}, the number of vertices follows the same
#' relationship.
#'
#' Furthermore, all methods except Still-White include the \code{Z} (nuisance
#' covariate) matrix when calculating the permuted statistics.
#'
#' @param X Numeric matrix or 3D array of the design matrix(es)
#' @param y Numeric matrix of outcome variables, with 1 column per region, or a
#'   single column if there is a different design matrix per region
#' @param ctype The contrast type
#' @param nC Integer; the number of contrasts
#' @param skip Integer vector indicating which (if any) contrasts to skip. Only
#'   used by \code{\link{NBS}}.
#' @param n,p,ny,dfR Integers for the number of observations, design matrix
#'   columns (its rank), number of regions/outcome variables, and residual
#'   degrees of freedom, respectively
#' @export
#' @inheritParams GLM
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#' @return A numeric array with dimensions \eqn{n_y \times N \times n_c};
#'   the number of rows equals number of regions/outcome variables, number of
#'   columns equals \code{N}, and the 3rd dimension is the number of contrasts
#' @rdname randomise
#' @references Draper, N.R. and Stoneman, D.M. (1966) Testing for the inclusion
#'   of variables in linear regression by a randomisation technique.
#'   \emph{Technometrics}. \bold{8(4)}, 695--699.
#' @references Freedman, D. and Lane, D. (1983) A nonstochastic interpretation
#'   of reported significance levels. \emph{J Bus Econ Stat}, \bold{1(4)},
#'   292--298. \url{https://dx.doi.org/10.1080/07350015.1983.10509354}
#' @references Manly B.F.J. (1986) Randomization and regression methods for
#'   testing for associations with geographical, environmental, and biological
#'   distances between populations. \emph{Res Popul Ecol}. \bold{28(2)},
#'   201--218.
#' @references Nichols, T.E. and Holmes, A.P. (2001) Nonparametric permutation
#'   tests for functional neuroimaging: A primer with examples. \emph{Human
#'   Brain Mapping}. \bold{15(1)}, 1--25.
#'   \url{https://dx.doi.org/10.1002/hbm.1058}
#' @references Smith, S.M. and Jenkinson, M. and Beckmann, C. and Miller, K. and
#'   Woolrich, M. (2007) Meaningful design and contrast estimability in fMRI.
#'   \emph{NeuroImage}. \bold{34(1)}, 127--36.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2006.09.019}
#' @references Still, A.W. and White, A.P. (1981) The approximate randomization
#'   test as an alternative to the F test in analysis of variance. \emph{Br J
#'   Math Stat Psychol}. \bold{34(2)}, 243--252.
#' @references ter Braak, C.J.F. 1992. Permutation versus bootstrap significance
#'   tests in multiple regression and ANOVA. \emph{Bootstrapping and related
#'   techniques}. Springer, Berlin, Heidelberg. 79--85.
#' @references Winkler, A.M. and Ridgway, G.R. and Webster, M.A. and Smith, S.M.
#'   and Nichols, T.E. (2014) Permutation inference for the general linear
#'   model. \emph{NeuroImage}. \bold{92}, 381--397.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2014.01.060}

randomise <- function(perm.method, part.method, N, perms, X, y, contrasts,
                      ctype, nC, skip=NULL, n=dim(X)[1L], p=qr.default(X)$rank,
                      ny=dim(y)[2L], dfR=n-p) {
  i <- NULL
  randMats <- setup_randomise(perm.method, part.method, X, contrasts, ctype, nC)
  Mp <- randMats$Mp; eC <- randMats$eC; CMtM <- randMats$CMtM; rkC <- randMats$rkC

  statfun <- switch(ctype, t=fastLmBG_t_rand, f=fastLmBG_f_rand)
  cxtxfun <- if (ctype == 't') function(con, x, k) diag_sq(tcrossprod(con %*% x, con), k)
    else function(con, x, k) inv(con %*% x, con, transpose=TRUE)
  combfun <- cbind
  doReshape <- FALSE
  null.dist <- array(0, dim=c(ny, N, nC))

  # If X is the same for all contrasts (Guttman), get all stats at one time
  #-----------------------------------------------------------------------------
  if (nC > 1L && all.identical(Mp)) {
    if (ctype == 't') {
      eC[[1L]] <- do.call(rbind, eC)
      rkC <- qr.default(eC[[1L]])$rank
      combfun <- rbind
      CMtM[[1L]] <- cxtxfun(eC[[1L]], randMats$MtM[[1L]], rkC[1L])
    }
    null.dist <- array(0, dim=c(ny * N, nC, 1L))
    nCorig <- nC; nC <- 1L; skip <- NULL  # Redefine so we only use the 1 combined contrast
    doReshape <- TRUE
  }
  #-----------------------------------------------------------------------------
  # Loop through contrasts
  #-----------------------------------------------------------------------------
  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  for (j in setdiff(seq_len(nC), skip)) {
    M <- Mp[[j]]; con <- eC[[j]]; rk <- rkC[j]

    if (perm.method %in% c('draperStoneman', 'smith', 'stillWhite')) {
      if (perm.method == 'stillWhite') {
        yz <- randMats$Rz[[j]] %*% y
        Mpfun <- function(M, porder, xcols) M[porder, , drop=FALSE]
        p <- rk
        df <- n - p
      } else {
        yz <- y; df <- dfR
        Mpfun <- function(M, porder, xcols) {
          M[, xcols] <- M[porder, xcols]
          M
        }
      }
      xcols <- seq_len(dim(randMats$Xp[[j]])[2L])
      I <- diag(1, n, p)  # Multiplier for "qr_Q2" (and "qr.qy")
      null.dist[, , j] <- foreach(i=seq_len(N), .combine=combfun) %dopar% {
        M <- Mpfun(M, perms[i, ], xcols)
        QR <- qr.default(M)
        Q <- qr_Q2(QR, I, n, p)
        R <- qr_R2(QR, p)
        MtM <- chol2inv(QR$qr, p)
        CMtM <- cxtxfun(con, MtM, rk)
        statfun(M, yz, con, QR, Q, R, n, p, ny, df, MtM, CMtM, rk)
      }

    # Freedman-Lane, Manly, & Ter Braak
    } else {
      Minv <- randMats$MtM[[j]]; CMinv <- CMtM[[j]]
      QR <- qr.default(M)
      Q <- qr_Q2(QR, n=n, p=p)
      R <- qr_R2(QR, p)
      yz <- randMats$Rz[[j]] %*% y
      if (ny == 1L) {
        yzMat <- matrix(0, n, N)
        for (i in seq_len(N)) yzMat[, i] <- yz[perms[i, ], ]
        null.dist[, , j] <- statfun(M, yzMat, con, QR, Q, R, n, p, N, dfR, Minv, CMinv, rk)
      } else {
        null.dist[, , j] <- foreach(i=seq_len(N), .combine=combfun) %dopar% {
          statfun(M, yz[perms[i, ], ], con, QR, Q, R, n, p, ny, dfR, Minv, CMinv, rk)
        }
      }
    }
  }
  if (isTRUE(doReshape)) null.dist <- array(null.dist, dim=c(ny, N, nCorig))
  return(null.dist)
}

#' @param runX Character vector of regions with non-singular designs
#' @export
#' @rdname randomise

randomise_3d <- function(perm.method, part.method, N, perms, X, y, contrasts,
                         ctype, nC, runX=dimnames(X)[[3L]], n=dim(X)[1L],
                         p=qr.default(X[, , 1L])$rank, ny=length(runX), dfR=n-p) {
  i <- NULL
  X <- X[, , runX, drop=FALSE]
  randMats <- setup_randomise(perm.method, part.method, X, contrasts, ctype, nC)
  Mp <- randMats$Mp; eC <- randMats$eC; CMtM <- randMats$CMtM; rkC <- randMats$rkC

  randfun <- switch(ctype, t=fastLmBG_t_rand_3d, f=fastLmBG_f_rand_3d)
  fitfun <- switch(perm.method, freedmanLane=, stillWhite=, terBraak=fastLmBG_3dY,
                   fastLmBG_3d)
  cxtxfun <- cxtxfun_3d(ctype)
  combfun <- cbind
  doReshape <- FALSE
  null.dist <- array(0, dim=c(ny, N, nC))

  # If X is the same for all contrasts (Guttman), get all stats at one time
  #-----------------------------------------------------------------------------
  if (nC > 1L && all.identical(Mp)) {
    if (ctype == 't') {
      eC[[1L]] <- do.call(rbind, eC); CMtM[[1L]] <- do.call(rbind, CMtM)
      rkC <- qr.default(eC[[1L]])$rank
      combfun <- rbind
      cxtxfun <- cxtxfun_3d(ctype, transpose=FALSE)
    }
    null.dist <- array(0, dim=c(ny * N, nC, 1L))
    nCorig <- nC; nC <- 1L  # Re-define so we don't loop past the 1 combined contrast
    doReshape <- TRUE
  }
  #-----------------------------------------------------------------------------
  # Loop through contrasts
  #-----------------------------------------------------------------------------
  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  for (j in seq_len(nC)) {
    M <- Mp[[j]]; con <- eC[[j]]; rk <- rkC[j]

    if (perm.method %in% c('draperStoneman', 'smith', 'stillWhite')) {
      if (perm.method == 'stillWhite') {
        yz <- apply(randMats$Rz[[j]], 3L,  `%*%`, y)
        dimnames(yz)[[2L]] <- runX
        Mpfun <- function(M, porder, xcols) M[porder, , , drop=FALSE]
        p <- rk
        fitfun <- if (p == 1L) fastLmBG_3dY_1p else fastLmBG_3dY
        dfR <- n - p
      } else {
        yz <- y
        Mpfun <- function(M, porder, xcols) {
          M[, xcols, ] <- M[porder, xcols, ]
          M
        }
      }
      xcols <- seq_len(dim(randMats$Xp[[j]])[2L])
      vnames <- dimnames(Mp[[j]])[[2L]]
      I <- diag(1, n, p)  # Multiplier for "qr_Q2" (and "qr.qy")
      null.dist[, , j] <- foreach(i=seq_len(N), .combine=combfun) %dopar% {
        M <- Mpfun(M, perms[i, ], xcols)
        QR <- qr(M)
        Q <- lapply(QR, qr_Q2, I, n, p)
        R <- lapply(QR, qr_R2, p)
        MtM <- inv(QR, p, ny, vnames, runX)
        CMtM <- cxtxfun(con, MtM, rk, ny)
        randfun(fitfun, M, yz, con, runX, QR, Q, R, n, p, ny, dfR, MtM, CMtM, rk)
      }

    # Freedman-Lane, Manly, & Ter Braak
    } else {
      Minv <- randMats$MtM[[j]]; CMinv <- CMtM[[j]]
      QR <- qr(M)
      Q <- lapply(QR, qr_Q2, n=n, p=p)
      R <- lapply(QR, qr_R2, p)
      yz <- if (perm.method == 'manly') y else apply(randMats$Rz[[j]], 3L, `%*%`, y)
      if (dim(yz)[2L] > 1L) dimnames(yz)[[2L]] <- runX

      null.dist[, , j] <- foreach(i=seq_len(N), .combine=combfun) %dopar% {
        randfun(fitfun, M, yz[perms[i, ], ], con, runX, QR, Q, R, n, p, ny, dfR, Minv, CMinv, rk)
      }
    }
  }
  if (isTRUE(doReshape)) null.dist <- array(null.dist, dim=c(ny, N, nCorig))
  return(null.dist)
}
