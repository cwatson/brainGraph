#' Network-based statistic for brain MRI data
#'
#' Calculates the \emph{network-based statistic (NBS)}, which allows for
#' family-wise error (FWE) control over network data, introduced for brain MRI
#' data by Zalesky et al. Requires a three-dimensional array of all subjects'
#' connectivity matrices and a \code{data.table} of covariates, in addition to a
#' contrast matrix or list. A null distribution of the largest connected
#' component size is created by fitting a GLM to permuted data. For details, see
#' \code{\link{GLM}}.
#'
#' When printing a \code{summary}, you can include arguments to
#' \code{\link[stats]{printCoefmat}}.
#'
#' @note It is assumed that the order of the subjects in \code{covars} matches
#' that of the input array \code{A}. You will need to ensure that this is the
#' case. Prior to \code{v3.0.0}, the \code{covars} table was sorted by
#' \code{Study.ID} before creating the design matrix.
#'
#' @param A Three-dimensional array of all subjects' connectivity matrices
#' @param p.init Numeric; the initial p-value threshold (default: \code{0.001})
#' @param symm.by Character string; how to create symmetric off-diagonal
#'   elements. Default: \code{max}
#' @inheritParams GLM
#' @export
#' @importFrom permute shuffleSet
#'
#' @return An object of class \code{NBS} with some input arguments in addition
#'   to:
#'   \item{X}{The design matrix}
#'   \item{removed.subs}{Character vector of subject ID's removed due to
#'     incomplete data (if any)}
#'   \item{T.mat}{3-d array of (symmetric) numeric matrices containing the
#'     statistics for each edge}
#'   \item{p.mat}{3-d array of (symmetric) numeric matrices containing the
#'     P-values}
#'   \item{components}{List containing data tables of the observed and permuted
#'     connected component sizes and P-values}
#'   \item{rank,df.residual,qr,cov.unscaled}{The rank, residual degrees of
#'     freedom, QR decomposition, and unscaled covariance matrix of the design
#'     matrix}
#'
#' @family Group analysis functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Zalesky, A. and Fornito,  A. and Bullmore, E.T. (2010)
#'   Network-based statistic: identifying differences in brain networks.
#'   \emph{NeuroImage}, \bold{53(4)}, 1197--1207.
#'   \url{https://dx.doi.org/10.1016/j.neuroimage.2010.06.041}
#' @examples
#' \dontrun{
#' max.comp.nbs <- NBS(A.norm.sub[[1]], covars.dti, N=5e3)
#' }

NBS <- function(A, covars, contrasts, con.type=c('t', 'f'), X=NULL, con.name=NULL,
                p.init=0.001, perm.method=c('freedmanLane', 'terBraak', 'smith',
                                            'draperStoneman', 'manly', 'stillWhite'),
                part.method=c('beckmann', 'guttman', 'ridgway'), N=1e3,
                perms=NULL, symm.by=c('max', 'min', 'avg'),
                alternative=c('two.sided', 'less', 'greater'), long=FALSE, ...) {
  contrast <- p.perm <- csize <- NULL
  dimA <- dim(A)
  stopifnot(dimA[3L] == dim(covars)[1L])

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- if (ctype == 'f') 'two.sided' else match.arg(alternative)

  # Remove subjects w/ incomplete data, and create design matrix
  sID <- getOption('bg.subject_id')
  covars <- droplevels(copy(covars))
  if (!hasName(covars, sID)) covars[, eval(sID) := seq_len(dim(covars)[1L])]
  covars[, eval(sID) := check_sID(get(sID))]
  incomp <- covars[!complete.cases(covars), which=TRUE]
  names(incomp) <- covars[incomp, get(sID)]
  if (length(incomp) > 0L) {
    A <- A[, , -incomp, drop=FALSE]
    dimA <- dim(A)
    covars <- covars[-incomp]
  }
  if (is.null(X)) X <- brainGraph_GLM_design(covars, ...)
  tmp <- contrast_names(contrasts, ctype, con.name, X)
  contrasts <- tmp$contrasts; nC <- tmp$nC

  dimX <- dim(X)
  n <- dimX[1L]; p <- dimX[2L]; dfR <- n - p; nV <- dimA[1L]

  # Symmetrize the array and use only the upper triangle's data
  A <- symmetrize(A, symm.by)
  inds.high <- which(upper.tri(A[, , 1L]), arr.ind=TRUE)
  ny <- dim(inds.high)[1L]
  yMat <- matrix(0, n, ny)
  for (r in seq_len(ny)) yMat[, r] <- A[inds.high[r, 1L], inds.high[r, 2L], ]

  # Remove region pairs where the value is 0 for all subjects
  sums <- .colSums(yMat, n, ny)
  if (any(sums == 0)) {
    v <- which(sums == 0)
    inds.high <- inds.high[-v, , drop=FALSE]
    yMat <- yMat[, -v, drop=FALSE]
    ny <- dim(yMat)[2L]
  }

  # Do the model fitting/estimation and filter based on "p.init"
  fits <- fastLmBG(X, yMat)
  obs_stats <- if (ctype == 't') fastLmBG_t(fits, contrasts, alt) else fastLmBG_f(fits, contrasts)
  sig_inds <- which(obs_stats[, 'p', , drop=FALSE] < p.init, arr.ind=TRUE)[, -2L, drop=FALSE]

  # Create stat and p-val matrices, and get observed components
  comps.obs <- vector('list', nC)
  T.max <- p.mat <- array(0, dim=c(nV, nV, nC))
  nrows <- rep.int(0L, nC)
  if (length(sig_inds) > 0L) {
    nrows <- vapply(seq_len(nC), function(x) sum(sig_inds[, 2L] == x), integer(1L))
  }
  skip <- which(nrows == 0L)
  noSkip <- setdiff(seq_len(nC), skip)
  if (length(skip) > 0L) for (k in skip) comps.obs[[k]] <- data.table(csize=0)
  for (j in noSkip) {
    rowInds <- sig_inds[sig_inds[, 2L] == j, 1L]
    T.max[, , j][inds.high[rowInds, , drop=FALSE]] <- obs_stats[rowInds, 'stat', j]
    p.mat[, , j][inds.high[rowInds, , drop=FALSE]] <- obs_stats[rowInds, 'p', j]

    clusts <- components(graph_from_adjacency_matrix(T.max[, , j], diag=FALSE,
                                                     mode='undirected', weighted=TRUE))
    comps.obs[[j]] <- data.table(csize=with(clusts, sort(csize[csize > 1L], decreasing=TRUE)))
  }
  comps.obs <- rbindlist(comps.obs, idcol='contrast')

  out <- list(covars=covars, X=X, con.type=ctype, contrasts=contrasts, con.name=tmp$con.name,
              alt=alt, p.init=p.init, removed.subs=incomp, T.mat=T.max, p.mat=p.mat, N=N)
  out <- c(out, fits[c('rank', 'df.residual', 'qr', 'cov.unscaled')])
  if (length(skip) == nC) {
    out <- c(out, list(components=list(observed=comps.obs, permuted=NULL)))
    class(out) <- c('NBS', class(out))
    return(out)
  }

  # Create a null distribution of maximum component sizes
  #---------------------------------------------------------
  part.method <- match.arg(part.method); perm.method <- match.arg(perm.method)
  if (is.null(perms) || dim(perms)[2L] != n) perms <- shuffleSet(n=n, nset=N)
  null.dist <- randomise(perm.method, part.method, N, perms, X, yMat,
                         contrasts, ctype, nC, skip=skip, n=n, p=p, ny=ny, dfR=dfR)

  if (ctype == 't') {
    thresh <- switch(alt, two.sided=qt(p.init / 2, dfR, lower.tail=FALSE),
                     less=qt(p.init, dfR), greater=qt(p.init, dfR, lower.tail=FALSE))
    statfun <- switch(alt, two.sided=function(stat) abs(stat) > thresh, less=function(stat) stat < thresh,
                      greater=function(stat) stat > thresh)
    sig_inds_null <- which(statfun(null.dist), arr.ind=TRUE)
  } else {
    rkC <- vapply(contrasts, function(x) qr.default(x)$rank, integer(1L))
    thresh <- qf(p.init / 2, rkC, dfR, lower.tail=FALSE)
    sig_inds_null <- vector('list', nC)
    for (j in noSkip) {
      sig_inds_null[[j]] <- which(null.dist[, , j] > thresh[j], arr.ind=TRUE)
      sig_inds_null[[j]] <- cbind(sig_inds_null[[j]], dim3=j)
    }
    sig_inds_null <- do.call(rbind, sig_inds_null)
  }

  # Get the maximum component for each contrast & permutation
  comps.perm <- matrix(0, N, nC)
  compsJ <- rep.int(0, N)
  comps.perm[, noSkip] <- foreach(j=noSkip, .combine=cbind) %dopar% {
    sigJ <- sig_inds_null[sig_inds_null[, 3L] == j, 1L:2L]
    kNumEdges <- tabulate(sigJ[, 2L])

    # If there are 1-2 edges, the component size can only be 2 or 3
    singleEdge <- which(kNumEdges == 1L)
    doubleEdge <- which(kNumEdges == 2L)
    sigJ2 <- sigJ[sigJ[, 2L] %in% doubleEdge, ]
    compsDouble <- vapply(doubleEdge, function(x)
                          sum(tabulate(inds.high[sigJ2[sigJ2[, 2L] == x, 1L], ]) != 0),
                          numeric(1L))
    compsDouble[compsDouble == 4] <- 2
    compsJ[singleEdge] <- 2
    compsJ[doubleEdge] <- compsDouble

    doPerms <- which(kNumEdges > 2L)
    sigJ <- sigJ[sigJ[, 2L] %in% doPerms, ]
    for (i in doPerms) {
      rowInds <- sigJ[sigJ[, 2L] == i, 1L]
      compsJ[i] <- max(components(graph(t(inds.high[rowInds, , drop=FALSE]), directed=FALSE))$csize)
    }
    compsJ
  }

  p.perm <- unlist(lapply(noSkip, function(j)
                          vapply(comps.obs[contrast == j, csize], function(x)
                                 sum(comps.perm[, j] >= x) + 1L, numeric(1L))))
  p.perm <- p.perm / (N + 1L)
  comps.obs[contrast %in% noSkip, p.perm := p.perm]

  comps.out <- list(observed=comps.obs)
  if (isTRUE(long)) comps.out$permuted <- comps.perm
  out <- c(out, list(part.method=part.method, perm.method=perm.method, components=comps.out))
  class(out) <- c('NBS', class(out))
  return(out)
}

################################################################################
# METHODS
################################################################################

#' Print a summary of NBS analysis
#'
#' @param object,x A \code{NBS} object
#' @inheritParams summary.bg_GLM
#' @export
#' @rdname NBS

summary.NBS <- function(object, contrast=NULL, digits=max(3L, getOption('digits') - 2L), ...) {

  # Observed component sizes
  #--------------------------------------
  nbs.dt <- with(object, data.table(alt=alt, N=N, components$observed, ecount=0))
  g.nbs <- make_brainGraphList(object, guess_atlas(object$T.mat), set.attrs=FALSE, .progress=FALSE)
  for (j in seq_along(g.nbs[])) {
    nbs.dt[contrast == j, ecount := vapply(seq_len(.N), function(x)
                                           ecount(subset_graph(g.nbs[j], paste('comp ==', x))$g),
                                           numeric(1L))]
  }

  nbs.sum <- c(object, list(printCon=contrast, DT.sum=nbs.dt, digits=digits), list(...))
  class(nbs.sum) <- c('summary.NBS', class(object))
  return(nbs.sum)
}

#' @aliases summary.NBS
#' @method print summary.NBS
#' @export

print.summary.NBS <- function(x, ...) {
  print_title_summary('Network-based statistic results')
  print_model_summary(x)
  print_permutation_summary(x)

  cat('Initial p-value: ', x$p.init, '\n\n')

  print_contrast_type_summary(x)
  print_subs_summary(x)

  # Print results for each contrast
  message('\n', 'Statistics', '\n', rep.int('-', getOption('width') / 4))
  print_contrast_stats_summary(x, ...)

  invisible(x)
}

#' @export
#' @rdname NBS
#' @include glm_methods.R
nobs.NBS <- nobs.bg_GLM

#' @export
#' @rdname NBS
terms.NBS <- function(x, ...) {
  x$outcome <- x$measure <- 'weight'
  x$DT.Xy <- x$covars
  terms(structure(x, class='bg_GLM'))
}

#' @export
#' @rdname NBS
formula.NBS <- function(x, ...) {
  x$outcome <- x$measure <- 'weight'
  x$DT.Xy <- x$covars
  formula(structure(x, class='bg_GLM'))
}

#' @export
#' @rdname NBS
#' @include glm_methods.R
labels.NBS <- labels.bg_GLM

#' @method case.names NBS
#' @export
#' @rdname NBS
#' @include glm_methods.R
case.names.NBS <- case.names.bg_GLM

#' @method variable.names NBS
#' @export
#' @rdname NBS
#' @include glm_methods.R
variable.names.NBS <- variable.names.bg_GLM

#' @method df.residual NBS
#' @export
#' @rdname NBS
#' @include glm_stats.R
df.residual.NBS <- df.residual.bg_GLM
