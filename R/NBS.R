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
#' @param A Three-dimensional array of all subjects' connectivity matrices
#' @param p.init Numeric; the initial p-value threshold (default: \code{0.001})
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM_design}}
#' @inheritParams GLM
#' @inheritParams symmetrize_mats
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
#'
#' @family Group analysis functions
#' @seealso \code{\link{brainGraph_GLM_design}, \link{brainGraph_GLM_fit_t}}
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
                p.init=0.001, perm.method=c('freedmanLane', 'terBraak', 'smith'),
                part.method=c('beckmann', 'guttman', 'ridgway'), N=1e3,
                perms=NULL, symm.by=c('max', 'min', 'avg'),
                alternative=c('two.sided', 'less', 'greater'), long=FALSE, ...) {
  i <- value <- Var1 <- Var2 <- Var3 <- p <- stat <- V1 <- contrast <- p.perm <- csize <- perm <- NULL
  dimA <- dim(A)
  stopifnot(dimA[3L] == dim(covars)[1L])
  Nv <- dimA[1L]

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- match.arg(alternative)
  if (ctype == 'f') alt <- 'two.sided'
  glmSetup <- setup_glm(covars, X, contrasts, ctype, con.name, measure=NULL, outcome=NULL, DT.y.m=NULL, level=NULL, ...)
  X <- glmSetup$X; incomp <- glmSetup$incomp; contrasts <- glmSetup$contrasts; nC <- glmSetup$nC
  if (length(incomp) > 0) A <- A[, , -glmSetup$covars[, which(get(getOption('bg.subject_id')) %in% incomp)]]

  # Get the outcome variables into a data.table; symmetrize and 0 the lower triangle for speed
  A <- symmetrize_array(A, symm.by)
  inds.low <- which(lower.tri(A[, , 1], diag=TRUE), arr.ind=TRUE)
  for (k in seq_len(dimA[3L])) A[cbind(inds.low, k)] <- 0
  A.m <- as.data.table(A)
  setnames(A.m, c('Var1', 'Var2', 'Var3', 'value'))
  pos.vals <- A.m[, sum(value) > 0, by=list(Var1, Var2)][V1 == 1, !'V1']
  A.m.sub <- A.m[pos.vals]

  # Do the model fitting/estimation and filter based on "p.init"
  DT.lm <- glm_fit_helper(A.m.sub, X, ctype, contrasts, alt, 'value', 'Var1,Var2')
  DT.lm <- DT.lm[p < p.init, list(contrast, Var1, Var2, stat, p)]

  # Create stat and p-val matrices, and get observed components
  comps.obs <- comps.perm <- vector('list', nC)
  T.max <- p.mat <- array(0, dim=c(Nv, Nv, nC))
  nrows <- DT.lm[, .N, by=contrast]$N
  skip <- which(nrows == 0L)
  if (length(skip) > 0L) comps.obs[skip] <- data.table(csize=0)
  for (j in setdiff(seq_len(nC), skip)) {
    T.max[, , j][DT.lm[contrast == j, cbind(Var1, Var2)]] <- DT.lm[contrast == j, stat]
    p.mat[, , j][DT.lm[contrast == j, cbind(Var1, Var2)]] <- DT.lm[contrast == j, p]

    clusts <- components(graph_from_adjacency_matrix(T.max[, , j], diag=FALSE, mode='undirected', weighted=TRUE))
    comps.obs[[j]] <- data.table(csize=with(clusts, sort(csize[csize > 1L], decreasing=TRUE)))
  }
  comps.obs <- rbindlist(comps.obs, idcol='contrast')

  part.method <- match.arg(part.method)
  perm.method <- match.arg(perm.method)
  out <- list(covars=glmSetup$covars, X=X, con.type=ctype, contrasts=contrasts, con.name=glmSetup$con.name,
              alt=alt, p.init=p.init, removed.subs=incomp, T.mat=T.max, p.mat=p.mat,
              N=N, perm.method=perm.method, part.method=part.method)
  if (length(skip) == nC) {
    out <- c(out, list(components=list(observed=comps.obs, permuted=NULL)))
    class(out) <- c('NBS', class(out))
    return(out)
  }

  # Create a null distribution of maximum component sizes
  #---------------------------------------------------------
  dimX <- dim(X)
  if (is.null(perms)) perms <- shuffleSet(n=dimX[1L], nset=N)
  null.dist <- randomise(perm.method, part.method, N, perms, contrasts, ctype, nC, skip,
                         A.m.sub, outcome='value', X, mykey='Var1,Var2')
  dfR <- dimX[1L] - dimX[2L]
  eqn <- if (ctype == 't') 'gamma / se' else 'numer / (se / dfR)'
  null.dist[, stat := eval(parse(text=eqn))]
  if (ctype == 't') {
    statfun <- switch(alt,
                      two.sided=function(stat, df) abs(stat) > qt(p.init / 2, df, lower.tail=FALSE),
                      less=function(stat, df) stat < qt(p.init, df),
                      greater=function(stat, df) stat > qt(p.init, df, lower.tail=FALSE))
    null.dist <- null.dist[statfun(stat, dfR), list(contrast, Var1, Var2, stat, perm)]
  } else {
    statfun <- function(stat, dfN, dfD) stat > qf(p.init / 2, dfN, dfD, lower.tail=FALSE)
    rkC <- vapply(contrasts, function(x) qr(x)$rank, integer(1))
    null.dist <- split(null.dist, by='contrast')
    for (j in setdiff(seq_len(nC), skip)) {
      null.dist[[j]] <- null.dist[[j]][statfun(stat, rkC[j], dfR), list(Var1, Var2, stat, perm)]
    }
    null.dist <- rbindlist(null.dist, idcol='contrast')
  }

  # Get the maximum component for each contrast & permutation
  for (j in setdiff(seq_len(nC), skip)) {
    comps.perm[[j]] <- foreach(i=null.dist[contrast == j, unique(perm)], .combine='c') %dopar% {
      T.mat.tmp <- matrix(0, Nv, Nv)
      T.mat.tmp[null.dist[contrast == j & perm == i, cbind(Var1, Var2)]] <- null.dist[contrast == j & perm == i, stat]
      max(components(graph_from_adjacency_matrix(T.mat.tmp, diag=F, mode='undirected', weighted=TRUE))$csize)
    }
    if (length(comps.perm[[j]]) < N) comps.perm[[j]] <- c(comps.perm[[j]], rep(0, N - length(comps.perm[[j]])))
    comps.perm[[j]] <- data.table(perm=comps.perm[[j]])
  }
  comps.perm <- rbindlist(comps.perm, idcol='contrast')
  kNumComps <- comps.obs[, .N, by=contrast]$N
  for (j in seq_len(nC)) {
    comps.obs[contrast == j,
              p.perm := (sum(comps.perm[contrast == j, perm] >= csize) + 1) / (N + 1),
              by=csize]
  }

  comps.out <- list(observed=comps.obs)
  if (isTRUE(long)) comps.out$permuted <- comps.perm
  out <- c(out, list(components=comps.out))
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
                                           numeric(1))]
  }

  nbs.sum <- c(object, list(printCon=contrast, DT.sum=nbs.dt, digits=digits))
  class(nbs.sum) <- c('summary.NBS', class(nbs.sum))
  return(nbs.sum)
}

#' @aliases summary.NBS
#' @method print summary.NBS
#' @export

print.summary.NBS <- function(x, ...) {
  print_title_summary('Network-based statistic results')
  print_permutation_summary(x)

  cat('Initial p-value: ', x$p.init, '\n\n')

  print_contrast_type_summary(x)
  print_subs_summary(x)

  # Print results for each contrast
  message('\n', 'Statistics', '\n', rep('-', getOption('width') / 4))
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
  terms(structure(x, class='bg_GLM'))
}

#' @export
#' @rdname NBS
formula.NBS <- function(x, ...) {
  x$outcome <- x$measure <- 'weight'
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
