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
  value <- Var1 <- Var2 <- Var3 <- p <- stat <- V1 <- contrast <- p.perm <- csize <- perm <- Study.ID <- NULL
  stopifnot(dim(A)[3] == nrow(covars))
  Nv <- nrow(A)

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- match.arg(alternative)
  if (ctype == 'f') alt <- 'two.sided'
  glmSetup <- setup_glm(covars, X, contrasts, ctype, con.name, measure=NULL, outcome=NULL, DT.y.m=NULL, level=NULL, ...)
  X <- glmSetup$X; incomp <- glmSetup$incomp; contrasts <- glmSetup$contrasts; nC <- glmSetup$nC
  if (length(incomp) > 0) A <- A[, , -glmSetup$covars[, which(Study.ID %in% incomp)]]

  # Get the outcome variables into a data.table; symmetrize and 0 the lower triangle for speed
  A <- symmetrize_array(A, symm.by)
  inds.low <- which(lower.tri(A[, , 1], diag=TRUE), arr.ind=TRUE)
  for (k in seq_len(dim(A)[3])) A[cbind(inds.low, k)] <- 0
  A.m <- setDT(melt(A))
  setkey(A.m, Var1, Var2, Var3)
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
    comps.obs[[j]] <- data.table(csize=sort(unique(clusts$csize), decreasing=TRUE))
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
  if (is.null(perms)) perms <- shuffleSet(n=nrow(X), nset=N)
  null.dist <- randomise(perm.method, part.method, N, perms, contrasts, ctype, nC, skip,
                         A.m.sub, outcome='value', X, mykey='Var1,Var2')
  dfR <- nrow(X) - ncol(X)
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
    comps.obs[contrast == j, p.perm := mapply(function(x, y) (sum(y >= x) + 1) / (N + 1),
                                              csize,
                                              rep(list(comps.perm[contrast == j, perm]), kNumComps[j]))]
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
#' @param object A \code{NBS} object
#' @inheritParams summary.bg_GLM
#' @export
#' @method summary NBS
#' @rdname NBS

summary.NBS <- function(object, contrast=NULL, digits=max(3L, getOption('digits') - 2L), ...) {
  contrast <- csize <- NULL

  # Observed component sizes
  #--------------------------------------
  ecounts <- vector('list', length(object$con.name))
  for (j in seq_along(object$con.name)) {
    if (sum(object$T.mat[, , j]) == 0) next  # No edges met initial criteria
    g.nbs <- graph_from_adjacency_matrix(object$T.mat[, , j], diag=F, mode='undirected', weighted=TRUE)
    clusts <- components(g.nbs)
    comps <- sort(unique(clusts$csize), decreasing=TRUE)
    z <- clusts$membership
    zsort <- match(z, order(table(z), decreasing=TRUE))
    z.ind <- unique(table(zsort))
    ecounts[[j]] <- rep(0, sum(z.ind > 1))
    for (i in seq_along(ecounts[[j]])) {
      ecounts[[j]][i] <- ecount(induced_subgraph(g.nbs, which(zsort == i)))
    }
  }

  nbs.dt <- with(object,
                 data.table(alt=alt, N=N, components$observed))
  nbs.dt[, ecount := 0]
  for (j in seq_along(object$con.name)) {
    if (sum(object$T.mat[, , j]) == 0) next  # No edges met initial criteria
    nbs.dt[contrast == j & csize > 1, ecount := ecounts[[j]]]
  }
  nbs.sum <- c(object, list(contrast=contrast, DT.sum=nbs.dt, digits=digits))
  class(nbs.sum) <- c('summary.NBS', class(nbs.sum))
  return(nbs.sum)
}

#' @aliases summary.NBS
#' @method print summary.NBS

print.summary.NBS <- function(x, ...) {
  `p-value` <- `# edges` <- csize <- p.perm <- NULL
  print_title_summary('Network-based statistic results')
  print_permutation_summary(x)

  cat('Initial p-value: ', x$p.init, '\n\n')

  print_contrast_type_summary(x)
  print_subs_summary(x)

  xdt <- x$DT.sum[csize > 1]
  setnames(xdt, c('csize', 'ecount'), c('# vertices', '# edges'))
  xdt[, `p-value` := signif(p.perm)]
  xdt[, c('alt', 'N', 'p.perm') := NULL]

  # Print results for each contrast
  message('\n', 'Statistics', '\n', rep('-', getOption('width') / 4))
  if (is.null(x$contrast)) {
    contrast <- xdt[, unique(contrast)]
  } else {
    contrast <- x$contrast
  }

  for (i in contrast) {
    message(x$con.name[i])
    if (nrow(xdt[contrast == i]) == 0) {
      message('\tNo signficant results!\n')
    } else {
      printCoefmat(xdt[contrast == i, !'contrast'], tst.ind=2, P.values=TRUE, has.Pvalue=TRUE, digits=x$digits, ...)
      cat('\n')
    }
  }
  invisible(x)
}
