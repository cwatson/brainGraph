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
  skip <- value <- Var1 <- Var2 <- Var3 <- p <- stat <- V1 <- contrast <- p.perm <- csize <- perm <- Study.ID <- NULL
  stopifnot(dim(A)[3] == nrow(covars))
  Nv <- nrow(A)

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- match.arg(alternative)
  if (ctype == 'f') alt <- 'two.sided'
  glmSetup <- setup_glm(covars, X, contrasts, ctype, con.name, measure=NULL, outcome=NULL, DT.y.m=NULL, level=NULL, ...)
  X <- glmSetup$X; incomp <- glmSetup$incomp; contrasts <- glmSetup$contrasts; nC <- glmSetup$nC

  # Get the outcome variables into a data.table; symmetrize and 0 the lower triangle for speed
  if (length(incomp) > 0) A <- A[, , -glmSetup$covars[, which(Study.ID %in% incomp)]]
  A <- symmetrize_array(A, symm.by)
  for (k in seq_len(dim(A)[3])) {
    x <- A[, , k]
    x[lower.tri(x)] <- 0
    A[, , k] <- x
  }
  A.m <- setDT(melt(A))
  inds.upper <- as.data.table(which(upper.tri(A[, , 1]), arr.ind=TRUE))
  setnames(inds.upper, c('Var1', 'Var2'))
  setkey(inds.upper, Var1, Var2)
  setkey(A.m, Var1, Var2)
  A.m <- A.m[inds.upper]
  setkey(A.m, Var1, Var2, Var3)
  pos.vals <- A.m[, sum(value) > 0, by=list(Var1, Var2)][V1 == 1, !'V1']
  A.m.sub <- A.m[pos.vals]

  # Do the model fitting/estimation
  DT.lm <- glm_fit_helper(A.m.sub, X, ctype, contrasts, alt, outcome='value', mykey='Var1,Var2')

  # Filter based on "p.init", and create stat and p-val matrices
  DT.lm <- DT.lm[p < p.init, list(Var1, Var2, stat, p, contrast)]
  comps.obs <- vector('list', nC)
  T.max <- p.mat <- array(0, dim=c(Nv, Nv, nC))
  compfun <- switch(alt,
                    two.sided=function(x) abs(x) > t(abs(x)),
                    less=function(x) x < t(x),
                    greater=function(x) x > t(x))
  for (j in seq_len(nC)) {
    T.mat <- matrix(0, Nv, Nv)
    if (nrow(DT.lm[contrast == j]) == 0) {
      warning(sprintf('No significant differences observed for contrast %i!', j))
      comps.obs[[j]] <- data.table(csize=0)
      skip <- c(skip, j)
      next
    }
    T.mat[DT.lm[contrast == j, cbind(Var1, Var2)]] <- DT.lm[contrast == j, stat]
    p.mat[, , j][DT.lm[contrast == j, cbind(Var1, Var2)]] <- DT.lm[contrast == j, p]

    inds.tr <- which(compfun(T.mat), arr.ind=TRUE)
    T.max[, , j] <- ifelse(compfun(T.mat), T.mat, t(T.mat))
    for (i in seq_len(nrow(inds.tr))) {
      p.mat[inds.tr[i, 2], inds.tr[i, 1], j] <- p.mat[inds.tr[i, 1], inds.tr[i, 2], j]
    }

    clusts <- components(graph_from_adjacency_matrix(T.max[, , j], diag=F, mode='undirected', weighted=TRUE))
    comps.obs[[j]] <- data.table(csize=sort(unique(clusts$csize), decreasing=TRUE))
  }
  comps.obs <- rbindlist(comps.obs, idcol='contrast')

  part.method <- match.arg(part.method)
  perm.method <- match.arg(perm.method)
  out <- list(covars=glmSetup$covars, X=X, p.init=p.init, con.type=ctype, contrasts=contrasts, con.name=glmSetup$con.name,
              alt=alt, N=N, removed.subs=incomp, T.mat=T.max, p.mat=p.mat,
              perm.method=perm.method, part.method=part.method)
  if (length(skip) == nC) {
    out <- c(out, list(components=list(observed=comps.obs, permuted=NULL)))
    class(out) <- c('NBS', class(out))
    return(out)
  }

  # Create a null distribution of maximum component sizes
  #---------------------------------------------------------
  if (is.null(perms)) perms <- shuffleSet(n=nrow(X), nset=N)

  comps.perm <- randomise_nbs(perm.method, part.method, ctype, N, perms, A.m.sub, nC, skip, p.init, X, contrasts, alt, Nv)
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
