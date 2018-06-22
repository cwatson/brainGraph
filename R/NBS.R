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
#' @param p.init Numeric; the initial p-value threshold (default: \code{0.001})
#' @param ... Other arguments passed to \code{\link{brainGraph_GLM_design}}
#' @inheritParams brainGraph_GLM
#' @inheritParams symmetrize_mats
#' @export
#' @importFrom permute shuffleSet
#'
#' @return An object of class \code{NBS} with some input arguments in addition
#'   to:
#'   \item{X}{The design matrix}
#'   \item{removed}{Character vector of subject ID's removed due to incomplete
#'     data (if any)}
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
#' @references Zalesky A., Fornito A., Bullmore E.T. (2010) \emph{Network-based
#'   statistic: identifying differences in brain networks}. NeuroImage,
#'   53(4):1197-1207.
#' @examples
#' \dontrun{
#' max.comp.nbs <- NBS(A.norm.sub[[1]], covars.dti, N=5e3)
#' }

NBS <- function(A, covars, con.mat, con.type=c('t', 'f'), X=NULL, con.name=NULL,
                p.init=0.001, N=1e3, perms=NULL, symm.by=c('max', 'min', 'avg'),
                alternative=c('two.sided', 'less', 'greater'), long=FALSE, ...) {
  skip <- value <- Var1 <- Var2 <- Var3 <- p <- stat <- V1 <- contrast <- p.perm <- csize <- perm <- Study.ID <- NULL
  stopifnot(dim(A)[3] == nrow(covars))
  Nv <- nrow(A)

  # Initial GLM setup
  ctype <- match.arg(con.type)
  alt <- match.arg(alternative)
  if (ctype == 'f') alt <- 'two.sided'
  glmSetup <- setup_glm(covars, X, con.mat, ctype, con.name, ...)
  covars <- glmSetup$covars; X <- glmSetup$X; incomp <- glmSetup$incomp
  con.mat <- glmSetup$con.mat; con.name <- glmSetup$con.name; nC <- glmSetup$nC

  # Get the outcome variables into a data.table; symmetrize and 0 the lower triangle for speed
  if (length(incomp) > 0) A <- A[, , -covars[, which(Study.ID %in% incomp)]]
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
  DT.lm <- glm_fit_helper(A.m.sub, X, ctype, con.mat, alt, 'value', 'Var1,Var2')$DT.lm

  # Filter based on "p.init", and create stat and p-val matrices
  DT.lm <- DT.lm[p < p.init, list(Var1, Var2, stat, p, contrast)]
  comps.obs <- vector('list', length=length(con.name))
  T.max <- p.mat <- array(0, dim=c(Nv, Nv, nC))
  compfun <- switch(alt,
                    two.sided=function(x) {abs(x) > t(abs(x))},
                    less=function(x) {x < t(x)},
                    greater=function(x) {x > t(x)})
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

  out <- list(X=X, p.init=p.init, con.type=ctype, con.mat=con.mat, con.name=con.name,
              alt=alt, N=N, removed=incomp, T.mat=T.max, p.mat=p.mat)
  if (length(skip) == length(con.name)) {
    out <- c(out, list(components=list(observed=comps.obs, permuted=NULL)))
    class(out) <- c('NBS', class(out))
    return(out)
  }

  # Create a null distribution of maximum component sizes
  #---------------------------------------------------------
  if (is.null(perms)) perms <- shuffleSet(n=nrow(X), nset=N)

  randMats <- setup_randomise(X, con.mat, nC)
  comps.perm <- randomise_nbs(ctype, N, perms, A.m.sub, nC, skip, randMats, p.init, alt, Nv)

  kNumComps <- comps.obs[, .N, by=contrast]$N
  for (j in seq_along(con.name)) {
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
# HELPER FUNCTION
################################################################################
randomise_nbs <- function(ctype, N, perms, DT, nC, skip, randMats, p.init, alt, Nv) {
  se <- perm <- Var1 <- Var2 <- i <- value <- stat <- numer <- NULL
  Mp <- randMats$Mp; Rz <- randMats$Rz; MtM <- randMats$MtM; eC <- randMats$eC
  dfR <- nrow(Mp[[1]]) - ncol(Mp[[1]])
  if (ctype == 't') {
    statfun <- switch(alt,
                      two.sided=function(stat, df) {abs(stat) > qt(p.init / 2, df, lower.tail=FALSE)},
                      less=function(stat, df) {stat <- qt(p.init, df)},
                      greater=function(stat, df) {stat > qt(p.init, df, lower.tail=FALSE)})
  } else {
    statfun <- function(stat, dfN, dfD) stat > qf(p.init / 2, dfN, dfD, lower.tail=FALSE)
    CMtM <- solve(eC[[1]] %*% MtM[[1]] %*% t(eC[[1]]))
    rkC <- qr(eC[[1]])$rank
  }
  maxfun.mat <- switch(alt,
                       two.sided=function(mat) {ifelse(abs(mat) > t(abs(mat)), mat, t(mat))},
                       less=function(mat) {ifelse(mat < t(mat), mat, t(mat))},
                       greater=function(mat) {ifelse(mat > t(mat), mat, t(mat))})
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
    if (sum(object$T.mat[[j]]) == 0) next  # No edges met initial criteria
    g.nbs <- graph_from_adjacency_matrix(object$T.mat[[j]], diag=F, mode='undirected', weighted=TRUE)
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
    if (sum(object$T.mat[[j]]) == 0) next  # No edges met initial criteria
    nbs.dt[contrast == j & csize > 1, ecount := ecounts[[j]]]
  }
  nbs.sum <- list(contrast=contrast, res.nbs=object, DT.sum=nbs.dt, digits=digits)
  class(nbs.sum) <- c('summary.NBS', class(nbs.sum))
  return(nbs.sum)
}

#' @aliases summary.NBS
#' @method print summary.NBS

print.summary.NBS <- function(x, ...) {
  `p-value` <- `# edges` <- csize <- p.perm <- NULL
  title <- 'Network-based statistic results'
  message('\n', title, '\n', rep('-', getOption('width') / 2))
  cat('Number of permutations: ', prettyNum(x$res.nbs$N, ','), '\n')
  cat('Initial p-value: ', x$res.nbs$p.init, '\n\n')

  cat('Contrast type: ', paste(toupper(x$res.nbs$con.type), 'contrast'), '\n')
  alt <- switch(x$res.nbs$alt,
                two.sided='C != 0',
                greater='C > 0',
                less='C < 0')
  cat('Alternative hypothesis: ', alt, '\n')
  cat('Contrast matrix: ', '\n')
  print(x$res.nbs$con.mat)

  if (length(x$res.nbs$removed) != 0) cat('\nSubjects removed due to incomplete data:\n', x$res.nbs$removed, '\n')

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
    message(x$res.nbs$con.name[i])
    if (nrow(xdt[contrast == i]) == 0) {
      message('\tNo signficant results!\n')
    } else {
      printCoefmat(xdt[contrast == i, !'contrast'], tst.ind=2, P.values=TRUE, has.Pvalue=TRUE, digits=x$digits, ...)
      cat('\n')
    }
  }
  invisible(x)
}
