#' Simulate N random graphs w/ same clustering and degree sequence as the input
#'
#' \code{sim.rand.graph.par} simulates \code{N} simple random graphs with the
#' same clustering (optional) and degree sequence as the input. Essentially a
#' wrapper for \code{\link[igraph]{sample_degseq}} (or, if you want to match by
#' clustering, \code{\link{sim.rand.graph.clust}}) and
#' \code{\link{make_brainGraph}}. It uses \code{\link[foreach]{foreach}} for
#' parallel processing.
#'
#' If you do not want to match by clustering, then simple rewiring of the input
#' graph is performed (the number of rewires equaling the larger of \code{1e4}
#' and \eqn{10 \times m}, where \eqn{m} is the graph's edge count).
#'
#' @param g A graph object
#' @param N Integer; the number of random graphs to simulate. Default: 100
#' @param clustering Logical; whether or not to control for clustering. Default:
#'   \code{FALSE}
#' @param rewire.iters Integer; number of rewiring iterations for the initial
#'   graph randomization. Default: 1e4
#' @param cl The clustering measure. Default: \emph{transitivity}
#' @param max.iters The maximum number of iterations to perform; choosing a
#'   lower number may result in clustering that is further away from the
#'   observed graph's. Default: 100
#' @param ... Other arguments passed to \code{\link{make_brainGraph}}
#' @inheritParams Creating_Graphs
#' @export
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#'
#' @return \code{sim.rand.graph.par} - a \emph{list} of \emph{N} random graphs
#'   with some additional vertex and graph attributes
#'
#' @rdname random_graphs
#' @seealso \code{\link[igraph]{rewire}, \link[igraph]{sample_degseq},
#'   \link[igraph]{keeping_degseq}}
#' @examples
#' \dontrun{
#' rand1 <- sim.rand.graph.par(g[[1]][[N]], N=1e3)
#' rand1.cl <- sim.rand.graph.par(g[[1]][[N]], N=1e2,
#'   clustering=T, max.iters=1e3)
#' }

sim.rand.graph.par <- function(g, level=c('subject', 'group'), N=100L,
                               clustering=FALSE, rewire.iters=max(10*ecount(g), 1e4L),
                               cl=g$transitivity, max.iters=100L, ...) {
  stopifnot(is_igraph(g))
  level <- match.arg(level)
  if (!getDoParRegistered()) {
    clust <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(clust)
  }
  if (isTRUE(clustering)) {
    r <- foreach(i=seq_len(N), .packages=c('igraph', 'brainGraph'), .export=g$atlas) %dopar% {
      tmp <- sim.rand.graph.clust(g, rewire.iters, cl, max.iters)
      tmp$atlas <- g$atlas
      tmp <- make_brainGraph(tmp, type='random', level=level, set.attrs=TRUE, ...)
      tmp
    }
  } else {
    if (is_connected(g)) {
      V(g)$degree <- degree(g)
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- sample_degseq(V(g)$degree, method='vl')
        tmp$Group <- g$Group
        tmp$name <- g$name
        tmp$threshold <- g$threshold
        tmp$density <- graph.density(tmp)
        tmp <- make_brainGraph(tmp, type='random', level=level, set.attrs=TRUE, ...)
        tmp
      }
    } else {
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- rewire(g, keeping_degseq(loops=FALSE, rewire.iters))
        tmp <- make_brainGraph(tmp, type='random', level=level, set.attrs=TRUE, ...)
        tmp
      }
    }
  }
}

#' Simulate a random graph with given degree sequence and clustering
#'
#' \code{sim.rand.graph.clust} simulates a random graph with a given degree
#' sequence \emph{and} clustering coefficient. Increasing the \code{max.iters}
#' value will result in a closer match of clustering with the observed graph.
#'
#' @export
#' @return \code{sim.rand.graph.clust} - A single \code{igraph} graph object
#'
#' @rdname random_graphs
#' @seealso \code{\link[igraph]{transitivity}}
#' @references Bansal, S. and Khandelwal, S. and Meyers, L.A. (2009) Exploring
#'   biological network structure with clustered random networks. \emph{BMC
#'   Bioinformatics}, \bold{10}, 405--421.
#'   \doi{10.1186/1471-2105-10-405}

sim.rand.graph.clust <- function(g, rewire.iters=1e4, cl=g$transitivity, max.iters=100) {
  g <- rewire(g, keeping_degseq(loops=FALSE, rewire.iters))
  degs <- degree(g)
  degs.large <- which(degs > 1)

  cur.iter <- 0
  while ((transitivity(g) < cl) & (cur.iter < max.iters)) {
    repeat {
      g.cand <- g

      # If E(y1, y2) and E(z1, z2) don't exist, rewire 2 edges
      repeat {
        e <- choose.edges(g.cand, degs.large)
        if (!are_adjacent(g.cand, e$y1, e$y2) &&
            !are_adjacent(g.cand, e$z1, e$z2)) break
      }

      edges_to_remove <- get.edge.ids(g.cand, c(e$y1, e$z1, e$y2, e$z2), directed=FALSE)
      edges_to_add <- c(e$y1, e$y2, e$z1, e$z2)
      g.cand <- delete_edges(g.cand, edges_to_remove)
      g.cand <- add_edges(g.cand, edges_to_add)

      if (transitivity(g.cand) > transitivity(g)) break
    }

    cur.iter <- cur.iter + 1
    g <- g.cand
  }

  return(g)
}

#' Select edges for re-wiring
#'
#' \code{choose.edges} selects edges to be re-wired when simulating random
#' graphs while controlling for \emph{clustering}. It is based on the algorithm
#' by Bansal et al. (2009), BMC Bioinformatics.
#'
#' @param g A graph object
#' @param degs.large Integer vector of vertex numbers with degree greater than
#'   one
#'
#' @return A data frame with four elements; two edges \code{(y1, z1)} and
#'   \code{(y2, z2)} will be removed, and two edges \code{(y1, y2)} and
#'   \code{(z1, z2)} will be added between the four vertices.
#' @keywords internal
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

choose.edges <- function(g, degs.large) {
  # Uniformly select a random vertex with degree > 1 (x),
  # and 2 of its neighbors (y1 & y2)
  #-----------------------------------------------------------------------------
  get_random_v <- function(g, degs.large) {
    repeat {
      x <- degs.large[sample.int(length(degs.large), 1L)]
      neighb <- intersect(neighbors(g, x), degs.large)
      if (length(neighb) >= 2L) return(list(x, neighb[sample.int(length(neighb), 2L)]))
    }
  }
  #-----------------------------------------------------------------------------
  repeat {
    tmp <- get_random_v(g, degs.large)
    x <- tmp[[1L]]
    y <- tmp[[2L]]

    # Try to select random neighbors (z1 & z2) of y1 & y2 s.t. z1 != z2 != x
    y1.neighb <- neighbors(g, y[1L])
    choices1 <- setdiff(y1.neighb, c(x, y[2L]))
    if (length(choices1) == 0L) {
      next
    } else {
      z1 <- choices1[sample.int(length(choices1), 1L)]
    }

    y2.neighb <- neighbors(g, y[2L])
    choices2 <- setdiff(y2.neighb, c(x, y[1L], z1))
    if (length(choices2) > 0L) break
  }

  z2 <- choices2[sample.int(length(choices2), 1L)]
  return(data.frame(y1=y[1L], y2=y[2L], z1, z2))
}

#' Hirschberger-Qi-Steuer null covariance matrix generation
#'
#' \code{sim.rand.graph.hqs} generates a number of random covariance matrices
#' using the Hirschberger-Qi-Steuer (HQS) algorithm, and create graphs from
#' those matrices.
#'
#' \code{sim.rand.graph.hqs} - The first step is to create the observed
#' covariance of residuals (or whatever matrix/data.table is provided). Then
#' random covariance matrices are created with the same distributional
#' properties as the observed matrix, they are converted to correlation
#' matrices, and finally graphs from these matrices. By default, weighted graphs
#' will be created in which the edge weights represent correlation values. If
#' you want binary matrices, you must provide a correlation threshold.
#'
#' @param resids A \code{brainGraph_resids} object, a \code{data.table} of
#'   residuals, or a numeric matrix
#' @param weighted Logical indicating whether to create weighted graphs. If
#'   true, a threshold must be provided.
#' @param r.thresh Numeric value for the correlation threshold, if
#'   \code{weighted=FALSE}.
#' @export
#' @inheritParams Creating_Graphs
#' @importFrom foreach getDoParRegistered
#' @importFrom doParallel registerDoParallel
#' @return \code{sim.rand.graph.hqs} - A list of random graphs from the null
#'   covariance matrices
#'
#' @rdname random_graphs
#' @references Hirschberger M., Qi Y., Steuer R.E. (2007) Randomly generating
#'   portfolio-selection covariance matrices with specified distributional
#'   characteristics. \emph{European Journal of Operational Research}.
#'   \bold{177}, 1610--1625. \doi{10.1016/j.ejor.2005.10.014}

sim.rand.graph.hqs <- function(resids, level=c('subject', 'group'), N=100L,
                               weighted=TRUE, r.thresh=NULL, ...) {
  level <- match.arg(level)
  if (inherits(resids, 'brainGraph_resids')) {
    sID <- getOption('bg.subject_id')
    gID <- getOption('bg.group')
    resids <- resids$resids.all[, !c(sID, gID), with=FALSE]
  }
  A <- cov.wt(resids)$cov
  a <- 1 / sqrt(diag(A))

  e <- mean(c(A[upper.tri(A)], A[lower.tri(A)]), na.rm=TRUE)
  v <- var(c(A[upper.tri(A)], A[lower.tri(A)]), na.rm=TRUE)
  ebar <- mean(diag(A), na.rm=TRUE)
  n <- dim(A)[1L]

  m <- max(2, floor((ebar^2 - e^2) / v))
  ehat <- sqrt(e / m)
  vhat <- sqrt(ehat^4 + v / m) - ehat^2

  # Returns a "n X n" matrix
  hqs <- function(ehat, vhat, n, m) {
    u <- matrix(0, nrow=n, ncol=m)
    for (i in seq_len(m)) {
      u[, i] <- runif(n)
    }
    Q <- qnorm(u)
    f <- ehat + sqrt(vhat)*Q
    sig <- tcrossprod(f)
    return(sig)
  }

  if (!getDoParRegistered()) {
    cl <- makeCluster(getOption('bg.ncpus'))
    registerDoParallel(cl)
  }
  if (isTRUE(weighted)) {
    g <- foreach(i=seq_len(N)) %dopar% {
      S <- hqs(ehat, vhat, n, m)
      R <- outer(a, a) * S  # Convert to a correlation matrix
      make_brainGraph(R, type='random', level=level, set.attrs=TRUE, weighted=TRUE, ...)
    }
  } else {
    if (is.null(r.thresh)) stop('Please provide a value for thresholding')
    g <- foreach(i=seq_len(N)) %dopar% {
      S <- hqs(ehat, vhat, n, m)
      R <- outer(a, a) * S  # Convert to a correlation matrix
      R.thresh <- ifelse(R > r.thresh, 1, 0)
      make_brainGraph(R.thresh, type='random', level=level, set.attrs=TRUE, ...)
    }
  }
  return(g)
}
