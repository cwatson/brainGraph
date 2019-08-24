#' Rich club calculations
#'
#' \code{rich_club_coeff} calculates the \emph{rich club} of a graph, returning
#' the rich-club coefficient, \eqn{\phi}, and the subgraph of rich club
#' vertices.
#'
#' @param g An \code{igraph} graph object
#' @param k Integer; the minimum degree for including a vertex (default: 1)
#' @param weighted Logical indicating whether or not edge weights should be
#'   used (default: \code{FALSE})
#' @export
#'
#' @return \code{\link{rich_club_coeff}} - a list with components:
#'   \item{phi}{The rich club coefficient, \eqn{\phi}.}
#'   \item{graph}{A subgraph containing only the rich club vertices.}
#'   \item{Nk}{The number of vertices in the rich club graph.}
#'   \item{Ek}{The number of edges in the rich club graph.}
#'
#' @name RichClub
#' @aliases rich_club_coeff
#' @rdname rich_club
#' @family Rich-club functions
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Zhou, S. and Mondragon, R.J. (2004) The rich-club phenomenon
#'   in the internet topology. \emph{IEEE Comm Lett}, \bold{8}, 180--182.
#'   \url{https://dx.doi.org/10.4018/978-1-59140-993-9.ch066}
#' @references Opsahl, T. and Colizza, V. and Panzarasa, P. and Ramasco, J.J.
#'   (2008) Prominence and control: the weighted rich-club effect. \emph{Physical
#'   Review Letters}, \bold{101.16}, 168702.
#'   \url{https://dx.doi.org/10.1103/PhysRevLett.101.168702}

rich_club_coeff <- function(g, k=1, weighted=FALSE) {
  stopifnot(is_igraph(g))
  degs <- check_degree(g)
  Nv <- vcount(g)
  Nk <- sum(degs > k)
  if (Nk == 0) {
    return(list(phi=NaN, graph=make_empty_graph(), Nk=0, Ek=0))
  } else {
    rich.club.verts <- order(degs)[(Nv - Nk + 1):Nv]
    rich.club.graph <- induced_subgraph(g, rich.club.verts)
    Ek <- ecount(rich.club.graph)

    if (isTRUE(weighted)) {
      Wr <- sum(E(rich.club.graph)$weight)
      weights <- sort(E(g)$weight, decreasing=TRUE)[1:Ek]
      phi <- Wr / sum(weights)
    } else {
      phi <- graph.density(rich.club.graph)
    }

    return(list(phi=phi, graph=rich.club.graph, Nk=Nk, Ek=Ek))
  }
}

#' Rich club coefficient for all degrees in a graph
#'
#' \code{rich_club_all} is a wrapper for \code{\link{rich_club_coeff}} that
#' calculates the rich-club coefficient for all degrees present in the graph. It
#' returns a \code{data.table} with the coefficients and vertex and edge counts
#' for each successive rich club.
#' @export
#' @return \code{\link{rich_club_all}} - a \code{data.table} with components:
#'   \item{k}{A vector of all vertex degrees present in the original graph}
#'   \item{phi}{The rich-club coefficient}
#'   \item{Nk}{The number of vertices in the rich club for each successive
#'     \emph{k}}
#'   \item{Ek}{The number of edges in the rich club for each successive
#'     \emph{k}}
#'
#' @aliases rich_club_all
#' @rdname rich_club

rich_club_all <- function(g, weighted=FALSE) {
  stopifnot(is_igraph(g))
  k <- check_degree(g)
  R <- lapply(1:max(k), function(x) rich_club_coeff(g, x, weighted))
  phi <- vapply(R, with, numeric(1), phi)
  Nk <- vapply(R, with, numeric(1), Nk)
  Ek <- vapply(R, with, numeric(1), Ek)
  dt.rich <- data.table(k=1:max(k), phi=phi, Nk=Nk, Ek=Ek)
  return(dt.rich)
}

#' Calculate the normalized rich club coefficient
#'
#' \code{rich_club_norm} will (optionally) generate a number of random graphs,
#' calculate their rich club coefficients (\eqn{\phi}), and return
#' \eqn{\phi_{norm}} of the graph of interest, which is the observed rich-club
#' coefficient divided by the mean across the random graphs.
#'
#' If random graphs have already been generated, you can supply a list as an
#' argument (since graph generation is time consuming).
#'
#' @param N Integer; the number of random graphs to generate (default: 100)
#' @param rand A list of \code{igraph} graph objects, if random graphs have
#'   already been generated (default: \code{NULL})
#' @param ... Other parameters (passed to \code{\link{sim.rand.graph.par}})
#' @export
#'
#' @return \code{\link{rich_club_norm}} - a data table with columns:
#'   \item{k}{Sequence of degrees}
#'   \item{rand}{Rich-club coefficients for the random graphs}
#'   \item{orig}{Rich-club coefficients for the original graph.}
#'   \item{norm}{Normalized rich-club coefficients.}
#'   \item{p}{The P-values based on the distribution of rich-club coefficients
#'     from the random graphs.}
#'   \item{p.fdr}{The FDR-adjusted P-values}
#'   \item{density}{The observed graph's density}
#'   \item{threshold}{}
#'   \item{Group}{}
#'   \item{name}{}
#'
#' @aliases rich_club_norm
#' @rdname rich_club
#' @family Random graph functions
#'
#' @references Colizza, V. and Flammini, A. and Serrano, M.A. and Vespignani, A.
#'   (2006) Detecting rich-club ordering in complex networks. \emph{Nature
#'   Physics}, \bold{2}, 110--115. \url{https://dx.doi.org/10.1038/nphys209}

rich_club_norm <- function(g, N=1e2, rand=NULL, ...) {
  k <- orig <- p <- p.fdr <- NULL
  stopifnot(is_igraph(g))
  degs <- check_degree(g)
  if (!'rich' %in% graph_attr_names(g)) g$rich <- rich_club_all(g)
  if (is.null(rand)) {
    rand <- sim.rand.graph.par(g, N, ...)
  } else {
    if (!all(vapply(rand, is_igraph, logical(1)))) {
      stop('Argument "rand" must be a list of igraph graph objects!')
    }
    N <- length(rand)
  }

  phi.rand <- t(sapply(rand, function(x) x$rich$phi))
  max.deg <- max(degs)
  DT <- data.table(k=rep(seq_len(max.deg), each=N),
                   rand=c(phi.rand),
                   orig=rep(g$rich$phi, each=N))
  DT[, norm := unique(orig) / mean(rand), by=k]
  DT[, p := (sum(rand >= unique(orig)) + 1) / (N + 1), by=k]
  dt.phi <- DT[, .SD[1], by=k]
  dt.phi[, p.fdr := p.adjust(p, 'fdr')]
  dt.phi[, rand := NULL]
  dt.phi$rand <- DT[, mean(rand), by=k]$V1
  if ('density' %in% graph_attr_names(g)) {
    dt.phi$density <- g$density
  } else {
    dt.phi$density <- graph.density(g)
  }
  if ('threshold' %in% graph_attr_names(g)) dt.phi$threshold <- g$threshold
  if ('Group' %in% graph_attr_names(g)) dt.phi$Group <- g$Group
  if (!is.null(g$name)) dt.phi$Study.ID <- g$name
  return(dt.phi)
}

#' Assign graph attributes based on rich-club analysis
#'
#' This function will assign vertex- and edge-level attributes based on the
#' results of a \emph{rich-club} analysis, based on a range of vertex degrees in
#' which the rich-club coefficient was determined to be significantly greater
#' than that of a set of random graphs (see \code{\link{rich_club_norm}}).
#'
#' Vertices which are in the rich club will be assigned an attribute
#' \code{rich}, taking on a binary value. Their colors (attribute
#' \code{color.rich}) will be either \emph{red} or \emph{gray}. Their sizes
#' (attribute \code{size.rich}) will either be 10 or will be proportional to
#' their degree.
#'
#' Edge attribute \code{type.rich} takes on three values: \emph{rich-club} (if
#' it connects two rich-club vertices), \emph{feeder} (if it connects a rich- to
#' a non-rich-club vertex), and \emph{local} (if it connects two non-rich-club
#' vertices). They will also be given a \code{color.rich} attribute (either
#' \emph{red}, \emph{orange}, or \emph{green}). Edge sizes (\code{size.rich})
#' will be largest for \emph{rich-club} connections, then smaller for
#' \emph{feeder}, and smallest for \emph{local}.
#'
#' @param g An \code{igraph} graph object
#' @param deg.range Integer vector of the range of degrees indicating
#'   inclusion in the rich-club; if the default \emph{NULL}, it will be from 1
#'   to the maximum degree in the graph
#' @param adj.vsize Logical indicating whether to adjust vertex size
#'   proportional to degree (default: FALSE)
#' @export
#'
#' @return An \code{igraph} graph object with additional attributes:
#'   \item{rich}{Binary indicating membership in the rich-club}
#'   \item{type.rich}{Edge attribute indicating the type of connection}
#'   \item{color.rich}{Edge and vertex attributes}
#'   \item{size.rich}{Edge and vertex attributes}
#'
#' @family Rich-club functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' g <- rich_club_attrs(g, rich.dt[density == densities[N] & p.fdr < .01,
#'                                 range(k)])
#' }

rich_club_attrs <- function(g, deg.range=NULL, adj.vsize=FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(deg.range)) deg.range <- c(1, max(V(g)$degree))
  V(g)$rich <- 0
  V(g)[V(g)$degree >= deg.range[1] & V(g)$degree <= deg.range[2]]$rich <- 1
  E(g)[which(V(g)$rich == 1) %--% which(V(g)$rich == 1)]$type.rich <- 'rich-club'
  E(g)[which(V(g)$rich == 1) %--% which(V(g)$rich == 0)]$type.rich <- 'feeder'
  E(g)[which(V(g)$rich == 0) %--% which(V(g)$rich == 0)]$type.rich <- 'local'

  V(g)[V(g)$rich == 1]$color.rich <- 'red'
  V(g)[V(g)$rich == 0]$color.rich <- 'gray'
  E(g)[E(g)$type.rich == 'rich-club']$color.rich <- 'red'
  E(g)[E(g)$type.rich == 'feeder']$color.rich <- 'orange'
  E(g)[E(g)$type.rich == 'local']$color.rich <- 'green'

  if (isTRUE(adj.vsize)) {
    V(g)[V(g)$rich == 1]$size.rich <- vec.transform(V(g)[which(V(g)$rich == 1)]$degree, 3, 15)
  } else {
    V(g)[V(g)$rich == 1]$size.rich <- 10
  }
  V(g)[V(g)$rich == 0]$size.rich <- 0
  E(g)[E(g)$type.rich == 'rich-club']$size.rich <- 3.5
  E(g)[E(g)$type.rich == 'feeder']$size.rich <- 1.5
  E(g)[E(g)$type.rich == 'local']$size.rich <- 0.5

  return(g)
}

#' Calculate the rich core of a graph
#'
#' \code{rich_core} finds the boundary of the rich core of a graph, based on the
#' decreasing order of vertex degree. It also calculates the degree that
#' corresponds to that rank, and the core size relative to the total number of
#' vertices in the graph.
#'
#' For weighted graphs, the degree is substituted by a normalized weight:
#' \deqn{ceiling(A / w_{min})}
#' where \eqn{w_{min}} is the minimum weight (that is greater than 0), and
#' \eqn{ceiling()} is the \emph{ceiling} function that rounds up to the nearest
#' integer.
#'
#' @importFrom Matrix colSums
#' @export
#' @return \code{\link{rich_core}} - a data table with columns:
#'   \item{density}{The density of the graph.}
#'   \item{rank}{The rank of the boundary for the rich core.}
#'   \item{k.r}{The degree/strength of the vertex at the boundary.}
#'   \item{core.size}{The size of the core relative to the graph size.}
#'   \item{weighted}{Whether or not weights were used}
#'
#' @aliases rich_core
#' @rdname rich_club
#'
#' @references Ma, A and Mondragon, R.J. (2015) Rich-cores in networks.
#'   \emph{PLoS One}, \bold{10(3)}, e0119678.
#'   \url{https://dx.doi.org/10.1371/journal.pone.0119678}

rich_core <- function(g, weighted=FALSE) {
  stopifnot(is_igraph(g))
  if (isTRUE(weighted)) {
    A <- as_adj(g, names=FALSE, attr='weight')
    w_min <- min(summary(A)$x)
    A <- ceiling(A / w_min)
  } else {
    A <- as_adj(g, names=FALSE, sparse=FALSE)
  }
  degs <- colSums(A)
  vorder <- order(degs, decreasing=TRUE)
  kplus <- vapply(vorder, function(x) sum(A[x, which(degs > degs[x])]), double(1))

  dens <- ifelse('density' %in% graph_attr_names(g), g$density, graph.density(g))
  r <- max(which(kplus == max(kplus)))
  k.r <- as.integer(degs[vorder][r])
  core.size <- r / ncol(A)
  return(data.table(density=dens, rank=r, k.r=k.r, core.size=core.size, weighted=weighted))
}
