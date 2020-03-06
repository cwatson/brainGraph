#' Calculate correlation matrix and threshold
#'
#' \code{corr.matrix} calculates the correlation between all column pairs of a
#' given data frame, and thresholds the resultant correlation matrix based on a
#' given density (e.g., \code{0.1} if you want to keep only the 10\% strongest
#' correlations). If you want to threshold by a specific correlation coefficient
#' (via the \code{thresholds} argument), then the \code{densities} argument is
#' ignored.
#'
#' If you wish to exclude regions from your analysis, you can give the indices
#' of their columns with the \code{exclude.reg} argument.
#'
#' By default, the Pearson correlation coefficients are calculated, but you can
#' return Spearman by changing the \code{type} argument.
#'
#' @param resids An object of class \code{brainGraph_resids} (the output from
#'   \code{\link{get.resid}})
#' @param densities Numeric vector indicating the resultant network
#'   density(ies); keeps the top \emph{X}\% of correlations
#' @param thresholds Numeric; absolute correlation value to threshold by
#'   (default: \code{NULL})
#' @param what Character string indicating whether to correlate the residuals or
#'   the raw structural MRI values (default: \code{'resids'})
#' @param exclude.reg Character vector of regions to exclude (default:
#'   \code{NULL})
#' @param type Character string indicating which type of correlation coefficient
#'   to calculate (default: \code{'pearson'})
#' @param rand Logical indicating whether the function is being called for
#'   permutation testing; not intended for general use (default: \code{FALSE})
#' @export
#'
#' @return A \code{corr_mats} object containing the following components:
#'   \item{R,P}{Numeric arrays of correlation coefficients and P-values. The
#'     length of the 3rd dimension equals the number of groups}
#'   \item{r.thresh}{A list of 3-d binary arrays indicating correlations that
#'     are above a certain threshold. The length of the list equals the number
#'     of groups, and the length of the 3rd dimension equals the number of
#'     thresholds/densities.}
#'   \item{thresholds}{Numeric matrix of the thresholds supplied. The number of
#'     columns equals the number of groups.}
#'   \item{what}{Residuals or raw values}
#'   \item{exclude.reg}{Excluded regions (if any)}
#'   \item{type}{Pearson or Spearman}
#'   \item{atlas}{The brain atlas used}
#'   \item{densities}{Numeric vector; the densities of the resulting graphs, if
#'     you chose to threshold each group to have equal densities.}
#'
#' @rdname correlation_matrices
#' @family Structural covariance network functions
#' @seealso \code{\link[Hmisc]{rcorr}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' myResids <- get.resid(lhrh, covars)
#' corrs <- corr.matrix(myResids, densities=densities)))
#' }

corr.matrix <- function(resids, densities, thresholds=NULL, what=c('resids', 'raw'),
                        exclude.reg=NULL, type=c('pearson', 'spearman'), rand=FALSE) {
  if (!requireNamespace('Hmisc', quietly=TRUE)) stop('You must install "Hmisc" to use this function.')
  stopifnot(inherits(resids, 'brainGraph_resids'))
  gID <- getOption('bg.group')
  N <- nregions(resids)
  regions <- region.names(resids)

  sID <- getOption('bg.subject_id')
  # Different behavior if called for permutation testing
  if (isTRUE(rand)) {
    res.all <- as.matrix(resids$resids.all[, !c(sID, gID), with=FALSE])
    corrs <- Hmisc::rcorr(res.all)
    r <- corrs$r
    emax <- N  * (N - 1) / 2
    thresholds <- get_thresholds(r, densities, emax)
    r.thresh <- array(0, dim=c(N, N, length(thresholds)), dimnames=list(regions, regions, NULL))
    for (i in seq_along(thresholds)) r.thresh[, , i] <- ifelse(r > thresholds[i], 1, 0)
    return(list(list(R=r, r.thresh=r.thresh)))
  }

  what <- match.arg(what)
  type <- match.arg(type)
  if (what == 'resids') {
    res.all <- resids$resids.all[, !sID, with=FALSE]
  } else if (what == 'raw') {
    res.all <- dcast(resids$all.dat.long, paste(sID, '+', gID, '~ Region'))
    setkeyv(res.all, c(gID, sID))
    res.all <- res.all[, !sID, with=FALSE]
  }
  if (!is.null(exclude.reg)) {
    res.all <- res.all[, -exclude.reg, with=FALSE]
    regions <- setdiff(regions, exclude.reg)
    N <- N - length(exclude.reg)
  }

  # Loop through the groups
  grps <- unique(groups(resids))
  kNumGroups <- length(grps)
  r <- p <- array(0, dim=c(N, N, kNumGroups), dimnames=list(regions, regions, grps))
  r.thresh <- setNames(vector('list', kNumGroups), grps)
  if (!is.null(thresholds)) {
    kNumThresh <- length(thresholds)
    thresh.mat <- matrix(rep.int(thresholds, kNumGroups),
                         nrow=kNumThresh, ncol=kNumGroups,
                         dimnames=list(NULL, grps))
  } else {
    kNumThresh <- length(densities)
    thresh.mat <- matrix(0, nrow=kNumThresh, ncol=kNumGroups,
                         dimnames=list(NULL, grps))
  }
  for (g in grps) {
    corrs <- Hmisc::rcorr(as.matrix(res.all[g, !gID, with=FALSE]), type=type)
    r[, , g] <- corrs$r
    p[, , g] <- corrs$P

    # Calculate a threshold so that "density"% of possible connections are present
    if (hasArg('densities')) {
      tmp <- corrs$r
      emax <- N  * (N - 1) / 2
      if (any(densities == 1)) {
        i <- which(densities == 1)
        thresh.mat[i, g] <- min(tmp, na.rm=TRUE)
        thresh.mat[-i, g] <- get_thresholds(tmp, densities[-i], emax)
      } else {
        thresh.mat[, g] <- get_thresholds(tmp, densities, emax)
      }
    }
    r.thresh[[g]] <- array(0, dim=c(N, N, kNumThresh), dimnames=list(regions, regions))
    for (i in seq_along(thresh.mat[, g])) {
      r.thresh[[g]][, , i] <- ifelse(r[, , g] > thresh.mat[i, g], 1, 0)
    }
  }
  out <- list(R=r, P=p, r.thresh=r.thresh, thresholds=thresh.mat, what=what,
              exclude.reg=exclude.reg, type=type, atlas=resids$atlas)
  if (hasArg('densities')) out <- c(out, list(densities=densities))
  class(out) <- c('corr_mats', class(out))
  return(out)
}

#' Subset correlation matrix objects
#'
#' @param x,object A \code{corr_mats} object
#' @param i Integer for subsetting by density/threshold
#' @param g Integer, character, or logical for subsetting by group
#' @export
#'
#' @name Extract.corr_mats
#' @rdname correlation_matrices

`[.corr_mats` <- function(x, i, g=NULL) {
  if (!is.null(g)) {
    grps <- groups(x)
    kNumGroups <- length(grps)
    if (is.logical(g) && length(g) != kNumGroups) {
      stop('Logical indexing vector must be of the same length as the number of groups.')
    }
    if (length(g) > kNumGroups) {
      warning('Indexing vector cannot have length greater than the number of groups.')
      g <- g[seq_len(kNumGroups)]
    }
    g <- switch(class(g),
                numeric=, logical=grps[g],
                character=which(grps %in% g))
    x$R <- x$R[, , g, drop=FALSE]
    x$P <- x$P[, , g, drop=FALSE]
    x$r.thresh <- x$r.thresh[g]
    x$thresholds <- x$thresholds[, g, drop=FALSE]
  }
  if (missing(i)) i <- seq_len(dim(x$thresholds)[1L])
  for (g in groups(x)) x$r.thresh[[g]] <- x$r.thresh[[g]][, , i, drop=FALSE]
  x$thresholds <- x$thresholds[i, , drop=FALSE]
  if (hasName(x, 'densities')) x$densities <- x$densities[i]
  return(x)
}

#' Plot raw or thresholded correlation matrices
#'
#' The \code{plot} method will plot \dQuote{heat maps} of the correlation
#' matrices.
#'
#' @section Plotting correlation matrices:
#' There are several ways to control the plot appearance. First, you may plot
#' the \dQuote{raw} correlations, or only those of the thresholded (binarized)
#' matrices. Second, you may order the vertices by a given vertex attribute; by
#' default, they will be ordered by \emph{lobe}, but you may also choose to
#' order by, e.g., \emph{network} (for the \code{dosenbach160} atlas) or by
#' \emph{community membership}. In the latter case, you need to pass a
#' \code{brainGraphList} object to the \code{graphs} argument; each graph in the
#' object must have a vertex attribute specified in \code{order.by}. Finally,
#' you can control the legend text with \code{grp.names}.
#'
#' @param mat.type Character string indicating whether to plot raw or thresholded
#'   (binarized) matrices. Default: \code{'raw'}
#' @param thresh.num Integer specifying which threshold to plot (if
#'   \code{mat.type='thresholded'}). Default: \code{1L}
#' @param ordered Logical indicating whether to order the vertices by some
#'   grouping. Default: \code{TRUE}
#' @param order.by Character string indicating how to group vertices. Default:
#'   \code{'lobe'}
#' @param graphs A \code{brainGraphList} object containing graphs with the
#' vertex-level attribute of interest. Default: \code{NULL}
#' @param grp.names Character vector specifying the names of each group of
#'   vertices. Default: \code{NULL}
#' @param legend.title Character string for the legend title. Default is to
#'   leave blank
#' @param ... Unused
#' @export
#' @rdname correlation_matrices
#' @examples
#' \dontrun{
#' corrs <- corr.matrix(myResids, densities)
#' plot(corrs, order.by='comm', graphs=g.list, grp.names='Community')
#' }

plot.corr_mats <- function(x, mat.type=c('thresholded', 'raw'), thresh.num=1L,
                           ordered=TRUE, order.by='lobe', graphs=NULL,
                           grp.names=NULL, legend.title=NULL, ...) {
  Var1 <- Var2 <- Group <- value <- memb1 <- memb2 <- memb <- col.text <- NULL
  grps <- groups(x)
  regions <- region.names(x)
  kNumVertices <- length(regions)
  base_size <- if (kNumVertices > 90L) 7.5 else 9
  if (is.null(legend.title) && isTRUE(ordered)) {
    legend.title <- switch(order.by,
      comp='Connected\nComponent', hemi='Hemisphere', comm=, comm.wt='Community (#)',
      class='Tissue class', simpleCap(order.by))
  }

  type <- match.arg(mat.type)
  legend.pos <- if (type == 'raw' || isTRUE(ordered)) 'right' else 'none'
  mytheme <- theme(legend.position=legend.pos,
    axis.text.x=element_text(size=0.7*base_size, angle=45),
    axis.text.y=element_text(size=0.7*base_size),
    axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
    plot.title=element_text(hjust=0.5, face='bold'))
  matplots <- setNames(vector('list', length(grps)), grps)
  if (type == 'raw') {
    legend.title <- paste0('Corr. coeff.\n(', simpleCap(x$type), ')')
    mats <- as.data.table(x$R, sorted=FALSE)
    setnames(mats, c('Var1', 'Var2', 'Group', 'value'))
    mats[, Var1 := factor(Var1, levels=regions)]
    mats[, Var2 := factor(Var2, levels=regions)]
    for (g in grps) {
      mats.m <- mats[Group == g]
      matplots[[g]] <- ggplot(mats.m, aes(Var1, Var2, fill=value)) +
        geom_tile() +
        scale_fill_gradient2(low='white', high='red') + mytheme +
        labs(title=g, fill=legend.title) + ylim(rev(levels(mats.m$Var2)))
    }
    return(matplots)
  } else {
    mats <- sapply(x$r.thresh, function(m) m[, , thresh.num], simplify='array')
  }

  if (isTRUE(ordered)) {
    if (!hasName(get(x$atlas), order.by)) {
      if (is.brainGraphList(graphs)) {
        g.list <- graphs[]
      } else if (is.brainGraphList(graphs[[thresh.num]])) {
        g.list <- graphs[[thresh.num]][]
      } else {
        stop('Please provide a (list of) "brainGraphList" object(s).')
      }
      order.mat <- vapply(g.list, vertex_attr, numeric(kNumVertices), order.by)
      if (is.null(grp.names)) grp.names <- ''
      if (length(grp.names) == 1L) grp.names <- paste(grp.names, seq_len(max(order.mat)))
    } else {
      tmp <- get(x$atlas)[, get(order.by)]
      order.mat <- matrix(rep.int(as.integer(tmp), 2L), ncol=2L)
      grp.names <- levels(tmp)
    }

    kNumCols <- apply(order.mat, 2L, max)
    cols.new <- lapply(kNumCols, function(y) c(group.cols[seq_len(y)], 'gray50', 'white'))
    if (order.by == 'comp') cols.new <- lapply(cols.new, function(y) setdiff(y, 'gray50'))
    grp.names <- c(grp.names, 'Inter', '')
    new.order <- apply(order.mat, 2L, order)

  } else {
    new.order <- order.mat <- matrix(rep.int(seq_len(kNumVertices), 2L), ncol=2L)
    cols.new <- lapply(grps, function(y) c('white', 'gray50'))
  }
  dimnames(order.mat) <- dimnames(new.order) <- list(regions, grps)
  names(cols.new) <- grps

  mats.ord <- lapply(grps, function(y) array(mats[new.order[, y], new.order[, y], y],
                                             dim=c(kNumVertices, kNumVertices, 1L)))
  names(mats.ord) <- grps
  for (g in grps) {
    dimnames(mats.ord[[g]]) <- list(regions[new.order[, g]], regions[new.order[, g]], g)
    mats.m <- as.data.table(mats.ord[[g]], sorted=FALSE)
    setnames(mats.m, c('Var1', 'Var2', 'Var3', 'value'))
    mats.m[, Var1 := factor(Var1, levels=regions[new.order[, g]])]
    mats.m[, Var2 := factor(Var2, levels=regions[new.order[, g]])]
    if (isTRUE(ordered)) {
      mats.m[, memb := '']
      mats.m[, memb1 := grp.names[order.mat[as.character(Var1), g]]]
      mats.m[, memb2 := grp.names[order.mat[as.character(Var2), g]]]
      mats.m[value == 1, memb := ifelse(memb1 == memb2, memb1, 'Inter')]
      mats.m[, memb := factor(memb, levels=grp.names)]
      mats.m[, memb1 := factor(memb1, levels=grp.names)]
      mats.m[, col.text := cols.new[[g]][as.integer(memb1)]]
    } else {
      mats.m[, memb := factor(value)]
      mats.m[, col.text := 'black']
    }

    mytheme$axis.text.x <- element_text(size=0.7*base_size, angle=45,
                                        color=mats.m[seq_len(kNumVertices), col.text])
    mytheme$axis.text.y <- element_text(size=0.7*base_size,
                                        color=mats.m[rev(seq_len(kNumVertices)), col.text])
    matplots[[g]] <- ggplot(mats.m, aes(Var1, Var2, fill=memb)) +
      geom_tile() +
      scale_fill_manual(values=cols.new[[g]]) + mytheme +
      labs(title=g, fill=legend.title) + ylim(rev(levels(mats.m$Var2)))
  }
  return(matplots)
}
