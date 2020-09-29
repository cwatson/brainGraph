#' Plot normalized rich club coefficients against degree threshold
#'
#' Returns a line plot of the normalized rich club coefficient. Optionally, can
#' include a shaded region demarcating the \code{\link{rich_core}} cutoff (if
#' you supply a list of graph objects to the \code{g} argument).
#'
#' @param rich.dt A \code{data.table} with rich-club coefficients
#' @param facet.by A character string indicating whether the variable of
#' interest is \dQuote{density} or \dQuote{threshold} (e.g. with DTI data)
#' @param densities A numeric vector of the densities to plot
#' @param alpha The significance level. Default: \code{0.05}
#' @param fdr A logical, indicating whether or not to use the FDR-adjusted
#'   p-value for determining significance. Default: \code{TRUE}
#' @param g.list A list \code{brainGraphList} objects; required if you want to
#'   plot a shaded region demarcating the \code{\link{rich_core}}
#' @param smooth Logical indicating whether or not to plot a smooth curve when
#'   data from multiple subjects (per group) are present. Default: \code{TRUE}.
#'   Ignored for group-level data.
#' @export
#'
#' @return A \code{trellis} or \code{ggplot} object
#'
#' @family Rich-club functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot_rich_norm(rich.dt, facet.by='density', densities[N:(N+1)], g=g)
#' }

plot_rich_norm <- function(rich.dt, facet.by=c('density', 'threshold'),
                           densities, alpha=0.05, fdr=TRUE, g.list=NULL, smooth=TRUE) {
  yloc <- norm <- xstart <- xend <- panel.num <- NULL
  gID <- getOption('bg.group')

  facet.by <- match.arg(facet.by)
  subDT <- rich.dt[round(get(facet.by), 2L) %in% round(densities, 2L)]
  pvar <- if (isTRUE(fdr)) 'p.fdr' else 'p'
  subDT[, star := ifelse(get(pvar) < alpha, '*', '')]
  subDT[, yloc := extendrange(norm)[1L], by=facet.by]
  subDT[, eval(gID) := as.factor(get(gID))]
  grps <- subDT[, levels(get(gID))]
  if (length(grps) > 1L) {
    for (i in 2L:length(grps)) {
      subDT[get(gID) == grps[i], yloc := yloc - i * 0.05 * diff(range(norm, na.rm=TRUE)), by=facet.by]
    }
  }
  setkeyv(subDT, c(facet.by, gID))

  rects <- data.table(density=subDT[, unique(get(facet.by))], xstart=0L, xend=0L,
                      Group=rep(subDT[, unique(get(gID))],
                                each=length(densities)))
  setnames(rects, 'Group', gID)
  if (!is.null(g.list)) {
    # Check if components are 'brainGraphList' objects
    matches <- vapply(g.list, is.brainGraphList, logical(1L))
    if (any(!matches)) stop("Input must be a list of 'brainGraphList' objects.")

    densities.g <- round(vapply(g.list, function(x) graph_attr(x[1L], 'density'), numeric(1L)),  2L)
    densities.g <- which(densities.g %in% round(densities, 2L))
    g.list <- g.list[densities.g]
    k <- vapply(g.list, function(y) vapply(y[], function(x) rich_core(x)$k.r, integer(1L)), integer(nobs(g.list[[1L]])))
    max.k <- apply(as.matrix(k), 1L, max)

    rects[, xstart := rep(max.k, each=length(densities))]
    rects[, xend := rep(subDT[, max(k), by=facet.by]$V1, each=length(densities))]
  }
  setnames(rects, 'density', facet.by)
  setkeyv(rects, key(subDT))
  rects[, eval(facet.by) := factor(get(facet.by), labels=densities)]
  subDT[, eval(facet.by) := factor(get(facet.by), labels=densities)]
  setkeyv(subDT, c(facet.by, gID))
  setkeyv(rects, c(facet.by, gID))
  rects <- subDT[rects]

  sID <- getOption('bg.subject_id')
  smooth <- FALSE
  if (hasName(subDT, sID) && subDT[, !all.equal(as.character(get(gID)), get(sID))]) {
    smooth <- TRUE
  }

  # 'base' plotting
  #-----------------------------------------------------------------------------
  if (!requireNamespace('ggplot2', quietly=TRUE)) {
    ymin <- rects[, min(yloc)]
    ymax <- rects[, max(norm)]
    ylimit <- extendrange(c(ymin, ymax))
    ylabel <- expression(phi[norm])
    rects[, panel.num := as.numeric(get(facet.by))]
    if (isTRUE(smooth)) {
      plot.type <- 'smooth'
      panelfun <- function(x, y, groups, xleft, xright, ...) {
        panel.num <- NULL
        panel.abline(h=1, lty=2)
        panel.xyplot(x, y, groups, ...)
        rects[panel.num == panel.number() & get(gID) == grps[1L] & star == '*',
              panel.text(k, yloc, labels='*', col=plot.cols[1L])]
        rects[panel.num == panel.number() & get(gID) == grps[2L] & star == '*',
              panel.text(k, yloc, labels='*', col=plot.cols[2L])]
      }
      p <- xyplot(norm ~ k | get(facet.by), groups=get(gID), data=rects,
                  ylim=ylimit, xlab='Degree (k)', ylab=ylabel,
                  type=plot.type, panel=panelfun)
    } else {
      plot.type <- 'l'
      panelfun <- function(x, y, groups, xleft, xright, ...) {
        panel.num <- NULL
        panel.rect(xleft, ybottom=-999, xright, ytop=999, col='lightgrey', border=NULL)
        panel.abline(h=1, lty=2)
        panel.xyplot(x, y, groups, ...)
        rects[panel.num == panel.number() & get(gID) == grps[1L] & star == '*',
              panel.text(k, yloc, labels='*', col=plot.cols[1L])]
        rects[panel.num == panel.number() & get(gID) == grps[2L] & star == '*',
              panel.text(k, yloc, labels='*', col=plot.cols[2L])]
      }
      p <- xyplot(norm ~ k | get(facet.by), groups=get(gID), data=rects,
                  ylim=ylimit, xlab=list(label='Degree (k)', cex=0.8), ylab=ylabel,
                  xleft=rects$xstart, xright=rects$xend, type=plot.type, panel=panelfun,
                  auto.key=list(space='bottom', title=gID, cex.title=1, columns=length(grps),
                                lines=TRUE, points=FALSE))
    }

  # 'ggplot2' plotting
  #-----------------------------------------------------------------------------
  } else {
    if (isTRUE(smooth)) {
      p <- ggplot2::ggplot(data=rects, ggplot2::aes(x=k, y=norm, group=get(sID))) +
        ggplot2::stat_smooth(ggplot2::aes(col=get(gID), group=get(gID)))
    } else {
      p <- ggplot2::ggplot(data=rects, ggplot2::aes(x=k, y=norm)) +
        ggplot2::geom_line(ggplot2::aes(col=get(gID)), size=1.25) +
        ggplot2::geom_rect(data=rects,
                  ggplot2::aes(x=NULL, y=NULL, xmin=xstart, xmax=xend, ymin=-Inf, ymax=Inf),
                  alpha=0.02, fill='lightpink')
    }
    p <- p +
      ggplot2::geom_hline(yintercept=1, size=0.5, lty=2) +
      ggplot2::geom_text(ggplot2::aes(y=yloc, col=get(gID), label=star), size=6, show.legend=FALSE) +
      ggplot2::facet_wrap(as.formula(paste('~', facet.by)), scales='free') +
      ggplot2::labs(x='Degree (k)', y=expression(phi[norm])) +
      ggplot2::theme(legend.position='bottom')
  }

  return(p)
}
