#' Plot normalized rich club coefficients against degree threshold
#'
#' Returns a \code{\link[ggplot2]{ggplot}} object of a line plot of the
#' normalized rich club coefficient. Optionally will include a shaded region
#' demarcating the \code{\link{rich_core}} cutoff (if you supply a list of graph
#' objects to the \code{g} argument).
#'
#' @param rich.dt A \code{data.table} with rich-club coefficients
#' @param facet.by A character string indicating whether the variable of
#' interest is "density" or "threshold" (e.g. with DTI data)
#' @param densities A numeric vector of the densities to plot
#' @param alpha The significance level (default: 0.05)
#' @param fdr A logical, indicating whether or not to use the FDR-adjusted
#'   p-value for determining significance (default: TRUE)
#' @param g A list (of lists) of \code{igraph} graph objects; required if you
#'   want to plot a shaded region demarcating the \code{\link{rich_core}}
#' @param smooth Logical indicating whether or not to use
#'   \code{\link[ggplot2]{stat_smooth}} when data from multiple subjects (per
#'   group) are present (default: \code{TRUE}). Ignored for group-level data.
#' @export
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @family Rich-club functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot_rich_norm(rich.dt, facet.by='density', densities[N:(N+1)], g=g)
#' }

plot_rich_norm <- function(rich.dt, facet.by=c('density', 'threshold'),
                           densities, alpha=0.05, fdr=TRUE, g=NULL, smooth=TRUE) {
  p.fdr <- yloc <- Group <- norm <- xstart <- xend <- Study.ID <- NULL

  facet.by <- match.arg(facet.by)
  subDT <- rich.dt[as.factor(get(facet.by)) %in% as.factor(densities)]
  if (isTRUE(fdr)) {
    subDT[, star := ifelse(p.fdr < alpha, '*', '')]
  } else {
    subDT[, star := ifelse(p < alpha, '*', '')]
  }
  subDT[, yloc := min(norm, na.rm=TRUE) - 0.05 * diff(range(norm, na.rm=TRUE)), by=facet.by]
  subDT[, Group := as.factor(Group)]
  groups <- subDT[, levels(Group)]
  if (length(groups) > 1) {
    for (i in 2:length(groups)) {
      subDT[Group == groups[i], yloc := yloc - i * 0.05 * diff(range(norm, na.rm=TRUE)), by=facet.by]
    }
  }
  setkeyv(subDT, c(facet.by, 'Group'))

  rects <- data.table(density=subDT[, unique(get(facet.by))], xstart=0, xend=0,
                      Group=rep(subDT[, unique(Group)],
                                each=length(densities)))
  if (!is.null(g)) {
    densities.g <- round(sapply(g[[1]], graph_attr, 'density'), 2)
    densities.g <- which(densities.g %in% round(densities, 2))
    g <- lapply(g, `[`, densities.g)
    k <- lapply(g, sapply, function(x) rich_core(x)$k.r)

    if (length(densities) == 1) {
      max.k <- apply(as.matrix(sapply(k, function(x) x), nrow=1), 2, max)
    } else {
      max.k <- apply(sapply(k, function(x) x), 1, max)
    }
    rects[, xstart := max.k]
    rects[, xend := subDT[, max(k), by=density]$V1]
  }
  setnames(rects, 'density', facet.by)
  setkeyv(rects, key(subDT))

  if ('Study.ID' %in% names(subDT)) {
    rects <- subDT[rects]
    p <- ggplot(data=rects, aes(x=k, y=norm, group=Study.ID))
  } else {
    p <- ggplot(data=subDT[rects], aes(x=k, y=norm))
    smooth <- FALSE
  }
  p <- p +
    geom_hline(yintercept=1, size=0.5, lty=2) +
    geom_text(aes(y=yloc, col=Group, label=star), size=6, show.legend=F)

  if (isTRUE(smooth)) {
    p <- p + stat_smooth(aes(col=Group, group=Group))
  } else {
    p <- p +
      geom_line(aes(col=Group)) +
      geom_rect(data=rects,
                aes(x=NULL, y=NULL, xmin=xstart, xmax=xend, ymin=-Inf, ymax=Inf),
                alpha=0.08, fill='red')
  }
  p <- p +
    facet_wrap(as.formula(paste('~', facet.by)), scales='free') +
    labs(x='Degree (k)', y=expression(phi[norm])) +
    theme(legend.position='bottom')

  return(p)
}
