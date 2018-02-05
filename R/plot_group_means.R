#' Plot group distributions of volumetric measures for a given brain region
#'
#' This function takes a "tidied" dataset of cortical volumetric measures
#' (thickness, volume, LGI, etc.) and plots a histogram or violin plot for 1 or
#' more groups, and of 1 or more brain regions.
#'
#' @param dat A data table of volumetric data; needs columns for 'Group',
#' 'region', and 'value'
#' @param regions A vector of character strings or integers of the brain
#' region(s) to plot; if integer, the region(s) is/are chosen from the input
#' data table based on the index
#' @param type A character string indicating the plot type; either 'histogram'
#' or 'violin'
#' @param all.vals A logical indicating whether or not to plot horizontal lines
#' for all observations (only valid for 'violin' plots) (default: TRUE)
#' @param modality A character string indicating the type of volumetric measure
#' ('thickness', 'volume', 'lgi', or 'area')
#' @export
#' @importFrom scales fullseq
#'
#' @return A ggplot object
#' @family Structural covariance network functions
#' @seealso \code{\link[ggplot2]{geom_histogram}, \link[ggplot2]{geom_vline}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_volumetric <- function(dat, regions, type=c('violin', 'histogram'),
                            all.vals=TRUE,
                            modality=c('thickness', 'volume', 'lgi', 'area')) {
  region <- value <- Group <- ..density.. <- avg <- group.mean <- bwidth <-
    breaks <- x <- width <- NULL
  if (!is.character(regions)) regions <- dat[, levels(region)][regions]
  stopifnot(all(regions %in% dat[, levels(region)]))

  modality <- match.arg(modality)
  ax.lab <- switch(modality,
                   thickness='Thickness (mm)',
                   volume=expression(paste('Volume (', mm^{3}, ')')),
                   lgi='Local Gyrification Index',
                   area=expression(paste('Surface area (', mm^{2}, ')')))

  subDT <- dat[region %in% regions]
  type <- match.arg(type)
  if (type == 'histogram') {
    # Allow for variable bin widths
    groups <- subDT[, levels(Group)]
    setkey(subDT, region, Group)
    breaksdt <- subDT[, list(breaks=pretty(range(value), n=nclass.FD(value))),
                      by=list(Group, region)]
    breaksdt[, bwidth := .SD[1:2, diff(breaks)], by=list(Group, region)]
    subDT[, bwidth := rep(breaksdt[, min(bwidth), by=region]$V1,
                          times=subDT[, .N, by=region]$N)]
    # A partial recreation of Hadley's ggplot2:::bin function
    create_bins <- function(x, binwidth) {
      breaks <- sort(scales::fullseq(range(x), binwidth, pad=TRUE))
      bins <- cut(x, breaks, include.lowest=TRUE, right=FALSE)
      left <- breaks[-length(breaks)]
      right <- breaks[-1]
      x <- (left + right) / 2
      width <- diff(breaks)

      out <- data.frame(count=as.numeric(tapply(rep(1, length(bins)), bins, sum,
                                                na.rm=TRUE)),
                        x=x,
                        width=width)
      out$count[is.na(out$count)] <- 0
      out$density <- out$count / out$width / sum(abs(out$count), na.rm=TRUE)
      out$ndensity <- out$density / max(abs(out$density), na.rm=TRUE)
      out$ncount <- out$count / max(abs(out$count), na.rm=TRUE)
      return(out)
    }
    my.df <- subDT[, create_bins(value, unique(bwidth)), by=list(Group, region)]

    meandt <- subDT[, list(avg=mean(value)), by=list(Group, region)]
    vol.plot <- ggplot(my.df) +
      geom_histogram(aes(x, y=density, width=width, fill=Group),
                     alpha=0.6, position='dodge', stat='identity') +
      geom_vline(data=meandt, aes(xintercept=avg, col=Group), lty=2, size=0.5) +
      geom_density(data=subDT, aes(x=value, col=Group), size=0.8) +
      facet_wrap(~ region, scales='free') +
      theme(legend.position=c(1, 1), legend.justification=c(1, 1),
            legend.background=element_rect(size=0.5)) +
      xlab(ax.lab)

  } else if (type == 'violin') {
    vol.plot <- ggplot(subDT, aes(x=Group, y=value, fill=Group)) +
      geom_violin(trim=FALSE) +
      facet_wrap(~ region, scales='free_y') +
      theme(legend.position='none') +
      ylab(ax.lab)
    if (isTRUE(all.vals)) {
      subDT[, group.mean := mean(value), by=list(Group, region)]
      vol.plot <- vol.plot +
        geom_segment(aes(x=as.numeric(Group)-0.1, xend=as.numeric(Group)+0.1,
                         y=value, yend=value), col='black') +
        geom_segment(aes(x=as.numeric(Group)-0.3, xend=as.numeric(Group)+0.3,
                         y=group.mean, yend=group.mean), col='black')
    } else {
      vol.plot <- vol.plot + geom_boxplot(width=0.1)
    }
  }
  return(vol.plot)
}

#' @inheritParams plot_volumetric
#' @export
#' @rdname plot_volumetric
plot_group_means <- function(dat, regions, type=c('violin', 'histogram'),
                             all.vals=TRUE,
                             modality=c('thickness', 'volume', 'lgi', 'area')) {
  .Deprecated('plot_volumetric')
  plot_volumetric(dat, regions, type, all.vals, modality)
}
