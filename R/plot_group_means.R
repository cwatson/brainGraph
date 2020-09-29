#' Plot group distributions of volumetric measures for a given brain region
#'
#' This function takes a \dQuote{tidied} dataset of cortical volumetric measures
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
#'
#' @return A \code{trellis} or \code{ggplot} object
#' @family Structural covariance network functions
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
    if (!all(vapply(c('scales', 'ggplot2'), requireNamespace, logical(1L), quietly=TRUE))) {
      stop('Please install the "scales" and "ggplot2" packages to plot with histograms.')
    }
    # Allow for variable bin widths
    groups <- subDT[, levels(Group)]
    setkey(subDT, region, Group)
    breaksdt <- subDT[, list(breaks=pretty(range(value), n=nclass.FD(value))),
                      by=list(Group, region)]
    breaksdt[, bwidth := .SD[1L:2L, diff(breaks)], by=list(Group, region)]
    subDT[, bwidth := rep(breaksdt[, min(bwidth), by=region]$V1,
                          times=subDT[, .N, by=region]$N)]
    # A partial recreation of Hadley's ggplot2:::bin function
    create_bins <- function(x, binwidth) {
      breaks <- sort(scales::fullseq(range(x), binwidth, pad=TRUE))
      bins <- cut(x, breaks, include.lowest=TRUE, right=FALSE)
      left <- breaks[-length(breaks)]
      right <- breaks[-1L]
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
    vol.plot <- ggplot2::ggplot(my.df) +
      ggplot2::geom_histogram(ggplot2::aes(x, y=density, width=width, fill=Group),
                     alpha=0.6, position='dodge', stat='identity') +
      ggplot2::geom_vline(data=meandt, ggplot2::aes(xintercept=avg, col=Group), lty=2, size=0.5) +
      ggplot2::geom_density(data=subDT, ggplot2::aes(x=value, col=Group), size=0.8) +
      ggplot2::facet_wrap(~ region, scales='free') +
      ggplot2::theme(legend.position=c(1, 1), legend.justification=c(1, 1),
            legend.background=ggplot2::element_rect(size=0.5)) +
      ggplot2::xlab(ax.lab)

  } else if (type == 'violin') {
    # 'base' plotting
    if (!requireNamespace('ggplot2', quietly=TRUE)) {
      gID <- getOption('bg.group')
      vol.plot <- bwplot(value ~ get(gID) | region, data=subDT, panel=panel.violin, xlab=gID, ylab=ax.lab)

    # 'ggplot2' plotting
    } else {
      vol.plot <- ggplot2::ggplot(subDT, ggplot2::aes(x=Group, y=value, fill=Group)) +
        ggplot2::geom_violin(trim=FALSE) +
        ggplot2::facet_wrap(~ region, scales='free_y') +
        ggplot2::theme(legend.position='none') +
        ggplot2::ylab(ax.lab)
      if (isTRUE(all.vals)) {
        subDT[, group.mean := mean(value), by=list(Group, region)]
        vol.plot <- vol.plot +
          ggplot2::geom_segment(ggplot2::aes(x=as.numeric(Group)-0.1, xend=as.numeric(Group)+0.1,
                                             y=value, yend=value), col='black') +
          ggplot2::geom_segment(ggplot2::aes(x=as.numeric(Group)-0.3, xend=as.numeric(Group)+0.3,
                                             y=group.mean, yend=group.mean), col='black')
      } else {
        vol.plot <- vol.plot + ggplot2::geom_boxplot(width=0.1)
      }
    }
  }
  return(vol.plot)
}
