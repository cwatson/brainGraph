#' Plot group distributions of cortical thickness for a given brain region
#'
#' This function takes a "tidied" dataset of cortical thickness and plots a
#' histogram or violin plot for 1 or more groups.
#'
#' @param dat A data table of cortical thicknesses; needs columns for 'Group',
#' 'region', and 'thickness'
#' @param regions A vector of character strings or integers of the brain
#' region(s) to plot; if integer, the region(s) is/are chosen from the atlas
#' data table based on the index
#' @param type A character string indicating the plot type; either 'histogram'
#' or 'violin'
#' @param all.vals A logical indicating whether or not to plot horizontal lines
#' for all observations (only valid for 'violin' plots) (default: TRUE)
#' @export
#'
#' @return A ggplot object
#' @seealso \code{\link[ggplot2]{geom_histogram}, \link[ggplot2]{geom_vline}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_thickness <- function(dat, regions, type=c('violin', 'histogram'),
                           all.vals=TRUE) {
  region <- thickness <- Group <- ..density.. <- avg <- mean.thickness <- NULL
  if(!is.character(regions)) {
    regions <- dat[, unique(region)][regions]
  } else {
    if(!all(regions %in% levels(dat$region))) {
      stop(paste('All region names must be valid!'))
    }
  }
  subDT <- dat[region %in% regions]

  type <- match.arg(type)
  if (type == 'histogram') {
    meandt <- subDT[, .(avg=mean(thickness)), by=Group]
    thick.plot <- ggplot(subDT, aes(x=thickness)) +
      geom_histogram(binwidth=0.05, alpha=0.4, position='dodge',
                     aes(y=..density.., fill=Group)) +
      geom_vline(data=meandt, aes(xintercept=avg, col=Group), lty=2, size=0.5) +
      geom_density(aes(col=Group), size=0.8) +
      facet_wrap(~ region, scales='free') +
      theme(legend.position=c(1, 1), legend.justification=c(1, 1),
            legend.background=element_rect(size=0.5)) +
      xlab('Thickness (mm)')
  } else if (type == 'violin') {
    thick.plot <- ggplot(subDT,
                         aes(x=Group, y=thickness, fill=Group)) +
      geom_violin(trim=FALSE) +
      facet_wrap(~ region, scales='free_y') +
      theme(legend.position='none') +
      ylab('Thickness (mm)')
    if (isTRUE(all.vals)) {
      thick.plot <- thick.plot +
        geom_segment(aes(x=as.numeric(Group)-0.1, xend=as.numeric(Group)+0.1,
                         y=thickness, yend=thickness), col='black') +
        geom_segment(aes(x=as.numeric(Group)-0.3, xend=as.numeric(Group)+0.3,
                         y=mean.thickness, yend=mean.thickness), col='black')
    } else {
      thick.plot <- thick.plot + geom_boxplot(width=0.1)
    }
  }
  return(thick.plot)
}
