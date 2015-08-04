#' Plot group distributions of cortical thickness for a given brain region
#'
#' This function takes a "tidied" dataset of cortical thickness and plots a
#' histogram or violin plot for 1 or more groups.
#'
#' @param dat A data table of cortical thicknesses; needs columns for 'Group',
#' 'region', and 'thickness'
#' @param cur.region A character string or integer of the brain region to plot;
#' if integer, the region is chosen from the data table based on the index
#' @param type A character string indicating the plot type; either 'histogram'
#' or 'violin'
#' @export
#'
#' @return A ggplot object
#' @seealso \code{\link[ggplot2]{geom_histogram}, \link[ggplot2]{geom_vline}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

plot_thickness <- function(dat, cur.region, type=c('histogram', 'violin')) {
  region <- thickness <- Group <- ..density.. <- avg <- NULL
  if(!is.character(cur.region)) {
    cur.region <- levels(dat$region)[cur.region]
  } else {
    if(!cur.region %in% levels(dat$region)) {
      stop(sprintf('%s is not an acceptable region!', cur.region))
    }
  }

  type <- match.arg(type)
  if (type == 'histogram') {
    meandt <- dat[region == cur.region, list(avg=mean(thickness)), by=Group]
    thick.plot <- ggplot(dat[region == cur.region], aes(x=thickness)) +
      geom_histogram(binwidth=0.05, alpha=0.4, position='dodge',
                     aes(y=..density.., fill=Group)) +
      geom_vline(data=meandt, aes(xintercept=avg, col=Group), lty=2, size=0.5) +
      geom_density(aes(col=Group), size=0.8)
  } else if (type == 'violin') {
    thick.plot <- ggplot(dat[region == cur.region],
                         aes(x=Group, y=thickness, fill=Group)) +
      geom_violin() +
      geom_boxplot(width=0.1)
  }
  thick.plot <- thick.plot + ggtitle(cur.region) + xlab('Cortical thickness (mm)')
}
