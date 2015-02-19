#' Plot group histograms of cortical thickness for a given brain region
#'
#' This function takes a "tidied" dataset of cortical thickness and plots a
#' histogram for 2 groups, including a dashed vertical line indicating the group
#' means.
#'
#' @param dat A data table of cortical thicknesses; needs columns for 'Group',
#' 'region', and 'thickness'
#' @param cur.region A character string or integer of the brain region to plot;
#' if integer, the region is chosen from the data table based on the index
#' @export
#'
#' @return A ggplot object for plotting the histograms
#' @seealso \code{\link[ggplot2]{geom_histogram}, \link[ggplot2]{geom_vline}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

hist_thickness <- function(dat, cur.region) {
  region <- thickness <- Group <- ..density.. <- avg <- NULL
  if(!is.character(cur.region)) {
    cur.region <- levels(dat$region)[cur.region]
  } else {
    if(!cur.region %in% levels(dat$region)) {
      stop(sprintf('%s is not an acceptable region!', cur.region))
    }
  }

  meandf <- dat[region == cur.region, list(avg=mean(thickness)), by=Group]
  ggplot(dat[region == cur.region], aes(x=thickness, fill=Group)) +
    geom_histogram(binwidth=0.05, alpha=0.4, position='dodge',
                   aes(y=..density..)) +
    geom_vline(data=meandf, aes(xintercept=avg, col=Group), lty=2, size=0.5) +
    ggtitle(cur.region)
}
