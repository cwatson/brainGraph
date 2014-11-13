#' Plot vectors for 2 groups along a range of some value, +/- std error.
#'
#' This function takes a vector of values from each of two groups, and plots
#' them over their respective ranges (e.g. plot the average degree over a range
#' of network densities). The default xlabel here is 'Density'. It also plots
#' plus/minus the standard error of the mean and fills in between.
#'
#' @param densities The densities to be plotted along the horizontal axis
#' @param y1 The first group's data
#' @param y2 The second group's data
#' @param xlabel The label to place on the x-axis (optional)
#' @param ylabel The label to place on the y-axis (optional)
#' @param group1 Character string identifying first subject group
#' @param group2 Character string identifying second subject group
#' @param y1se Standard error of the mean for group 1 (optional)
#' @param y2se Standard error of the mean for group 2 (optional)
#' @export

plot.se <- function(densities, y1, xlabel='Density', y2=NULL, ylabel=NULL,
                    group1='Control', group2=NULL, y1se=NULL, y2se=NULL) {
  if (length(y1se) == 0) {
    y1se <- 0
  }
  if (length(y2se) == 0) {
    y2se <- 0
  }
  if (length(y2) == 0) {
    ymin <- min(y1 - y1se)
    ymax <- max(y1 + y1se)
  } else {
    ymin <- min(y1 - y1se, y2 - y2se)
    ymax <- max(y1 + y1se, y2 + y2se)
  }

  plot(densities, y1, type='n', ylim=c(ymin, ymax), xlab=xlabel, ylab=ylabel)
  lines(densities, y1, col='blue')
  polygon(c(densities, rev(densities)), c(y1 - y1se, rev(y1 + y1se)),
          col=rgb(0, 0, 1, 0.5), border=NA)

  if (length(y2) != 0) {
    lines(densities, y2, col='red')
    polygon(c(densities, rev(densities)), c(y2 - y2se, rev(y2 + y2se)),
            col=rgb(1, 0, 0, 0.5), border=NA)
  }
  grid()

  legend('topright', c(group1, group2), cex=0.75, fill=c('blue', 'red'))
}
