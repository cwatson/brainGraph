#' Plot vectors for 2 groups along a range of some value, +/- std error.
#'
#' This function takes a vector of values from each of two groups, and plots
#' them over their respective ranges (e.g. plot the average degree over a range
#' of correlation thresholds). The default xlabel here is 'Cost'. It also plots
#' plus/minus the standard error of the mean and fills in between.
#'
#' @param y1 The first group's data
#' @param y2 The second group's data
#' @param xlabel The label to place on the x-axis (optional)
#' @param ylabel The label to place on the y-axis (optional)
#' @param group1 Character string identifying first subject group
#' @param group2 Character string identifying second subject group
#' @param y1se Standard error of the mean for group 1
#' @param y2se Standard error of the mean for group 2
#' @export

plot.se <- function(y1, y2, xlabel='Cost', ylabel=NULL, group1, group2,
                          y1se, y2se) {
  ymin <- min(y1 - y1se, y2 - y2se)
  ymax <- max(y1 + y1se, y2 + y2se)

  plot(threshes, y1, type='n', ylim=c(ymin, ymax), xlab=xlabel, ylab=ylabel)
  lines(threshes, y1, col='blue')
  polygon(c(threshes, rev(threshes)), c(y1 - y1se, rev(y1 + y1se)),
          col=rgb(0, 0, 1, 0.5), border=NA)
  lines(threshes, y2, col='red')
  polygon(c(threshes, rev(threshes)), c(y2 - y2se, rev(y2 + y2se)),
          col=rgb(1, 0, 0, 0.5), border=NA)
  grid()

  legend('topright', c(group1, group2), cex=0.75, fill=c('blue', 'red'))
}
