#' Plot vectors for 2 groups along a range of some value.
#'
#' This function takes a vector of values from each of two groups, and plots
#' them over their respective ranges (e.g. plot the average degree over a range
#' of correlation thresholds). The default xlabel here is 'Cost'.
#'
#' @param yvar.1 the first group's data
#' @param yvar.2 the second group's data
#' @param xlabel the label to place on the x-axis (optional)
#' @param ylabel the label to place on the y-axis (optional)
#' @param group1 Character string identifying first subject group
#' @param group2 Character string identifying second subject group
#' @export

plot.threshes <- function(yvar.1, yvar.2, xlabel='Cost', ylabel=NULL, group1, group2) {
  ymin <- min(yvar.1, yvar.2, na.rm=T)
  ymax <- max(yvar.1, yvar.2, na.rm=T)

  plot(threshes[1:length(yvar.1)], yvar.1, type='l', col='blue', ann=F,
       ylim=c(ymin, ymax))
  grid()
  lines(threshes[1:length(yvar.2)], yvar.2, col='red')
  legend('topright', c(group1, group2), cex=0.75, col=c('blue', 'red'),
         lty=c(1,1))
  title(xlab=xlabel, ylab=ylabel)
}
