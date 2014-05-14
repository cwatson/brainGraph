#' Draw an axial outline of a brain.
#'
#' This function draws an outline of a brain (axial view), to plot the nodes of
#' a graph over it. It will optionally write to a filename for output. It
#' requires that the \code{\link{png}} package be loaded. This specific png is
#' 412x480 in resolution; the outline goes from ~x=50 -> 350, y=60 -> 420.
#'
#' @param flag binary indicating whether or not a png file should be saved
#' @param fname the name of the file to be saved
#' @export

plot.over.brain <- function(flag=0, fname=NULL) {
    if (flag==1) {
      png(filename=fname)
    } else {
      if (length(dev.list() == 0)) {
        dev.new()
      }
    }
    plot.new()
    rasterImage(r, -0.12, -0.21, 1.12, 1.22)
    par(new=T)
}
