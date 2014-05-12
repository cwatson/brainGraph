#' Draw axial slice of the MNI152 T1 image.
#'
#' This function draws an axial slice from the MNI152 T1 image, to plot the
#' nodes of a graph over it. It will optionally write to a filename for output.
#'
#' @param flag binary indicating whether or not a png file should be saved
#' @param fname the name of the file to be saved
#' @param z The z-coordinate of the slice to use (defaults to 46, the center)
#' @export

plot.over.brain.axial <- function(flag=0, fname=NULL, z=46) {
    if (flag==1) {
      png(file=fname)
    } else {
      if (length(dev.list()) == 0) {
        dev.new()
      }
    }
    #plot.new()
    image(mni152, plot.type='single', z=z, zlim=c(3500, max(mni152[, , z])))
    par(new=T, mai=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
}
