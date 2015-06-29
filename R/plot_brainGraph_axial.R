#' Draw axial slice of the MNI152 T1 image.
#'
#' This function draws an axial slice from the MNI152 T1 image, to plot the
#' nodes of a graph over it. It will optionally write to a filename for output.
#'
#' @param save Binary indicating whether or not a png file should be saved
#' @param fname The name of the file to be saved
#' @param z The z-coordinate of the slice to use (defaults to 46, the center)
#' @export
#'
#' @seealso \code{\link{plot_brainGraph_sagittal}}

plot_brainGraph_axial <- function(save=0, fname=NULL, z=46) {
    if (save==1) {
      png(filename=fname)
    } else {
      if (length(dev.list()) == 0) {
        dev.new()
      }
    }
    image(mni152, plot.type='single', z=z, zlim=c(3500, max(mni152[, , z])))
    par(new=T, mai=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
}
