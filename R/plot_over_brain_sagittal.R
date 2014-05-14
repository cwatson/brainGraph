#' Draw sagittal slice of the MNI152 T1 image.
#'
#' This function draws an sagittal slice from the MNI152 T1 image, to plot the
#' nodes of a graph over it. It will optionally write to a filename for output.
#'
#' @param flag Binary indicating whether or not a png file should be saved
#' @param fname Character string; the name of the file to be saved
#' @param z The x-coordinate of the slice to use (defaults to 50, the center)
#' @param hemi The hemisphere (either 'left' or 'right')
#' @export

plot.over.brain.sagittal <- function(flag=0, fname=NULL, z=50, hemi) {
    if (hemi == 'right') {
      x <- mni152
    } else if (hemi == 'left') {
      tmp <- mni152@.Data
      x <- nifti(tmp[rev(seq_len(dim(tmp)[1])), rev(seq_len(dim(tmp)[2])), ])
    } else {
      stop('Argument "hemi" must be "left" or "right".')
    }

    if (flag==1) {
      png(filename=fname)
    } else {
      dev.new()
    }
    plot.new()
    image(x, plot.type='single', plane='sagittal',
        z=z, zlim=c(3500, max(mni152[z, , ])))
    par(new=T, mai=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
}
