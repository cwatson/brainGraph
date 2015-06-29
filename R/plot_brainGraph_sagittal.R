#' Draw sagittal slice of the MNI152 T1 image.
#'
#' This function draws an sagittal slice from the MNI152 T1 image, to plot the
#' nodes of a graph over it. It will optionally write to a filename for output.
#'
#' @param save Binary indicating whether or not a png file should be saved
#' @param fname Character string; the name of the file to be saved
#' @param z The x-coordinate of the slice to use (defaults to 50, the center)
#' @param hemi Character string indicating the hemisphere (either 'left' or
#' 'right')
#' @export
#'
#' @seealso \code{\link{plot_brainGraph_axial}}

plot_brainGraph_sagittal <- function(save=0, fname=NULL, z=50,
                                     hemi=c('left', 'right')) {
    hemi <- match.arg(hemi)
    if (hemi == 'right') {
      x <- mni152
    } else if (hemi == 'left') {
      tmp <- mni152@.Data
      x <- nifti(tmp[rev(seq_len(nrow(tmp))), rev(seq_len(ncol(tmp))), ])
    } else {
      stop('Argument "hemi" must be "left" or "right".')
    }

    if (save==1) {
      png(filename=fname)
    } else {
      if (length(dev.list()) == 0) {
        dev.new()
      }
    }
    image(x, plot.type='single', plane='sagittal',
        z=z, zlim=c(3500, max(mni152[z, , ])))
    par(new=T, mai=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
}
