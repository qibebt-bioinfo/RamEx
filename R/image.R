#' Generate Raman imaging
#'
#' This function generates a Raman imaging plot for a specified peak in a Ramanome object.
#'
#' @param object The Ramanome object
#' @param peak The specified peak for Raman imaging
#'
#' @importFrom grDevices jpeg
#' @importFrom graphics image
#' @importFrom rlist list.map
#' @importFrom grDevices dev.off
#' @export image.peak

image.peak <- function(object, peak) {
  data <- data.frame(
    x = as.numeric(object@meta.data$x),
    y = as.numeric(object@meta.data$y),
    value = object@interested.bands[[as.character(peak)]] * 100
  )
  data.matrix <- lapply(
    unique(data$y), function(loc) {
      rr.data <- data[data$y == loc,]
      return(rr.data$value[order(rr.data$x)])
    }
  )
  data.matrix <- do.call(rbind, rlist::list.map(data.matrix, .))
  grDevices::jpeg(
    filename = paste0('Raman_Imaging_', peak, ' .jpeg'),
    width = length(unique(data$x)) * 30,
    height = length(unique(data$x)) * 30,
    quality = 150,
    res = 200
  )
  graphics::image(data.matrix, axes = F)
  grDevices::dev.off()
}
