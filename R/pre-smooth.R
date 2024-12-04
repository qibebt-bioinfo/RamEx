#' Smooth Spectral Data Using Savitzky-Golay Filter
#'
#' This function applies the Savitzky-Golay filter to the spectral data in a Ramanome
#' object to reduce noise. The filter is applied to the central portion of the data,
#' preserving the first and last few columns as specified by the `w` parameter.
#'
#' @param object A Ramanome object containing the spectral data.
#' @param m The order of the polynomial used in the Savitzky-Golay filter. Defaults to 0.
#' @param p The number of points to use in the Savitzky-Golay filter. Defaults to 5.
#' @param w The window size used in the Savitzky-Golay filter. Defaults to 11.
#' @param delta.wav The change in wavelength between each data point. Defaults to 2.
#' @return The updated Ramanome object with the smoothed spectral data.
#' @export pre.smooth.sg
#' @importFrom prospectr savitzkyGolay

pre.smooth.sg <- function(object, m = 0, p = 5, w = 11, delta.wav = 2) {
  pred.data <- get.nearest.dataset(object)
  pred.data <- cbind(
    pred.data[, 1:((w - 1) / 2)],
    prospectr::savitzkyGolay(pred.data, m = m, p = p, w = w, delta.wav = delta.wav),
    pred.data[, (ncol(pred.data) - (w - 1) / 2 + 1):ncol(pred.data)]
  )
  object@datasets$smooth.data <- pred.data
  return(object)
}

#' Apply Standard Normal Variate (SNV) Transformation
#'
#' This function applies the Standard Normal Variate (SNV) transformation to the raw
#' spectral data in a Ramanome object. The SNV transformation is a type of preprocessing
#' technique used to normalize the data and reduce the impact of multiplicative scatter.
#' The function also reports the number of spectra and the wavelength range after
#' the transformation.
#'
#' @param object A Ramanome object containing the raw spectral data.
#' @return The updated Ramanome object with the SNV-transformed data.
#' @export pre.smooth.snv
#' @importFrom prospectr standardNormalVariate

pre.smooth.snv <- function(object) {
  #pred.data <- get.nearest.dataset(object)

  object@datasets$raw.data <- prospectr::standardNormalVariate(as.data.frame(object@datasets$raw.data))

  writeLines(paste("after snv there is ", length(object), " spec data", sep = ""))
  min_wave <- min(object@wavenumber)
  max_wave <- max(object@wavenumber)
  writeLines(paste("wave range is ", min_wave, "-", max_wave, sep = ""))
  return(object)
}
