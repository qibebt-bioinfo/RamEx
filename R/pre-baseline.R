#' Preprocess baseline
#'
#' This function performs baseline pre-processing on the input object.
#'
#' @param object The input object
#' @param order The order of the polynomial used in the baseline calculation (default: 1)
#'
#' @return The modified object with preprocessed baseline data
#'
#' @importFrom hyperSpec spc.fit.poly.below
#' @importFrom hyperSpec spc.fit.poly
#' @export pre.baseline.polyfit
pre.baseline.polyfit <- pre.baseline <- function(object, order = 1) {
  wavenumber <- object@wavenumber
  spc2hs <- new(
    "hyperSpec",
    spc = get.nearest.dataset(object),
    wavelength = wavenumber
  )
  
  data_hyperSpec <- spc2hs
  wave_max <- max(data_hyperSpec@wavelength)
  wave_min <- min(data_hyperSpec@wavelength)
  hyperspec <- data_hyperSpec
  if (wave_max >= 3050 & wave_min <= 600) {
    hyperspec_1 <- hyperspec[, , floor(wave_min) ~ 1790] -
      spc.fit.poly.below(
        hyperspec[, , floor(wave_min) ~ 1790],
        hyperspec[, , floor(wave_min) ~ 1790],
        poly.order = order
      )
    hyperspec_2 <- hyperspec[, , 1790 ~ floor(wave_max)] -
      spc.fit.poly(
        hyperspec[, , c(1790 ~ 2065, 2300 ~ 2633, 2783, floor(wave_max))],
        hyperspec[, , 1790 ~ floor(wave_max)],
        poly.order = 6
      )
    hyperspec_baseline <- cbind(hyperspec_1, hyperspec_2)
    print(paste0("The Ramanome contains ", length(hyperspec_baseline), " spectra"))
  }
  else if (wave_max > 1750 & wave_min < 600) {
    hyperspec_baseline <- hyperspec[, , 550 ~ 1750] -
      spc.fit.poly.below(
        hyperspec[, , 550 ~ 1750],
        hyperspec[, , 550 ~ 1750],
        poly.order = order
      )
    print(paste0("The Ramanome contains ", length(hyperspec_baseline), "spectra"))
  }
  else if (wave_max > 3050 &
           wave_min < 1800 &
           wave_min >
           600) {
    hyperspec_baseline <- hyperspec[, , 1800 ~ 3050] -
      spc.fit.poly(
        hyperspec[, , c(1800 ~ 2065, 2300 ~ 2633, 2783, 3050)],
        hyperspec[, , 1800 ~ 3050],
        poly.order = order
      )
    print(paste0("The Ramanome contains ", length(hyperspec_baseline), "spectra"))
  }
  else if (wave_max < 3050 & wave_min > 600) {
    hyperspec_baseline <- hyperspec[, , floor(wave_min) ~ floor(wave_max)] -
      spc.fit.poly(
        hyperspec[, , c(floor(wave_min) ~ floor(wave_min) + 10, 1800~floor(wave_max))],
        hyperspec[, ,floor(wave_min) ~ floor(wave_max) ],
        poly.order = order
      )}
  else if (wave_max < 1750 & wave_min > 600) {
    print("the spc is too small to baseline")
  }
  else {
    print("something is wrong in your spc data")
  }
  
  data_hyperSpec_baseline <- hyperspec_baseline

  data <- data_hyperSpec_baseline$spc[, !duplicated(colnames(data_hyperSpec_baseline$spc))]
  object@datasets$baseline.data <- data
  object@wavenumber <- as.numeric(colnames(data))
  
  return(object)
}


#' Remove baseline by bubble
#'
#' This function performs baseline pre-processing on the input object.
#'
#' @param object Ramanome object
#'
#' @return The baseline corrected spectra with bubble method
#' @export pre.baseline.bubble
pre.baseline.bubble <- function(object){
  band_1 <- which(object@wavenumber < 1800)
  band_2 <- which(object@wavenumber > 1800 & object@wavenumber < 2700)
  band_3 <- which(object@wavenumber > 2700)
  data_matrix <- get.nearest.dataset(object)
  data_matrix[,band_1] <- t(apply(data_matrix[,band_1], 1, function(x)bubblefill(x)$raman))
  data_matrix[,band_2] <- polyfit(data_matrix[,band_2], degree = 1, tol = 0.001, rep = 100)$corrected
  data_matrix[,band_3] <- polyfit(data_matrix[,band_3], degree = 1, tol = 0.001, rep = 100)$corrected
  data_matrix[, band_2] <- data_matrix[, band_2] - (data_matrix[,band_2[1]]-data_matrix[,band_2[1]-1])
  data_matrix[, band_3] <- data_matrix[, band_3] - (data_matrix[,band_3[1]]-data_matrix[,band_3[1]-1])
  object@datasets$baseline.data <- data_matrix
  return(object)
}
