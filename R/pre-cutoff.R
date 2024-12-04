#' Pre-process Raman data by cutting the wavenumber range
#'
#' This function cuts the wavenumber range of Raman data and updates the Raman object.
#'
#' @param object A Ramanome object.
#' @param from The lower bound of the wavenumber range to cut.
#' @param to The upper bound of the wavenumber range to cut.
#'
#' @return The updated Ramanome object with the wavenumber range cut and the corresponding dataset.
#'
#' @export pre.cutoff
pre.cutoff <- function(object, from, to) {
  waves <- object@wavenumber
  pred.data <- get.nearest.dataset(object)
  bands <- waves > from & waves < to
  pred.data <- pred.data[, bands]
  object@wavenumber <- waves[bands]
  object@datasets$cut.data <- pred.data
  return(object)
}
