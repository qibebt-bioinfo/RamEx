#' Pre-process Raman data by ch normalization method
#'
#' This function applies normalization to the Raman data based on the ch method and updates the Raman object.
#'
#' @param object A Ramanome object.
#'
#' @return The updated Ramanome object with the ch normalized Raman data.
#'
#' @export pre.normalize.ch
pre.normalize.ch <- function(object) {
  dataset <- get.nearest.dataset(object)
  range <- object@wavenumber > 2850 & object@wavenumber < 3000
  value <- apply(dataset[, range], 1, max)
  object@datasets$normalized.data <- dataset / value
  return(object)
}

#' Pre-process Raman data by max normalization
#'
#' This function applies normalization to the Raman data based on the max method and updates the Raman object.
#'
#' @param object A Ramanome object.
#'
#' @return The updated Ramanome object with the max normalized Raman data.
#'
#' @export pre.normalize.max
pre.normalize.max <- function(object) {
  dataset <- get.nearest.dataset(object)
  value <- apply(dataset, 1, max)
  object@datasets$normalized.data <- dataset / value
  return(object)
}


#' Pre-process Raman data by specific normalization
#'
#' This function applies normalization to the Raman data based on the specific method and updates the Raman object.
#'
#' @param object A Ramanome object.
#' @param wave The specific wavenumber to use for the "specific" normalization method.
#'
#' @return The updated Ramanome object with the specific normalized Raman data.
#'
#' @export pre.normalize.specific
pre.normalize.specific <- function(object, wave) {
  dataset <- get.nearest.dataset(object)
  if (is.null(wave)) {
    stop('Error! Please input the interested wavenumber!')
  }
  loc <- which.min(abs(object@wavenumber - wave))
  value <- dataset[, loc]
  object@datasets$normalized.data <- dataset / value
  return(object)
}


#' Pre-process Raman data by area normalization method
#'
#' This function applies normalization to the Raman data based on the area method and updates the Raman object.
#'
#' @param object A Ramanome object.
#'
#' @return The updated Ramanome object with the area normalized Raman data.
#'
#' @export pre.normalize.area
pre.normalize.area <- function(object) {
  dataset <- get.nearest.dataset(object)
  value <- rowSums(dataset)
  object@datasets$normalized.data <- dataset / value
  return(object)
}
