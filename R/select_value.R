#' Retrieve the Nearest Dataset from a Ramanome Object
#'
#' This function extracts the most recently added dataset from a Ramanome object.
#' It is assumed that the dataset is stored in the `datasets` slot of the Ramanome object.
#'
#' @param object A Ramanome object that contains datasets in its `datasets` slot.
#' @return The most recently added dataset from the Ramanome object.
#' @export get.nearest.dataset

get.nearest.dataset <- function(object) {
  dataset <- tail(names(object@datasets), 1)
  dataset <- object@datasets[dataset][[1]]
  return(dataset)
}

#' Select a single Raman intensity value for a given wavenumber
#'
#' This function selects a single Raman intensity value for a given wavenumber from a Ramanome object.
#'
#' @param object A Ramanome object.
#' @param wave The wavenumber for which to retrieve the intensity value.
#' @export select.value
#' @return The Raman intensity value for the given wavenumber.
select.value <- function(object, wave) {
  loc <- which.min(abs(object@wavenumber - wave))
  dataset <- get.nearest.dataset(object)
  return(dataset[, loc])
}

#' Select a range of Raman intensity values for a given wavenumber range
#'
#' This function selects a range of Raman intensity values for a given wavenumber range from a Ramanome object.
#'
#' @param object A Ramanome object.
#' @param waves A vector with two values representing the lower and upper bounds of the wavenumber range.
#' @export select.band
#' @return The Raman intensity values within the given wavenumber range.
select.band <- function(object, waves) {
  locs <- object@wavenumber <= waves[2] & object@wavenumber >= waves[1]
  dataset <- get.nearest.dataset(object)
  return(rowSums(dataset[, locs]))
}

#' string. If two values are provided, they are combined into a range string with a
#' tilde (~) separator. If more than two values are provided, the function stops with
#' an error.
#'
#' @param waves A single numeric value or a vector of two numeric values representing
#' wavelengths or a wavelength range.
#' @return A character string representing the formatted wavelength name or range.
#' @export confirm.name

confirm.name <- function(waves) {
  if (length(waves) == 1) {
    name <- as.character(waves)
  } else if (length(waves) == 2) {
    name <- paste(waves[1], waves[2], sep = '~')
  } else {
    stop('Error! Please input a wavenumber or a Raman band!')
  }
  return(name)
}

#' Confirm Selection and Retrieve Spectral Data
#'
#' This function determines the type of selection based on the input and retrieves the
#' corresponding spectral data from a Ramanome object. If a single wavelength is
#' provided, it uses `select.value` to retrieve the data at that wavelength. If a
#' range of wavelengths is provided, it uses `select.band` to retrieve the sum of
#' the data within that range.
#'
#' @param object A Ramanome object containing spectral data and wavelength information.
#' @param waves A single numeric value representing a wavelength or a vector of two
#' numeric values representing a wavelength range.
#' @return The selected spectral data based on the input wavelength or range.
#' @export confirm.select
#' @seealso select.value for retrieving data at a single wavelength.
#' @seealso select.band for retrieving the sum of data within a wavelength range.
#'
confirm.select <- function(object, waves) {
  if (length(waves) == 1) {
    values <- select.value(object, waves)
  } else if (length(waves) == 2) {
    values <- select.band(object, waves)
  }else { print('Error! Please input a wavenumber or a Raman band!') }
  return(values)
}


#' Extract Intensity Values at Specified Wavelengths
#'
#' This function retrieves the intensity values at specified wavelengths or wavelength
#' ranges from a Ramanome object. It uses `confirm.select` to get the spectral data and
#' `confirm.name` to format the names of the selected wavelengths or ranges.
#'
#' @param object A Ramanome object containing spectral data and wavelength information.
#' @param wavenumber A single numeric value, a vector of numeric values, or a list of
#' numeric values representing wavelengths or wavelength ranges.
#' @return The updated Ramanome object with the extracted intensity values appended to
#' the `interested.bands` slot.
#' @export intensity
intensity <- function(object, wavenumber) {
  wavenumber <- as.list(wavenumber)
  a <- lapply(wavenumber, function(x) confirm.select(object, x))
  name <- lapply(wavenumber, confirm.name)
  object@interested.bands <- append(object@interested.bands, a)
  names(object@interested.bands) <- c(names(object@interested.bands), unlist(name))
  return(object)
}
