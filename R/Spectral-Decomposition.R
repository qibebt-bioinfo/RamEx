#' Perform Multivariate Curve Resolution Alternating Least Squares (MCR-ALS)
#'
#' This function applies the MCR-ALS algorithm to the dataset stored in a Ramanome object.
#' MCR-ALS is a chemometric technique used to resolve overlapping spectra into their
#' constituent pure components. The function returns the estimated components and
#' concentrations.
#'
#' @param object A Ramanome object containing the dataset.
#' @param ncomp The number of components to estimate.
#' @return A list containing the estimated components and concentrations.
#' @export Spectral.Decomposition.Mcrals
#' @importFrom mdatools mcrals
#' @importFrom mdatools constraint
#' @examples
#' # Create sample spectral data
#' wavenumbers <- seq(500, 3500, by = 10)
#' n_samples <- 100
#' n_waves <- length(wavenumbers)
#' 

Spectral.Decomposition.Mcrals<- function(object, ncomp) {
  mcrals <- mdatools::mcrals(get.nearest.dataset(object), ncomp,cont.constraints =list(
    constraint("norm", params = list(type = "sum"))
  ), spec.constraints=list(
    constraint("angle", params = list(weight = 0.05)),
    constraint("norm", params = list(type = "length"))
  ) )
  results <- list(components=mcrals$resspec, concentration=mcrals$rescont)
  return(results)
}




#' Perform Independent Component Analysis (ICA)
#'
#' This function applies the FastICA algorithm to the dataset retrieved from the Ramanome
#' object using `get.nearest.dataset`. It returns the result of the ICA decomposition.
#'
#' @param object A Ramanome object containing the dataset.
#' @param ncomp The number of components to retain in the ICA decomposition.
#' @return A list containing the results of the ICA decomposition.
#' @export Spectral.Decomposition.Ica
#' @importFrom fastICA fastICA
#' @importFrom ica icafast
#' @examples
#' # Create sample spectral data with mixed sources
#' wavenumbers <- seq(500, 3500, by = 10)
#' n_samples <- 100
#' n_waves <- length(wavenumbers)
#' 

Spectral.Decomposition.Ica <- function(object, ncomp) {
  ica.result <- icafast(get.nearest.dataset(object), ncomp)
  return(ica.result)
}



#' Perform Non-negative Matrix Factorization (NMF)
#'
#' This function applies Non-negative Matrix Factorization to the dataset retrieved from
#' a Ramanome object. NMF is a technique used to decompose a non-negative matrix into
#' two non-negative matrices, often used for data dimensionality reduction and feature
#' extraction.
#'
#' @param object A Ramanome object containing the dataset.
#' @return A list containing the results of the NMF decomposition, including the basis
#' matrices and the reconstructed matrix.
#' @export Spectral.Decomposition.Nmf
#' @importFrom NMF nmf
#' @examples
#' # Create sample non-negative spectral data
#' wavenumbers <- seq(500, 3500, by = 10)
#' n_samples <- 100
#' n_waves <- length(wavenumbers)
#' 
#' # Generate two non-negative basis spectra
#' basis1 <- dnorm(wavenumbers, mean = 1000, sd = 100)
#' basis2 <- dnorm(wavenumbers, mean = 2000, sd = 150)
#' 
#' # Create non-negative coefficients
#' coef1 <- abs(rnorm(n_samples, 2, 0.5))
#' coef2 <- abs(rnorm(n_samples, 1, 0.3))
#' 
#' # Create mixture spectra (ensuring non-negativity)
#' spectra <- matrix(0, nrow = n_samples, ncol = n_waves)
#' for(i in 1:n_samples) {
#'   spectra[i,] <- coef1[i] * basis1 + coef2[i] * basis2
#' }
#' 
#' # Create Ramanome object
#' raman_obj <- new("Ramanome",
#'   datasets = list(raw = spectra),
#'   wavenumber = wavenumbers
#' )
#' 
#' # Apply NMF
#' nmf_result <- Spectral.Decomposition.Nmf(raman_obj)

Spectral.Decomposition.Nmf <- function(object) {
  dataset <- get.nearest.dataset(object)
  res <- nmf(dataset, 2, method="ns", seed=123456)
  return(res)
}
