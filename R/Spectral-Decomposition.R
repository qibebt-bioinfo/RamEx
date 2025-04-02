#' Multivariate Curve Resolution Alternating Least Squares (MCR-ALS)
#'
#' A chemometric method that decomposes a overlapping spectra matrix into concentration and pure spectra profiles iteratively under constraints like non-negativity.
#' This function applies the MCR-ALS algorithm to the dataset stored in a Ramanome object.
#' @param object A Ramanome object.
#' @param n_comp The number of components to estimate.
#' @param seed The seed for the random number generator.
#' @return A list containing: 
#' \item{components}{The estimated components (n_wavelength x n_comp).} 
#' \item{concentration}{The estimated concentrations of each cell (n_cell x n_comp).}
#' while the reconstructed matrix is the product of the concentration and component matrices.
#' @export Spectral.Decomposition.Mcrals
#' @importFrom mdatools mcrals
#' @importFrom mdatools constraint
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' decom_mcrals <- Spectral.Decomposition.Mcrals(data_processed, 2)

Spectral.Decomposition.Mcrals<- function(object, n_comp=2, seed=42) {
  set.seed(seed)
  mcrals <- mdatools::mcrals(get.nearest.dataset(object), n_comp,cont.constraints =list(
    constraint("norm", params = list(type = "sum"))
  ), spec.constraints=list(
    constraint("angle", params = list(weight = 0.05)),
    constraint("norm", params = list(type = "length"))
  ) )
  results <- list(components=mcrals$resspec, concentration=mcrals$rescont)
  return(results)
}




#' Independent Component Analysis (ICA)
#'
#' A statistical technique that separates a multivariate signal into additive, statistically independent components, assuming non-Gaussian sources.
#'
#' @param object A Ramanome object.
#' @param n_comp The number of components to retain in the ICA decomposition.
#' @param seed The seed for the random number generator.
#' @return A list containing the results of the ICA decomposition.
#' @export Spectral.Decomposition.Ica
#' @importFrom ica icafast
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' decom_ica <- Spectral.Decomposition.Ica(data_processed, 2)

Spectral.Decomposition.Ica <- function(object, n_comp=2, seed=42) {
  set.seed(seed)
  ica.result <- icafast(get.nearest.dataset(object), n_comp)
  return(ica.result)
}



#' Non-negative Matrix Factorization (NMF)
#'
#' A matrix factorization method that approximates a data matrix as the product of two non-negative matrices, ensuring interpretability through non-negativity constraints.
#'
#' @param object A Ramanome object.
#' @param ncomp The number of components to retain in the NMF decomposition.
#' @param seed The seed for the random number generator.
#' @return A list containing:
#' \item{basis}{The basis matrix of the NMF decomposition (n_cell x n_comp).} 
#' \item{coef}{The coefficient matrix of the NMF decomposition (n_comp x n_wavelength).}
#' while the reconstructed matrix is the product of the basis and coefficient matrices.
#' @export Spectral.Decomposition.Nmf
#' @importFrom NMF nmf
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' decom_nmf <- Spectral.Decomposition.Nmf(data_processed)


Spectral.Decomposition.Nmf <- function(object, n_comp=2, seed=42) {
  dataset <- get.nearest.dataset(object)
  dataset[dataset<0] <- 0
  dataset <- dataset[rowSums(dataset) > 0, ]
  dataset <- dataset[, colSums(dataset) > 0]
  res <- nmf(dataset, n_comp, method="ns", seed=seed)
  return(list(basis=res@fit@W, coef=res@fit@H))
}
