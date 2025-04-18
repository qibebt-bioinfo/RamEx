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
  row.names(mcrals$resspec) <- object@wavenumber
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
  row.names(ica.result$M) <- object@wavenumber
  return(ica.result)
}



#' Non-negative Matrix Factorization (NMF)
#'
#' A matrix factorization method that approximates a data matrix as the product of two non-negative matrices, ensuring interpretability through non-negativity constraints.
#' Lee & Seung's algorithm is used to implement NMF.
#' @param object A Ramanome object.
#' @param ncomp The number of components to retain in the NMF decomposition.
#' @param seed The seed for the random number generator.
#' @return A list containing:
#' \item{basis}{The basis matrix of the NMF decomposition (n_cell x n_comp).} 
#' \item{coef}{The coefficient matrix of the NMF decomposition (n_comp x n_wavelength).}
#' \item{error}{The reconstruction error of the NMF decomposition.}
#' while the reconstructed matrix is the product of the basis and coefficient matrices.
#' @export Spectral.Decomposition.Nmf
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' decom_nmf <- Spectral.Decomposition.Nmf(data_processed)


Spectral.Decomposition.Nmf <- function(object, n_comp=2, seed=42, max_iter = 100, tol = 1e-4, verbose = FALSE) {
  dataset <- get.nearest.dataset(object)
  dataset[dataset<0] <- 0
  decomposition <- nmf_lee(dataset, n_comp, max_iter = max_iter, tol = tol, seed = seed, verbose = verbose)
  return(list(basis=decomposition$W, coef=decomposition$H, error=decomposition$error))
}

nmf_lee <- function(X, rank, max_iter = 100, tol = 1e-4, seed = 42, verbose = FALSE) {  
  if (any(X < 0)) stop("Input matrix X must be non-negative.")  
  
  m <- nrow(X)  
  n <- ncol(X) 
  k <- rank     
  
  set.seed(seed)  
  W <- matrix(runif(m * k, min = 0, max = 1), nrow = m, ncol = k)  
  H <- matrix(runif(k * n, min = 0, max = 1), nrow = k, ncol = n)  
  
  for (i in 1:max_iter) {  
    WH <- W %*% H  
    
    H <- H * (t(W) %*% X) / (t(W) %*% WH + 1e-9)  
    W <- W * (X %*% t(H)) / (W %*% H %*% t(H) + 1e-9)  
    error <- norm(X - WH, "F")  
    
    if (verbose) cat("Iteration:", i, "Reconstruction Error:", error, "\n")  
    
    if (error < tol) break  
  }  
  
  list(W = W, H = H, error = error)  
}  