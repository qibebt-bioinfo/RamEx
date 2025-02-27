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
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' qc_icod <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[qc_icod$quality,]
#' #data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' decom_mcrals <- Spectral.Decomposition.Mcrals(data_cleaned, 2)

Spectral.Decomposition.Mcrals<- function(object, ncomp) {
  mcrals <- mdatools::mcrals(object@datasets$normalized.data, ncomp,cont.constraints =list(
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
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' qc_icod <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[qc_icod$quality,]
#' #data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' decom_ica <- Spectral.Decomposition.Ica(data_cleaned, 2)

Spectral.Decomposition.Ica <- function(object, ncomp) {
  ica.result <- icafast(object@datasets$normalized.data, ncomp)
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
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' qc_icod <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[qc_icod$quality,]
#' #data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' decom_nmf <- Spectral.Decomposition.Nmf(data_cleaned)


Spectral.Decomposition.Nmf <- function(object) {
  dataset <- object@datasets$normalized.data
  dataset[dataset<0] <- 0
  dataset <- dataset[rowSums(dataset) > 0, ]
  dataset <- dataset[, colSums(dataset) > 0]
  res <- nmf(dataset, 2, method="ns", seed=123456)
  return(res)
}
