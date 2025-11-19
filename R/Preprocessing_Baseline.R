#' Baseline Correction via ALS
#'
#' This function performs asymmetric least squares (ALS) baseline correction
#'
#' @param object Ramanome object.
#' @param lambda Log10 smoothing penalty used by the ALS solver (default: 6).
#' @param p Weight for positive residuals (default: 0.05).
#' @param max_iter Maximum ALS iterations per spectrum (default: 20).
#' @param n_threads Number of CPU cores to use. Use \code{NULL} or \code{1} for
#'   single-core execution (default: 1).
#'
#' @return Ramanome object with baseline-corrected spectra stored in
#'   \code{baseline.data}.
#' 
#' @useDynLib RamEx, .registration = TRUE
#' @import RcppEigen
#' @export Preprocessing.Baseline.ALS
#'
#' @examples
#' \dontrun{
#' data(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.ALS(RamEx_data, lambda = 6, p = 0.05)
#' }
Preprocessing.Baseline.ALS <- function(object,
                                       lambda = 6,
                                       p = 0.05,
                                       max_iter = 20,
                                       n_threads = 1) {
  stopifnot(inherits(object, "Ramanome"))

  spectra <- get.nearest.dataset(object)
  if (!is.matrix(spectra)) {
    spectra <- as.matrix(spectra)
  }

  if (is.null(n_threads) || n_threads <= 0) {
    n_threads <- 1
  }

  corrected <- ALSBaselineCpp(
    spectra,
    lambda = lambda,
    p = p,
    maxIter = max_iter,
    n_threads = n_threads
  )
  colnames(corrected) <- object@wavenumber
  object@datasets$baseline.data <- corrected

  object
}

#' Baseline Correction via SNIP
#'
#' This function performs Statistics-Sensitive Non-Linear Iterative
#' Peak-Clipping (SNIP) baseline correction.
#'
#' @param object Ramanome object.
#' @param iterations Maximum SNIP half-window (default: 100).
#' @param decreasing Whether to use the decreasing window variant
#'   (default: \code{TRUE}).
#' @param n_threads Number of CPU cores to use. Use \code{NULL} or \code{1}
#'   for single-core execution (default: 1).
#'
#' @return Ramanome object with baseline-corrected spectra stored in
#'   \code{baseline.data}.
#'
#' @export Preprocessing.Baseline.SNIP
#'
#' @examples
#' \dontrun{
#' data(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.SNIP(
#'   RamEx_data,
#'   iterations = 100,
#'   decreasing = TRUE
#' )
#' }
Preprocessing.Baseline.SNIP <- function(object,
                                        iterations = 100,
                                        decreasing = TRUE,
                                        n_threads = 1) {
  stopifnot(inherits(object, "Ramanome"))

  spectra <- get.nearest.dataset(object)
  if (!is.matrix(spectra)) {
    spectra <- as.matrix(spectra)
  }

  iterations <- as.integer(iterations)
  if (is.na(iterations) || iterations < 0) {
    iterations <- 0L
  }

  if (is.null(n_threads) || n_threads <= 0) {
    n_threads <- 1
  }

  corrected <- SNIPBaselineCpp(
    spectra,
    iterations = iterations,
    decreasing = decreasing,
    n_threads = n_threads
  )
  colnames(corrected) <- object@wavenumber
  object@datasets$baseline.data <- corrected

  object
}
