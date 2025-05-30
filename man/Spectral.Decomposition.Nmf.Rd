% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spectral-Decomposition.R
\name{Spectral.Decomposition.Nmf}
\alias{Spectral.Decomposition.Nmf}
\title{Non-negative Matrix Factorization (NMF)}
\usage{
Spectral.Decomposition.Nmf(
  object,
  n_comp = 2,
  seed = 42,
  max_iter = 100,
  tol = 1e-04,
  verbose = FALSE
)
}
\arguments{
\item{object}{A Ramanome object.}

\item{seed}{The seed for the random number generator.}

\item{ncomp}{The number of components to retain in the NMF decomposition.}
}
\value{
A list containing:
\item{basis}{The basis matrix of the NMF decomposition (n_cell x n_comp).} 
\item{coef}{The coefficient matrix of the NMF decomposition (n_comp x n_wavelength).}
\item{error}{The reconstruction error of the NMF decomposition.}
while the reconstructed matrix is the product of the basis and coefficient matrices.
}
\description{
A matrix factorization method that approximates a data matrix as the product of two non-negative matrices, ensuring interpretability through non-negativity constraints.
Lee & Seung's algorithm is used to implement NMF.
}
\examples{
data(RamEx_data)
data_processed <- Preprocessing.OneStep(RamEx_data)
decom_nmf <- Spectral.Decomposition.Nmf(data_processed)
}
