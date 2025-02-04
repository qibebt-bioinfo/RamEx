% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spectral-Decomposition.R
\name{Spectral.Decomposition.Mcrals}
\alias{Spectral.Decomposition.Mcrals}
\title{Perform Multivariate Curve Resolution Alternating Least Squares (MCR-ALS)}
\usage{
Spectral.Decomposition.Mcrals(object, ncomp)
}
\arguments{
\item{object}{A Ramanome object containing the dataset.}

\item{ncomp}{The number of components to estimate.}
}
\value{
A list containing the estimated components and concentrations.
}
\description{
This function applies the MCR-ALS algorithm to the dataset stored in a Ramanome object.
MCR-ALS is a chemometric technique used to resolve overlapping spectra into their
constituent pure components. The function returns the estimated components and
concentrations.
}
\examples{
data(RamEx_data)
data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
data_cleaned <- data_normalized[data_cleaned$index_good,]
data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
decom_mcrals <- Spectral.Decomposition.Mcrals(data_cleaned, 2)
}
