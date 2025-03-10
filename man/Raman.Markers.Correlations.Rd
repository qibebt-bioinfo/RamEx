% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Raman-Markers.R
\name{Raman.Markers.Correlations}
\alias{Raman.Markers.Correlations}
\title{Compute Correlations with a Target Variable}
\usage{
Raman.Markers.Correlations(X, y, min.cor = 0.8)
}
\arguments{
\item{X}{A data frame containing the dataset.}

\item{y}{A vector containing the target variable.}

\item{min.cor}{The minimum correlation threshold. Defaults to 0.8.}
}
\value{
A list containing the high correlation elements, individual correlations,
correlations for two-variable combinations, and correlations for three-variable combinations.
}
\description{
This function calculates the correlations between each column of a dataset X and a target
variable y. It also computes the correlations for combinations of two and three
variables, after applying certain transformations.
}
\examples{
data(RamEx_data)
data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
qc_icod <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
data_cleaned <- data_normalized[qc_icod$quality,]
data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
options(mc.cores = 2)
#cor_markers <- Raman.Markers.Correlations(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],as.numeric(data_cleaned@meta.data$group), min.cor = 0.6)
}
