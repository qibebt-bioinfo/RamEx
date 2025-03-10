% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Quality-Control.R
\name{Qualitycontrol.Dis}
\alias{Qualitycontrol.Dis}
\title{Detect outliers in a dataset based on distance from the mean spectrum}
\usage{
Qualitycontrol.Dis(x, min.dis = 1)
}
\arguments{
\item{x}{A matrix of spectral data where each row represents a single spectrum.}

\item{min.dis}{The minimum distance threshold for marking a spectrum as an outlier.
Defaults to 1.}
}
\value{
A data frame containing two columns: 'out' indicating whether each spectrum is an
        outlier, and 'dis' containing the corresponding distance values.
}
\description{
This function identifies outliers in a dataset by calculating the distance of each spectrum
from the mean spectrum of non-outlier spectra. The distance is calculated as the square
root of the sum of the squared differences between each spectrum and the mean spectrum.
Spectra with a distance greater than a specified threshold are marked as outliers.
}
\examples{
data(RamEx_data)
data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
qc_icod <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
qc_dis <- Qualitycontrol.Dis(data_normalized@datasets$normalized.data)
}
