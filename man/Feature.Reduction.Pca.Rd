% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Feature-Reduction.R
\name{Feature.Reduction.Pca}
\alias{Feature.Reduction.Pca}
\title{Perform PCA Reduction and Plotting}
\usage{
Feature.Reduction.Pca(object, draw = TRUE, save = FALSE)
}
\arguments{
\item{object}{A Ramanome object containing the dataset and metadata.}

\item{draw}{A logical value indicating whether to draw the PCA plot. Defaults to TRUE.}

\item{save}{A logical value indicating whether to save the plot as a file. Defaults to FALSE.}
}
\value{
The updated Ramanome object with the PCA results appended to the `reductions` slot.
}
\description{
This function performs Principal Component Analysis (PCA) on the dataset stored in
a Ramanome object and plots the first two principal components. The function can
also save the plot as a PNG file.
}
\examples{
data(RamEx_data)
data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
data_cleaned <- data_normalized[data_cleaned$index_good,]
data.reduction.pca <- Feature.Reduction.Pca(data_cleaned, draw=TRUE, save = FALSE)
}
