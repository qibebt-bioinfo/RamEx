% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Quality-Control.R
\name{Qualitycontrol.Mcd}
\alias{Qualitycontrol.Mcd}
\title{Detect outliers using the Minimum Covariance Determinant (MCD) method}
\usage{
Qualitycontrol.Mcd(object, h = 0.5, alpha = 0.01, na.rm = TRUE)
}
\arguments{
\item{object}{Ramanome Object}

\item{h}{The fraction of the data to use for computing the MCD.}

\item{alpha}{The significance level for the MCD.}

\item{na.rm}{A logical value indicating whether to remove NA values.}
}
\value{
A data frame with quality defined by MCD method and their mahalanobis distance from the center.
}
\description{
This function extends the outliers_mcdEst function to detect outliers in a dataset
using the Minimum Covariance Determinant method. It modifies the output to include
additional information such as the distance and center of the outliers.
}
\examples{
data(RamEx_data)
qc_mcd <- Qualitycontrol.Mcd(RamEx_data)
}
