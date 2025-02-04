% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RamEx-Basic-Funciton.R
\name{single}
\alias{single}
\title{Extract a single spectrum from Raman data}
\usage{
single(spec, x)
}
\arguments{
\item{spec}{A data frame containing the Raman spectrum data.}

\item{x}{The x value to search for.}
}
\value{
A data frame containing columns 3 and 4 of the input data frame for the rows where
        the concatenated values of V1 and V2 match the input x value.
}
\description{
This function extracts a single spectrum from Raman data based on a given x value.
}
