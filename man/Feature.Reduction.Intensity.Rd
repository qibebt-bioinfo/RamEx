% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Feature-Reduction.R
\name{Feature.Reduction.Intensity}
\alias{Feature.Reduction.Intensity}
\title{Extract Intensity Values at Specified Wavelengths or a Range of Wavelengths}
\usage{
Feature.Reduction.Intensity(object, bands)
}
\arguments{
\item{object}{A Ramanome object.}

\item{bands}{A single numeric value, a vector of numeric values, or a list of numeric values representing wavelengths or wavelength ranges.}
}
\value{
The updated Ramanome object with the extracted intensity values appended to the `interested.bands` slot.
}
\description{
Retrieves the intensity values at specified wavelengths or the area under the curve of intensity values at a range of wavelengths from a Ramanome object.
}
\examples{
data(RamEx_data)
data_processed <- Preprocessing.OneStep(RamEx_data)
data_cleaned <- Feature.Reduction.Intensity(data_processed, list(c(2000,2250),c(2750,3050), 1450, 1665)) 
}
