% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Preprocessing.R
\name{Preprocessing.Background.Remove}
\alias{Preprocessing.Background.Remove}
\title{Subtract background spectra from the cell spectra}
\usage{
Preprocessing.Background.Remove(object, cell.index, cal_mean = FALSE)
}
\arguments{
\item{object}{A Ramanome object.}

\item{cell.index}{The index of the cell component in the filenames, where the background spectra are recorded as 'bg' or 'BG'.}

\item{cal_mean}{Logical value indicating whether to calculate the mean spectra with in a cell when a cell contains multiple spectra.}
}
\value{
A Ramanome object with background-subtracted spectra, named as "sub_data" in the 'datasets' slot
}
\description{
This function performs background subtraction on a Ramanome object by subtracting the background spectra from the dataset.
}
