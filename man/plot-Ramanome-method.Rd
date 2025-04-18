% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RamEx-Class.R
\name{plot,Ramanome-method}
\alias{plot,Ramanome-method}
\title{Plot mean spectra for Ramanome objects}
\usage{
\S4method{plot}{Ramanome}(x, y, ...)
}
\arguments{
\item{x}{A Ramanome object}

\item{y}{Not used (required for plot generic)}

\item{...}{Additional arguments passed to mean.spec()}
}
\value{
A plot showing mean spectra and the diversity of each group
}
\description{
Mean spectra of each group in a given Ramanome
The latest spectral matrix in datasets slot and the group information in meta.data slot are used for plotting
}
