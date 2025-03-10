% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Preprocessing.R
\name{matrix2vector}
\alias{matrix2vector}
\title{Convert Image Matrix to Vector Form}
\usage{
matrix2vector(image, kernel_height, kernel_width, device)
}
\arguments{
\item{image}{A numeric matrix representing the input image}

\item{kernel_height}{Height of the kernel window}

\item{kernel_width}{Width of the kernel window}

\item{device}{Device to use for computation - either "CPU" or "GPU"}
}
\value{
A matrix where each column is a vectorized block from the input image.
        For CPU device, returns a regular matrix.
        For GPU device, returns a gpuMatrix.
}
\description{
This function converts a 2D image matrix into a vectorized form by sliding a kernel window
over the image and extracting blocks. Each block is reshaped into a column vector.
}
