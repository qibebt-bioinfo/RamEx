% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Raman-Markers.R
\name{save.rf.results}
\alias{save.rf.results}
\title{Save random forest classification results}
\usage{
save.rf.results(result, rf.opts, feature.ids, outdir)
}
\arguments{
\item{result}{The result object returned by the random forest classification.}

\item{rf.opts}{The options used for random forest classification.}

\item{feature.ids}{The ids of the features used for classification.}

\item{outdir}{The directory where the results will be saved.}
}
\value{
None. The function saves multiple files to the specified output directory:
  - A summary file with error metrics and model parameters
  - A file containing classification probabilities for each sample
  - A file with potential mislabeling information
  - A file with feature importance scores
  - A confusion matrix file showing classification performance
}
\description{
This function saves the various results of random forest classification to the specified output directory.
}
