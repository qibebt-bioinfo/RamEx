\name{RamEx-package}
\alias{RamEx}
\docType{package}
\title{RamEx Package}
\description{
By providing simple-to-use modules for computational tasks that otherwise would require substantial programming experience and algorithmic skills,
RamEx should greatly facilitate the computational mining of ramanomes.More complete information about how to use RamEx can be found in the [github](https://github.com/qibebt-bioinfo/RamEx),
but here you'll find a brief overview of the function and examples to mining of ramanomes:
}
\details{
  This package includes functions listed as below:\cr
 1. read.spec, which Load data from diverse instrument platforms\cr
 2. Ramanome, which Define a standardized format for Ramanomes\cr
 3. Qualitycontrol.ICOD, which Detect outliers among all samples via sudden intensity changes in spectra\cr
 4. Qualitycontrol.Mcd, which The minimum covariance determinant (MCD) method is a highly robust estimator of multivariate location and scatter\cr
 5. Qualitycontrol.T2, which Hoteling’s T2 test is a multivariate statistical test used to determine whether the mean of the two groups is significantly different.\cr
 6. Qualitycontrol.Dis, which The Euclidean distance is a measure of the straight-line distance between two points in Euclidean space.\cr
 7. Qualitycontrol.Snr, which The Signal-to-Noise Ratio (SNR) is a measure used in science and engineering to quantify how much a signal has been corrupted by noise.\cr
 8. Preprocesssing.Background.Remove, which Substract solute background from Raman spectra\cr
 9. Preprocesssing.Background.Spike, which Remove cosmic spikes from Raman spectra\cr
 10. Preprocesssing.Smooth.Sg, which General least-squares smoothing and differentiation by the convolution (Savitzky-Golay) method\cr
 11. Preprocesssing.Normalize, which Normalize the by maximum intensity, intensity at specified wavelength or total area under the curve\cr
 12. Preprocesssing.Cutoff, which Spectral truncation\cr
 13. Preprocesssing.Baseline.Polyfit, which Corrects baseline using polynomial fitting\cr
 14. Preprocesssing.Baseline.Bubble, which Corrects baseline using bubble method, more suitable for liquid phase acquisition spectra\cr
 15. Feature.Reduction.Umap, which A nonlinear dimensionality reduction method that preserves global structure of the data by modeling it as a manifold\cr
 16. Feature.Reduction.Tsne, which A nonlinear dimensionality reduction technique that focuses on preserving local relationships by minimizing the divergence between high-dimensional and low-dimensional pairwise similarities\cr
 17. Feature.Reduction.Pca, which A linear dimensionality reduction method that transforms data into a set of orthogonal components, maximizing variance in the lower-dimensional space\cr
 18. Feature.Reduction.Pcoa, which Similar to PCA but works directly with distance matrices, making it suitable for non-Euclidean data\cr
 19. Feature.Reduction.Intensity, which Extract peak intensities at specific wavelengths\cr
 20. Classification.Rf, which An ensemble learning method that builds multiple decision trees and combines their outputs (e.g., majority voting for classification)\cr
 21. Classification.Lda, which A linear classification method that projects data into a lower-dimensional space to maximize class separability, works better for linearly separable data and small datasets\cr
 22. Classification.Svm, which Finds the optimal hyperplane to separate classes by maximizing the margin between them, powerful for complex, high-dimensional data but computationally expensive\cr
 23. Classification.Gmm, which A probabilistic model that assumes data is generated from a mixture of Gaussian distributions and assigns probabilities to each class\cr
 24. Quantification.Pls, which A regression method that reduces predictors to a smaller set of latent variables while maximizing the covariance between predictors and response variables, ideal for situations where predictors are highly collinear or when the number of predictors exceeds the number of observations\cr
 25. Quantification.Mlr, which Simple and interpretable regression method, but assumes no multicollinearity among predictors and a linear relationship between predictors and response\cr
 26. Quantification.Glm, which Extends linear regression by allowing the dependent variable to follow distributions other tha n normal (e.g., binomial, Poisson) and uses a link function to relate predictors to the response\cr
 27. Raman.Markers.Rbcs, which Identifies important spectral regions using machine learning methods, bands with high contributions to the classifier are considered Raman barcodes\cr
 28. Raman.Markers.Roc, which Identifies Raman markers by measuring the discriminative power of specific bands in distinguishing between given labels\cr
 29. Raman.Markers.Correlations, which Identifies Raman markers by measuring the linear relationships of specific bands, given labels should be meaningful continuous values\cr
 30. Phenotype.Analysis.Louvaincluster, which A community detection algorithm primarily used for clustering in networks or graphs by maximizing modularity (a measure of cluster quality)\cr
 31. Phenotype.Analysis.Kmeans, which A centroid-based clustering algorithm that partitions data into a predefined number of clusters by assigning sample to the nearest center\cr
 32. Phenotype.Analysis.Hca, which Builds a hierarchy of clusters by either merging smaller clusters (agglomerative) or splitting larger ones (divisive) based on a distance metric\cr
 33. Spectral.Decomposition.Mcrals, which Decomposes mixed spectra into pure component spectra and their concentrations using constraints (e.g., non-negativity)\cr
 34. Spectral.Decomposition.Ica, which A statistical method that separates mixed signals into statistically independent components by maximizing non-Gaussianity, effectively separating overlapping signals\cr
 35. Spectral.Decomposition.Nmf, which Decomposes data into non-negative components, ensuring that both the basis and coefficients are non-negative, with less constrains compared to MCR-ALS\cr
 36. Intraramanome.Analysis.Irca.Global, which Analyzes correlations between Raman bands within a Ramanome to identify co-varying spectral features\cr
 37. Intraramanome.Analysis.Irca.Local, which Analyzes correlation between interested Raman bands record by defined attribution\cr
 38. Intraramanome.Analysis.2Dcos, which Captures both synchronous (simultaneous changes) and asynchronous (sequential changes) relationships, providing detailed insights into spectral changes\cr
}
\section{Getting Started}{
  To get started with the RamEx package, you can install it using the following command:
  \preformatted{
      library(devtools)
      install_github("qibebt-bioinfo/RamEx")
  }
  Load the package with:
  \preformatted{
    library(RamEx)
  }
  For more detailed information, see the vignettes:
  \preformatted{
    browseVignettes("RamEx")
  }
}

\section{Datasets}{
  \describe{
    \item{\code{\link{RamEx_data}}}{The "RamEx-data" dataset originates from continuous sampling and monitoring during the liquid fermentation process of Saccharomyces cerevisiae. Samples were taken at the following time points:4, 12, 20, 28 and 36 hours. After three rounds of centrifugation, resuspension, and washing, the cells were air-dried naturally. The spectra were collected using a CAST-R device. The 532 nm laser power was set at 100 mW, with an integration time of 1 second. A 100× objective lens (NA = xx, Olympus, JP) was used to focus the laser beam on the sample.
}
  }
}
\author{
  Gongchao Jing \email{jinggc@qibet.ac.cn}
}

\examples{
  # Example usage of the package
  library(RamEx)
  data(RamEx_data)
  data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
  data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
  data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
  data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
  data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
  data_cleaned <- data_normalized[data_cleaned$index_good,]
}
\keyword{Ramanome, Visualization, Clustering, Software}
