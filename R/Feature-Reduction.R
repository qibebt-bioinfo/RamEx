#' Perform PCA Reduction and Plotting
#'
#' Performs Principal Component Analysis (PCA) based on the last spectral matrix in the 'datasets' slot
#'
#' @param object A Ramanome object containing the dataset and metadata.
#' @param draw A logical value indicating whether to draw the PCA plot. Defaults to TRUE.
#' @param save A logical value indicating whether to save the plot as a file. Defaults to FALSE.
#' @param n_pc The number of principal components to save. Defaults to 2.
#' 
#' @return The updated Ramanome object with the PCA results appended to the `reductions` slot.
#' @export Feature.Reduction.Pca
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggsave

#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' data.reduction.pca <- Feature.Reduction.Pca(data_processed, draw=TRUE, save = FALSE)

Feature.Reduction.Pca <- function(object, draw = TRUE, save=FALSE, n_pc = 2) {
  dataset <- get.nearest.dataset(object)
  data.red <- data.frame(prcomp_irlba(dataset, n = n_pc, center = TRUE, scale. = TRUE)$x)
  names(data.red) <- paste0('PC ', 1:n_pc)
  
  object@reductions$PCA <- data.red
  if (draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point(alpha = 0.5) +
      labs(x = names[1], y = names[2]) +
      theme_classic()
    print(plot)
  }
  
  if (save) {
    cat('Saving PCA plot to the current working directory: ', getwd(), '\n')
    ggsave(paste("Reduction.pca.png"), plot, width = 8, height = 6)
  }
  
  return(object)
}

#' Perform t-SNE Dimensionality Reduction and Plotting
#'
#' Performs t-SNE (t-Distributed Stochastic Neighbor Embedding) on the last spectral matrix in the 'datasets' slot
#'
#' @param object A Ramanome object containing the dataset and metadata.
#' @param PCA Optional. A numeric value indicating the number of principal components to use. Defaults to 20 (if NULL,raw data matrix will be used to generate t-SNE plot).
#' @param perplexity Controls the balance between local and global structure in the t-SNE algorithm. Typical values range from 5 to 50, with higher values preserving more global structure.
#' @param theta A speed-up parameter for approximate computation (Barnes-Hut algorithm). It ranges from 0 to 1, where smaller values yield more precise results but at higher computational cost (default is 0.5)
#' @param max_iter Controls the balance between local and global structure in the t-SNE algorithm. Typical values range from 5 to 50, with higher values preserving more global structure.
#' @param draw A logical value indicating whether to draw the plot. Defaults to TRUE.
#' @param save A logical value indicating whether to save the plot as a PNG file. Defaults to FALSE.
#' @param seed A numeric value indicating the seed for the random number generator. Defaults to 42.
#' 
#' @return The updated Ramanome object with the t-SNE results appended to the reductions slot.
#' @export Feature.Reduction.Tsne
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggsave

#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' data.reduction.tsne <- Feature.Reduction.Tsne(data_processed, draw=TRUE, save = FALSE)
Feature.Reduction.Tsne <- function(object, PCA=20, perplexity=5, theta=0.5, max_iter=1000,draw = TRUE, save=FALSE, seed=42) {
  set.seed(seed)
  if (!is.null(PCA)) {
    if (is.null(object@reductions$PCA)) {object <- Feature.Reduction.Pca(object, PCA)} 
    dataset <- object@reductions$PCA}
  else {
    dataset <- get.nearest.dataset(object)
  }
  
  data.red <- data.frame(Rtsne::Rtsne(
    dataset,
    dims = 2,
    perplexity = perplexity,
    theta = theta,
    verbose = FALSE,
    max_iter = max_iter
  )$Y)
  colnames(data.red) <- c("tSNE 1", "tSNE 2")
  object@reductions$tSNE <- data.red
  
  if (draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point(alpha = 0.5) +
      labs(x = names[1], y = names[2]) +
      theme_classic()
    print(plot)
  }
  
  if (save) {
    cat('Saving t-SNE plot to the current working directory: ', getwd(), '\n')
    ggsave(paste("Reduction.tsne.png"), plot, width = 8, height = 6)
  }
  
  return(object)
}

#' Perform UMAP Dimensionality Reduction
#'
#' This function performs Uniform Manifold Approximation and Projection (UMAP) on the
#' dataset stored in a Ramanome object. It first applies PCA to reduce the dataset
#' to 20 dimensions and then applies UMAP to reduce it to two dimensions. The reduced
#' data is stored in the object and a plot is generated if specified.
#'
#' @param object A Ramanome object containing the dataset.
#' @param PCA A numeric value indicating the number of principal components to use. Defaults to 20 (if NULL,raw data matrix will be used to generate UMAP plot). 
#' @param n_neighbors Controls the balance between local and global structure in UMAP. Higher values capture more global structure, while lower values focus on local details (typical range: 5–50).
#' @param min.dist Determines the minimum distance between points in the embedding. Lower values (e.g., 0.1) produce tighter clusters, while higher values (e.g., 0.5) allow more spread.
#' @param spread Scales the effective scale of points in the embedding. Works with min.dist to control cluster density—higher values spread points apart, while lower values compress them.
#' @param draw A logical value indicating whether to draw the plot. Defaults to TRUE.
#' @param save A logical value indicating whether to save the plot as a PNG file. Defaults to FALSE.
#' @param seed A numeric value indicating the seed for the random number generator. Defaults to 42.
#' 
#' @return The updated Ramanome object with the UMAP reduction results.
#' @export Feature.Reduction.Umap
#' 
#' @importFrom uwot umap
#' @importFrom irlba prcomp_irlba
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggsave

#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' data.reduction.umap <- Feature.Reduction.Umap(data_processed, draw=TRUE, save = FALSE)
#'
Feature.Reduction.Umap <- function(object, PCA=20, n_neighbors=30, min.dist=0.01,spread=1, draw = TRUE, save=FALSE, seed=42) {
  set.seed(seed)
  if (!is.null(PCA)) {
    if (is.null(object@reductions$PCA)) {object <- Feature.Reduction.Pca(object, PCA)} 
    dataset <- object@reductions$PCA}
  else {
    dataset <- get.nearest.dataset(object)
  }

  data.red <- data.frame(uwot::umap(dataset, scale = FALSE,  n_threads = detectCores(),
  n_neighbors = n_neighbors, min.dist=min.dist, spread=spread, a=0.9922, b=1.112, metric = 'cosine',seed=seed))
  colnames(data.red) <- c("UMAP 1", "UMAP 2")
  
  object@reductions$UMAP <- data.red
  
  if(draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point(alpha = 0.5) +
      labs(x = names[1], y = names[2]) +
      theme_classic()
    print(plot)
  }
  
  if (save) {
    cat('Saving UMAP plot to the current working directory: ', getwd(), '\n')
    ggsave(paste("Reduction.umap.png"), plot, width = 8, height = 6)
  }
  return(object)
}


#' Perform Principal Coordinate Analysis (PCoA) and Plot Results
#'
#' This function performs Principal Coordinate Analysis (PCoA) on the dataset stored in
#' the Ramanome object. It calculates the distance matrix using the 'euclidean' method
#' and then applies PCoA to reduce the dimensionality to two principal components.
#' If drawing is enabled, it creates a plot of the PCoA results, colored by group.
#' The plot can also be saved as a PNG file if specified.
#'
#' @param object A Ramanome object containing the dataset and metadata.
#' @param draw A logical value indicating whether to draw the PCoA plot. Defaults to TRUE.
#' @param save A logical value indicating whether to save the plot as a PNG file. Defaults to FALSE.
#' @return The updated Ramanome object with the PCoA results appended to the `reductions` slot.
#' @export Feature.Reduction.Pcoa
#' @importFrom vegan vegdist
#' @importFrom stats aov
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggsave
#' @importFrom ade4 dudi.pco
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' qc_icod <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[qc_icod$quality,]
#' data.reduction.pcoa <- Feature.Reduction.Pcoa(data_cleaned, draw=TRUE, save = FALSE)

Feature.Reduction.Pcoa <- function(object, draw = TRUE, save=FALSE) {
  dataset <- get.nearest.dataset(object)
  distance.matrix <- vegdist(dataset, method = "euclidean")
  data.red <- data.frame(dudi.pco(distance.matrix, scannf = FALSE, nf=2)$li)
  names(data.red) <- c('PCoA 1', 'PCoA 2')
  
  object@reductions$PCoA <- data.red
  if (draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point(alpha = 0.5) +
      labs(x = names[1], y = names[2]) +
      theme_classic()
    print(plot)
  }
  
  if (save) {
    cat('Saving PCoA plot to the current working directory: ', getwd(), '\n')
    ggsave(paste("Reduction.pcoa.png"), plot, width = 8, height = 6)
  }
  
  return(object)
}


#' Extract Intensity Values at Specified Wavelengths or a Range of Wavelengths
#'
#' Retrieves the intensity values at specified wavelengths or the area under the curve of intensity values at a range of wavelengths from a Ramanome object.
#'
#' @param object A Ramanome object.
#' @param bands A single numeric value, a vector of numeric values, or a list of numeric values representing wavelengths or wavelength ranges.
#' @return The updated Ramanome object with the extracted intensity values appended to the `interested.bands` slot.
#' @export Feature.Reduction.Intensity
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' data_cleaned <- Feature.Reduction.Intensity(data_processed, list(c(2000,2250),c(2750,3050), 1450, 1665)) 
Feature.Reduction.Intensity <- function(object, bands) {
  bands <- as.list(bands)
  a <- lapply(bands, function(x) confirm.select(object, x))
  name <- lapply(bands, confirm.name)
  names(a) <- name
  object@interested.bands <- merge(object@interested.bands, a, by = "row.names", all = TRUE)
  return(object)
}
