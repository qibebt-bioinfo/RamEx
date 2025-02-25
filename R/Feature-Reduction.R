#' Perform PCA Reduction and Plotting
#'
#' This function performs Principal Component Analysis (PCA) on the dataset stored in
#' a Ramanome object and plots the first two principal components. The function can
#' also save the plot as a PNG file.
#'
#' @param object A Ramanome object containing the dataset and metadata.
#' @param draw A logical value indicating whether to draw the PCA plot. Defaults to TRUE.
#' @param save A logical value indicating whether to save the plot as a file. Defaults to FALSE.
#' @return The updated Ramanome object with the PCA results appended to the `reductions` slot.
#' @export Feature.Reduction.Pca
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr %>%
#' @importFrom ggthemes theme_tufte
#' @import ggtext
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data.reduction.pca <- Feature.Reduction.Pca(data_cleaned, draw=TRUE, save = FALSE)

Feature.Reduction.Pca <- function(object, draw = TRUE, save=FALSE) {
  dataset <- get.nearest.dataset(object)
  data.red <- data.frame(prcomp_irlba(dataset, n = 20, center = TRUE, scale. = TRUE)$x[, 1:2])
  names(data.red) <- c('PC 1', 'PC 2')
  
  object@reductions <- append(object@reductions, list(data.red))
  object@reductions$PCA <- data.red
  if (draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point() +
      theme_bw() +
      labs(x = names[1], y = names[2]) +
      theme(
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 15),
        legend.background = element_blank(),
        text = element_text(color = "black"),
        axis.text.x = element_text(size = 15, angle = 0),
        axis.text.y = element_text(size = 15),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15)
      )
    print(plot)
  }
  
  if (FALSE) {
    ggsave(paste("Reduction.pca.png"), plot, width = 8, height = 6)
  }
  
  return(object)
}

#' Perform t-SNE Dimensionality Reduction and Plotting
#'
#' This function performs t-SNE (t-Distributed Stochastic Neighbor Embedding) on a
#' dataset extracted from a Ramanome object. It reduces the dimensionality to two
#' dimensions and creates a plot of the results. The plot can be customized with various
#' themes and labels.
#'
#' @param object A Ramanome object containing the dataset and metadata.
#' @param draw A logical value indicating whether to draw the plot. Defaults to TRUE.
#' @param save A logical value indicating whether to save the plot as a PNG file.
#' Defaults to FALSE.
#' @return The updated Ramanome object with the t-SNE results appended to the reductions slot.
#' @export Feature.Reduction.Tsne
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggsave
#' @import ggtext
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data.reduction.tsne <- Feature.Reduction.Tsne(data_cleaned, draw=TRUE, save = FALSE)
Feature.Reduction.Tsne <- function(object, draw = TRUE, save=FALSE) {
  dataset <- get.nearest.dataset(object)
  
  data.red <- data.frame(Rtsne::Rtsne(
    dataset,
    dims = 2,
    perplexity = 5,
    theta = 0.5,
    verbose = FALSE,
    max_iter = 1000
  )$Y)
  colnames(data.red) <- c("tSNE 1", "tSNE 2")
  object@reductions <- append(object@reductions, list(data.red))
  object@reductions$tSNE <- data.red
  
  if (draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point() +
      theme_bw() +
      labs(x = names[1], y = names[2]) +
      theme(
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 15),
        legend.background = element_blank(),
        text = element_text(color = "black"),
        axis.text.x = element_text(size = 15, angle = 0),
        axis.text.y = element_text(size = 15),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15)
      )
    print(plot)
  }
  
  if (save) {
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
#' @param draw A logical value indicating whether to generate a plot. Defaults to TRUE.
#' @param save A logical value indicating whether to save the plot as a PNG file. Defaults to FALSE.
#' @return The updated Ramanome object with the UMAP reduction results.
#' @export Feature.Reduction.Umap
#' @importFrom uwot umap
#' @importFrom irlba prcomp_irlba
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggthemes theme_tufte
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggsave
#' @import ggtext
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data.reduction.umap <- Feature.Reduction.Umap(data_cleaned, draw=TRUE, save = FALSE)
#'
Feature.Reduction.Umap <- function(object, draw = TRUE, save=FALSE) {
  dataset <- get.nearest.dataset(object)
  data.red.pca <- data.frame(prcomp_irlba(dataset, n = 20, center = TRUE, scale. = TRUE)$x[, 1:20])
  data.red <- data.frame(uwot::umap(data.red.pca, scale = FALSE,  n_threads = detectCores(),n_neighbors = 30, a=0.9922, b=1.112, metric = 'cosine',seed=123))
  colnames(data.red) <- c("UMAP 1", "UMAP 2")
  
  object@reductions <- append(object@reductions, list(data.red))
  object@reductions$UMAP <- data.red
  
  if(draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point() +
      theme_bw() +
      labs(x = names[1], y = names[2]) +
      theme(
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 15),
        legend.background = element_blank(),
        text = element_text(color = "black"),
        axis.text.x = element_text(size = 15, angle = 0),
        axis.text.y = element_text(size = 15),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15)
      )
    print(plot)
  }
  
  if (save) {
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
#' @import ggtext
#' @importFrom stats aov
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggthemes theme_tufte
#' @importFrom ggplot2 ggsave
#' @importFrom ade4 dudi.pco
#' @import ggtext
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data.reduction.pcoa <- Feature.Reduction.Pcoa(data_cleaned, draw=TRUE, save = FALSE)

Feature.Reduction.Pcoa <- function(object, draw = TRUE, save=FALSE) {
  dataset <- get.nearest.dataset(object)
  distance.matrix <- vegdist(dataset, method = "euclidean")
  data.red <- data.frame(dudi.pco(distance.matrix, scannf = FALSE, nf=2)$li)
  names(data.red) <- c('PC 1', 'PC 2')
  
  object@reductions <- append(object@reductions, list(data.red))
  object@reductions$PCoA <- data.red
  if (draw){
    names <- colnames(data.red)
    plot <- ggplot(data.red, aes(get(names[1]), get(names[2]), color = as.factor(object@meta.data$group))) +
      geom_point() +
      theme_bw() +
      labs(x = names[1], y = names[2]) +
      theme(
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 15),
        legend.background = element_blank(),
        text = element_text(color = "black"),
        axis.text.x = element_text(size = 15, angle = 0),
        axis.text.y = element_text(size = 15),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(size = 15)
      )
    print(plot)
  }
  
  if (save) {
    ggsave(paste("Reduction.pcoa.png"), plot, width = 8, height = 6)
  }
  
  return(object)
}


#' Extract Intensity Values at Specified Wavelengths
#'
#' This function retrieves the Feature.Reduction.Intensity values at specified wavelengths or wavelength
#' ranges from a Ramanome object. It uses `confirm.select` to get the spectral data and
#' `confirm.name` to format the names of the selected wavelengths or ranges.
#'
#' @param object A Ramanome object containing spectral data and wavelength information.
#' @param wavenumber A single numeric value, a vector of numeric values, or a list of
#' numeric values representing wavelengths or wavelength ranges.
#' @return The updated Ramanome object with the extracted Feature.Reduction.Intensity values appended to
#' the `interested.bands` slot.
#' @export Feature.Reduction.Intensity
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
Feature.Reduction.Intensity <- function(object, wavenumber) {
  wavenumber <- as.list(wavenumber)
  a <- lapply(wavenumber, function(x) confirm.select(object, x))
  name <- lapply(wavenumber, confirm.name)
  object@interested.bands <- append(object@interested.bands, a)
  names(object@interested.bands) <- c(names(object@interested.bands), unlist(name))
  return(object)
}
