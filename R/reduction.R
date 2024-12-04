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
#' @export reductions.pca
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr %>%
#' @importFrom ggthemes theme_tufte
#' @import ggtext
reductions.pca <- function(object, draw = TRUE, save=FALSE) {
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
        legend.text = element_markdown(size = 15, family = "myFont"),
        legend.background = element_blank(),
        text = element_text(color = "black"),
        axis.text.x = element_text(size = 15, angle = 0, family = "myFont"),
        axis.text.y = element_text(size = 15, family = "myFont"),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(family = "myFont", size = 15)
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
#' @export reductions.tsne
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggsave
#' @import ggtext
reductions.tsne <- function(object, draw = TRUE, save=FALSE) {
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
      legend.text = element_markdown(size = 15, family = "myFont"),
      legend.background = element_blank(),
      text = element_text(color = "black"),
      axis.text.x = element_text(size = 15, angle = 0, family = "myFont"),
      axis.text.y = element_text(size = 15, family = "myFont"),
      axis.ticks = element_line(linewidth = 1),
      axis.ticks.length = unit(0.4, "lines"),
      axis.title = element_text(family = "myFont", size = 15)
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
#' @export reductions.umap
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
reductions.umap <- function(object, draw = TRUE, save=FALSE) {
  dataset <- get.nearest.dataset(object)
  data.red.pca <- data.frame(prcomp_irlba(dataset, n = 20, center = TRUE, scale. = TRUE)$x[, 1:20])
  data.red <- data.frame(uwot::umap(data.red.pca, scale = F,  n_threads = detectCores(),n_neighbors = 30, a=0.9922, b=1.112, metric = 'cosine',seed=123))
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
      legend.text = element_markdown(size = 15, family = "myFont"),
      legend.background = element_blank(),
      text = element_text(color = "black"),
      axis.text.x = element_text(size = 15, angle = 0, family = "myFont"),
      axis.text.y = element_text(size = 15, family = "myFont"),
      axis.ticks = element_line(linewidth = 1),
      axis.ticks.length = unit(0.4, "lines"),
      axis.title = element_text(family = "myFont", size = 15)
    )
  print(plot)
  }

if (draw) {

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
#' @export reductions.pcoa
#' @importFrom vegan vegdist
#' @import ggtext
#' @importFrom stats aov
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggthemes theme_tufte
#' @importFrom ggplot2 ggsave
#' @import ggtext
#'
reductions.pcoa <- function(object, draw = TRUE, save=FALSE) {
  dataset <- get.nearest.dataset(object)
  distance.matrix <- vegdist(dataset, method = "euclidean")
  data.red <- data.frame(dudi.pco(distance.matrix, scannf = F, nf=2)$li)
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
        legend.text = element_markdown(size = 15, family = "myFont"),
        legend.background = element_blank(),
        text = element_text(color = "black"),
        axis.text.x = element_text(size = 15, angle = 0, family = "myFont"),
        axis.text.y = element_text(size = 15, family = "myFont"),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(0.4, "lines"),
        axis.title = element_text(family = "myFont", size = 15)
      )
    print(plot)
  }

  if (FALSE) {
    ggsave(paste("Reduction.pcoa.png"), plot, width = 8, height = 6)
  }

  return(object)
}
