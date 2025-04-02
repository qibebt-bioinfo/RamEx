#' Calculate accuracy of prediction matrix
#'
#' This function calculates the accuracy of a prediction matrix.
#'
#' @param pred_matrix The prediction matrix
#'
#' @return A data frame containing the accuracy matrix with the following columns:
#'   \item{true_labels}{The true class labels}
#'   \item{pred_labels}{The predicted class labels}
#'   \item{Freq}{The normalized frequencies/probabilities for each true-predicted label pair}
#' @noRd


confusion.cal <- function(pred_matrix) {
  acc_matrix <- pred_matrix
  acc <- sum(pred_matrix[as.character(pred_matrix[, 1]) == as.character(pred_matrix[, 2]), 3]) / sum(pred_matrix[, 3])
  cat(acc, '\n')
  for (group in unique(pred_matrix[, 1])) {
    acc_matrix[which(pred_matrix[, 1] == group), 3] <- pred_matrix[which(pred_matrix[, 1] == group), 3] / sum(pred_matrix[which(pred_matrix[, 1] == group), 3])
  }
  return(acc_matrix)
}

#' Plot confusion matrix
#'
#' This function plots the confusion matrix using a tile plot.
#'
#' @param true_labels The true labels
#' @param pred_labels The predicted labels
#'
#' @return A ggplot2 object containing a tile plot visualization of the confusion matrix with:
#'   \item{x-axis}{Predicted labels}
#'   \item{y-axis}{True labels}
#'   \item{tile colors}{Normalized prediction frequencies (in percentages)}
#'   \item{text annotations}{Percentage values in each tile}
#' @import ggplot2
#' @import dplyr
confusion.plot <- function(true_labels, pred_labels) {
  pred_matrix <- as.data.frame(table(true_labels, pred_labels))
  colnames(pred_matrix) <- c('true_labels', 'pred_labels', 'Freq')
  pred_matrix <- confusion.cal(pred_matrix)
  pred_matrix <- pred_matrix[order(pred_matrix$true_labels), ]
  pred_matrix$Freq <- round(pred_matrix$Freq, 2) * 100
  pred_matrix <- pred_matrix %>%
    mutate(
      x = factor(pred_labels),
      y = factor(true_labels, levels = rev(unique(true_labels)))
    )
  text_color <- ifelse(pred_matrix$Freq > 1 / 2 * max(pred_matrix$Freq), "white", "black")

  plot_tile <- ggplot(
    pred_matrix,
    aes(x = x, y = y, fill = Freq)
  ) +
    geom_tile(color = "black") +
    theme_bw() +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    geom_text(aes(label = Freq), color = text_color) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  return(plot_tile)
}



#' Plot and Save the Mean Spectrum
#'
#' This function generates a plot of the mean spectrum for each group in the dataset
#' contained within a Ramanome object. It uses the `mean.spec` function to calculate
#' the mean spectrum and then saves the plot as a PNG file.
#'
#' @param object A Ramanome object containing the spectral data and metadata.
#' @param gap An optional numeric value specifying the gap between groups in the plot.
#' Defaults to 0.
#' @return NULL. The function saves the plot as a PNG file named 'mean spec.png' in the current working directory.
#' @importFrom ggplot2 ggsave
#' @noRd


spec.mean.draw <- function(object, gap=0) {
  plot <- mean.spec(get.nearest.dataset(object), object@meta.data$group, gap)
  ggsave(
    'mean spec.png',
    plot,
    width = 10,
    height = 8
  )
}

#' Calculate Mean and Standard Deviation of Spectral Data by Group and Plot
#'
#' Calculates the mean and standard deviation of spectral data for each group then creates a plot with ribbons representing the standard deviation and lines
#'
#' @param data A matrix or data frame containing the spectral data, colnames are wavenumbers.
#' @param group Factor indicating the group for each row in the data.
#' @param gap A numeric value representing the vertical gap between the mean lines of different groups.
#' @return A ggplot2 plot
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs scale_x_continuous theme
#' @export mean.spec
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' mean.spec(data_processed@datasets$normalized.data, data_processed@meta.data$group)

mean.spec <- function(data, group, gap = 0.3) {
  levels <- levels(group)
  group <- as.character(group)
  print(levels)
  data <- as.matrix(data)
  spec_mean <- aggregate(data, by = list(group), mean)
  spec_sd <- aggregate(data, by = list(group), sd)
  n <- nrow(spec_mean)
  i <- which(spec_mean[, 1] == levels[1])
  data_ggplot <- cbind(as.numeric(colnames(data)), t(spec_mean[i, -1] + gap * (n - 1)), t(spec_sd[i, -1]), spec_mean$Group.1[i])
  j <- 2
  for (level in levels[-1]) {
    i <- which(spec_mean[, 1] == level)
    data_ggplot <- rbind(
      data_ggplot,
      cbind(
        as.numeric(colnames(data)),
        t(spec_mean[i, -1] + gap * (n - j)),
        t(spec_sd[i, -1]),
        spec_mean$Group.1[i]
      )
    )
    j <- j + 1
  }

  data_ggplot <- as.data.frame(data_ggplot)
  colnames(data_ggplot) <- c('wavenumber', 'value', 'sd', 'Group')
  data_ggplot$wavenumber <- as.numeric(data_ggplot$wavenumber)
  data_ggplot$value <- as.numeric(data_ggplot$value)
  data_ggplot$sd <- as.numeric(data_ggplot$sd)
  data_ggplot$Group <- factor(data_ggplot$Group, levels = levels)
  plot <- ggplot(
    data = data_ggplot,
    aes(
      x = wavenumber,
      y = value,
      group = Group
    )
  ) +
    geom_ribbon(
      aes(
        ymin = value - sd,
        ymax = value + sd,
        fill = Group
      ),
      alpha = 0.3
    ) +
    geom_line(
      aes(color = Group),
      linewidth = 0.8
    ) +
    labs(y = "Normalized Intensity (a.u.)") +
    xlab(expression(paste("Wavenumber (cm"^{ -1 }, ")"))) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = c(500,1000,1500,2000, 2500,3000,3500)
    ) +
    theme_classic()

  return(plot)
}


#' Retrieve the final spectral matrix from a Ramanome Object
#'
#' This function extracts the most recently added dataset from a Ramanome object.
#' It is assumed that the dataset is stored in the `datasets` slot of the Ramanome object.
#'
#' @param object A Ramanome object.
#' @return A matrix or data frame containing the most recently added dataset from the Ramanome object.
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' data_matrix <- get.nearest.dataset(data_processed)

get.nearest.dataset <- function(object) {
  dataset <- tail(names(object@datasets), 1)
  dataset <- object@datasets[dataset][[1]]
  return(dataset)
}

#' Select a single Raman Feature.Reduction.Intensity value for a given wavenumber
#'
#' This function selects a single Raman Feature.Reduction.Intensity value for a given wavenumber from a Ramanome object.
#'
#' @param object A Ramanome object.
#' @param wave The wavenumber for which to retrieve the Feature.Reduction.Intensity value.
#' @return A numeric vector containing the Feature.Reduction.Intensity values at the specified wavenumber for all spectra in the dataset.
#' @examples
#' # Create a sample Ramanome object
#' wavenumbers <- seq(500, 3500, by = 10)
#' spectra <- matrix(rnorm(100 * length(wavenumbers)), nrow = 100)
#' raman_obj <- new("Ramanome",
#'   datasets = list(raw = spectra),
#'   wavenumber = wavenumbers
#' )
#' @noRd

select.value <- function(object, wave) {
  loc <- which.min(abs(object@wavenumber - wave))
  dataset <- get.nearest.dataset(object)
  return(dataset[, loc])
}

#' Select a range of Raman Feature.Reduction.Intensity values for a given wavenumber range
#'
#' This function selects a range of Raman Feature.Reduction.Intensity values for a given wavenumber range from a Ramanome object.
#'
#' @param object A Ramanome object.
#' @param waves A vector with two values representing the lower and upper bounds of the wavenumber range.
#' @return A numeric vector containing the summed Feature.Reduction.Intensity values across the specified wavenumber range for all spectra in the dataset.
#' @noRd

select.band <- function(object, waves) {
  locs <- object@wavenumber <= waves[2] & object@wavenumber >= waves[1]
  dataset <- get.nearest.dataset(object)
  return(rowSums(dataset[, locs]))
}

#' Format wavelength values into a string
#'
#' This function takes one or two wavelength values and formats them into a string.
#' If two values are provided, they are combined into a range string with a
#' tilde (~) separator. If more than two values are provided, the function stops with
#' an error.
#'
#' @param waves A single numeric value or a vector of two numeric values representing
#' wavelengths or a wavelength range.
#' @return A character string representing either:
#'   \item{Single wavelength}{The wavelength value as a string}
#'   \item{Wavelength range}{Two wavelength values separated by a tilde (~)}
#' @noRd 

confirm.name <- function(waves) {
  if (length(waves) == 1) {
    name <- as.character(waves)
  } else if (length(waves) == 2) {
    name <- paste(waves[1], waves[2], sep = '~')
  } else {
    stop('Error! Please input a wavenumber or a Raman band!')
  }
  return(name)
}

#' Confirm Selection and Retrieve Spectral Data
#'
#' This function determines the type of selection based on the input and retrieves the
#' corresponding spectral data from a Ramanome object. If a single wavelength is
#' provided, it uses `select.value` to retrieve the data at that wavelength. If a
#' range of wavelengths is provided, it uses `select.band` to retrieve the sum of
#' the data within that range.
#'
#' @param object A Ramanome object containing spectral data and wavelength information.
#' @param waves A single numeric value representing a wavelength or a vector of two
#' numeric values representing a wavelength range.
#' @return A numeric vector containing either:
#'   \item{Single wavelength}{Feature.Reduction.Intensity values at the specified wavelength}
#'   \item{Wavelength range}{Summed Feature.Reduction.Intensity values across the specified range}
#' @seealso select.value for retrieving data at a single wavelength.
#' @seealso select.band for retrieving the sum of data within a wavelength range.
#' @noRd 

confirm.select <- function(object, waves) {
  if (length(waves) == 1) {
    values <- select.value(object, waves)
  } else if (length(waves) == 2) {
    values <- select.band(object, waves)
  }else { print('Error! Please input a wavenumber or a Raman band!') }
  return(values)
}

#' Generate time series plots for Raman data
#'
#' This function generates time series plots for Raman data using a specified reduction method (e.g., UMAP).
#'
#' @param object A Ramanome object.
#' @param reduction The reduction method to use (default is "UMAP").
#'
#' @return A cowplot object containing multiple plots arranged in a grid, with each plot showing the time series
#'         visualization for a different group in the dataset. The function also saves the plot as a PNG file
#'         named 'Time.series_[reduction].png'.
#' @importFrom dbscan dbscan
#' @importFrom stats quantile
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggsave


time.series <- function(object, reduction = 'UMAP') {
  dataset <- get.nearest.dataset(object)
  groups <- unique(object@meta.data$group)
  clusters <- paste('Cluster_', dbscan::dbscan(dataset, eps = quantile(dist(dataset), 0.01))$cluster)
  data.red <- data.frame(
    object@reductions[[reduction]],
    group = object@meta.data$group,
    cluster = clusters
  )
  plots <- lapply(groups, function(x) return(cluster.color(data.red, x)))
  plot <- cowplot::plot_grid(plotlist = plots, ncol = length(groups))
  ggsave(
    paste('Time.series_', reduction, '.png'),
    plot,
    width = 8,
    height = 6
  )
  return(plot)
}


#' Generate a color-coded plot for cluster visualization
#'
#' This function generates a color-coded plot to visualize clusters in a dataset.
#'
#' @param data A data frame containing the dataset.
#' @param group The group to focus on.
#' @return A ggplot2 object containing:
#'   \item{Points}{Individual data points colored by cluster}
#'   \item{Hulls}{Convex hulls around each cluster}
#'   \item{Colors}{Distinct colors for each cluster}
#'   \item{Theme}{Minimalist theme with removed grid and axis elements}
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggforce geom_mark_hull
#' @importFrom dplyr filter
#' @importFrom scales hue_pal
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual

cluster.color <- function(data, group) {
  ncluster <- length(unique(data$cluster))
  colors <- scales::hue_pal(direction = -1)(ncluster)
  names <- names(data)

  ggplot(data, aes_string(names[1], names[2])) +
    geom_point(data = dplyr::filter(data, group != group), color = 'grey90', size = 0.9) +
    ggforce::geom_mark_hull(
      aes(fill = cluster, color = cluster),
      alpha = 0.1,
      expand = unit(2, 'mm'),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_fill_manual(values = colors) +
    geom_point(
      data = dplyr::filter(data, group == group),
      aes(color = cluster, fill = cluster)
    ) +
    scale_color_manual(values = colors) +
    labs(x = '', y = '', title = group) +
    theme_classic()
}


#' Extract a single spectrum from Raman data
#'
#' This function extracts a single spectrum from Raman data based on a given x value.
#'
#' @param spec A data frame containing the Raman spectrum data.
#' @param x The x value to search for.
#'
#' @return A data frame containing columns 3 and 4 of the input data frame for the rows where
#'         the concatenated values of V1 and V2 match the input x value.
#' @noRd

single <- function(spec, x) {
  return(spec[paste(spec$V1, spec$V2) == x, 3:4])
}

#' Save Renishaw Raman data
#'
#' This function saves the Renishaw Raman data after processing.
#'
#' @param dir_path The directory path where the Raman data is located.
#' @return NULL. The function processes the input text files and saves new text files in the same directory
#'         with modified names based on the grouping of the data. The original input files are removed.
#' @importFrom purrr map
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom dplyr %>%

save.renishaw <- function(dir_path) {
  setwd(dir_path)
  filenames <- list.files(
    dir_path,
    pattern = ".txt",
    full.names = TRUE,
    include.dirs = TRUE,
    recursive = TRUE
  )
  for (j in seq_along(filenames)) {
    filename <- filenames[j] %>% gsub('.txt', '', .)
    spec <- fread(filenames[j], header = FALSE, sep = "\t")
    groups <- base::unique(paste(spec$V1, spec$V2))
    aa <- purrr::map(as.list(groups), single, spec = spec)
    names(aa) <- seq_along(aa)
    for (i in seq_along(aa)) {
      fwrite(
        aa[[i]][order(aa[[i]][, 1])],
        file = paste(filename, groups[i], i, '.txt', sep = '_'),
        sep = '\t',
        col.names = FALSE
      )
    }
    file.remove(filenames[j])
  }
}

