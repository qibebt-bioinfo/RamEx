#' Calculate accuracy of prediction matrix
#'
#' This function calculates the accuracy of a prediction matrix.
#'
#' @param pred_matrix The prediction matrix
#'
#' @return The accuracy matrix

confusion.cal <- function(pred_matrix) {
  acc_matrix <- pred_matrix
  acc <- sum(pred_matrix[as.character(pred_matrix[, 1]) == as.character(pred_matrix[, 2]), 3]) / sum(pred_matrix[, 3])
  cat(acc, '\n')
  for (group in unique(pred_matrix[, 1])) {
    acc_matrix[which(pred_matrix[, 1] == group), 3] <- pred_matrix[which(pred_matrix[, 1] == group), 3] / sum(pred_matrix[which(pred_matrix[, 1] == group), 3])
  }
  return(acc_matrix)
}

#' Plot prediction results
#'
#' This function plots the prediction results using a tile plot.
#'
#' @param true_labels The true labels
#' @param pred_labels The predicted labels
#'
#' @return A tile plot of the prediction results
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
#' @importFrom ggplot2 ggsave

spec.mean.draw <- function(object, gap=0) {
  plot <- mean.spec(get.nearest.dataset(object), object@meta.data$group, gap)
  ggsave(
    'mean spec.png',
    plot,
    width = 10,
    height = 8
  )
}

#' Calculate Mean and Standard Deviation of Spectral Data by Group
#'
#' This function calculates the mean and standard deviation of spectral data for each group.
#' It then creates a ggplot2 plot with ribbons representing the standard deviation and lines
#' representing the mean for each group.
#'
#' @param data A matrix or data frame containing the spectral data.
#' @param group A factor or character vector indicating the group for each row in the data.
#' @param gap A numeric value representing the vertical gap between the mean lines of different groups.
#' @return A ggplot object representing the plot.
#' @importFrom hyperSpec aggregate
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line theme_bw labs scale_x_continuous theme
#' @importFrom stats sd

mean.spec <- function(data, group, gap = 0.3) {
  levels <- levels(group)
  group <- as.character(group)
  print(levels)
  data <- as.matrix(data)
  spec_mean <- hyperSpec::aggregate(data, by = list(group), mean)
  spec_sd <- hyperSpec::aggregate(data, by = list(group), sd)
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
    theme_bw() +
    labs(y = "Normalized Intensity (a.u.)") +
    xlab(expression(paste("Wavenumber (cm"^{ -1 }, ")"))) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = c(500,1000,1500,2000, 2500,3000,3500)
    ) +
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15),
      legend.background = element_blank(),
      text = element_text(color = "black"),
      axis.title.y = element_text(size = 20, angle = 90),
      axis.text.x = element_text(size = 15, angle = 0),
      axis.text.y = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 20),
      axis.ticks = element_line(linewidth = 1),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.4, "lines"),
      axis.title = element_text(size = 20)
    )

  return(plot)
}

#' Generate Raman imaging
#'
#' This function generates a Raman imaging plot for a specified peak in a Ramanome object.
#'
#' @param object The Ramanome object
#' @param peak The specified peak for Raman imaging
#'
#' @importFrom grDevices jpeg
#' @importFrom graphics image
#' @importFrom rlist list.map
#' @importFrom grDevices dev.off

image.peak <- function(object, peak) {
  data <- data.frame(
    x = as.numeric(object@meta.data$x),
    y = as.numeric(object@meta.data$y),
    value = object@interested.bands[[as.character(peak)]] * 100
  )
  data.matrix <- lapply(
    unique(data$y), function(loc) {
      rr.data <- data[data$y == loc,]
      return(rr.data$value[order(rr.data$x)])
    }
  )
  data.matrix <- do.call(rbind, rlist::list.map(data.matrix, .))
  grDevices::jpeg(
    filename = paste0('Raman_Imaging_', peak, ' .jpeg'),
    width = length(unique(data$x)) * 30,
    height = length(unique(data$x)) * 30,
    quality = 150,
    res = 200
  )
  graphics::image(data.matrix, axes = F)
  grDevices::dev.off()
}


#' Retrieve the Nearest Dataset from a Ramanome Object
#'
#' This function extracts the most recently added dataset from a Ramanome object.
#' It is assumed that the dataset is stored in the `datasets` slot of the Ramanome object.
#'
#' @param object A Ramanome object that contains datasets in its `datasets` slot.
#' @return The most recently added dataset from the Ramanome object.


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

#' @return The Raman Feature.Reduction.Intensity value for the given wavenumber.
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

#' @return The Raman Feature.Reduction.Intensity values within the given wavenumber range.
select.band <- function(object, waves) {
  locs <- object@wavenumber <= waves[2] & object@wavenumber >= waves[1]
  dataset <- get.nearest.dataset(object)
  return(rowSums(dataset[, locs]))
}

#' string. If two values are provided, they are combined into a range string with a
#' tilde (~) separator. If more than two values are provided, the function stops with
#' an error.
#'
#' @param waves A single numeric value or a vector of two numeric values representing
#' wavelengths or a wavelength range.
#' @return A character string representing the formatted wavelength name or range.


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
#' @return The selected spectral data based on the input wavelength or range.

#' @seealso select.value for retrieving data at a single wavelength.
#' @seealso select.band for retrieving the sum of data within a wavelength range.
#'
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
#' @return A combined ggplot object of the time series plots.

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

#' Custom theme for time series plots
#'
#' This function defines a custom theme for time series plots.
#'
SelfTheme <- function() {
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    strip.background = element_rect(fill = "white")
  )
}

#' Generate a color-coded plot for cluster visualization
#'
#' This function generates a color-coded plot to visualize clusters in a dataset.
#'
#' @param data A data frame containing the dataset.
#' @param group The group to focus on.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggforce geom_mark_hull
#' @importFrom dplyr filter
#' @importFrom scales hue_pal
#' @importFrom ggthemes theme_tufte
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @return A ggplot object.
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
    ggthemes::theme_tufte(ticks = FALSE) +
    scale_color_manual(values = colors) +
    labs(x = '', y = '', title = group) +
    theme_bw() +
    SelfTheme
}


#' Extract a single spectrum from Raman data
#'
#' This function extracts a single spectrum from Raman data based on a given x value.
#'
#' @param spec A data frame containing the Raman spectrum data.
#' @param x The x value to search for.
#'
#' @return A data frame containing the extracted spectrum.
single <- function(spec, x) {
  return(spec[paste(spec$V1, spec$V2) == x, 3:4])
}

#' Save Renishaw Raman data
#'
#' This function saves the Renishaw Raman data after processing.
#'
#' @param dir_path The directory path where the Raman data is located.
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
    include.dirs = T,
    recursive = T
  )
  for (j in seq_along(filenames)) {
    filename <- filenames[j] %>% gsub('.txt', '', .)
    spec <- fread(filenames[j], header = FALSE, sep = "\t")
    groups <- base::unique(paste(spec$V1, spec$V2))
    aa <- map(as.list(groups), single, spec = spec)
    names(aa) <- seq_along(aa)
    for (i in seq_along(aa)) {
      fwrite(
        aa[[i]][order(aa[[i]][, 1])],
        file = paste(filename, groups[i], i, '.txt', sep = '_'),
        sep = '\t',
        col.names = F
      )
    }
    file.remove(filenames[j])
  }
}


#' Cut Spectral Data Based on Wavelength Range
#'
#' This function extracts a subset of spectral data based on a specified wavelength range.
#' It filters the data to include only the rows where the wavelength (V1) is within the
#' lower and upper bounds of the cutoff range.
#'
#' @param single_spec A data frame containing spectral data with wavelengths in the first column (V1).
#' @param cutoff A numeric vector of two values specifying the lower and upper bounds of the wavelength range.
#' @return A data frame containing the subset of spectral data within the specified wavelength range.
#' @export cut.spec


cut.spec <- function(single_spec, cutoff) {
  return(single_spec[single_spec$V1 < cutoff[2] & single_spec$V1 > cutoff[1]])
}
