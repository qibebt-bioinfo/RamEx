#' Generate time series plots for Raman data
#'
#' This function generates time series plots for Raman data using a specified reduction method (e.g., UMAP).
#'
#' @param object A Ramanome object.
#' @param reduction The reduction method to use (default is "UMAP").
#'
#' @return A combined ggplot object of the time series plots.
#' @export time.series
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
#' @export cluster.color
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
