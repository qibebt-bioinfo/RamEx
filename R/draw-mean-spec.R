#' Plot and Save the Mean Spectrum
#'
#' This function generates a plot of the mean spectrum for each group in the dataset
#' contained within a Ramanome object. It uses the `mean.spec` function to calculate
#' the mean spectrum and then saves the plot as a PNG file.
#'
#' @param object A Ramanome object containing the spectral data and metadata.
#' @param gap An optional numeric value specifying the gap between groups in the plot.
#' Defaults to 0.
#' @export spec.mean.draw
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
#' @export mean.spec
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
