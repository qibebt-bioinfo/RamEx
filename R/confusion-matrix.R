#' Calculate accuracy of prediction matrix
#'
#' This function calculates the accuracy of a prediction matrix.
#'
#' @param pred_matrix The prediction matrix
#'
#' @return The accuracy matrix
#' @export confusion.cal
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
#' @export confusion.plot
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
