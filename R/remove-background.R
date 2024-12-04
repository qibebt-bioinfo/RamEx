#' Perform background subtraction on Ramanome object
#'
#' This function performs background subtraction on a Ramanome object by subtracting the background spectra from the dataset.
#'
#' @param object A Ramanome object.
#' @param cell.index The index of the cell component in the filenames.
#' @param cal_mean Logical value indicating whether to calculate the mean background spectrum.
#'
#' @return The modified Ramanome object.
#' @importFrom stringr str_extract
#' @importFrom stringr str_split_i
#' @importFrom stringr str_detect
#' @importFrom hyperSpec aggregate
#' @importFrom rlist list.map
#' @export background.remove
background.remove <- function(object, cell.index, cal_mean = FALSE) {

  # Extract group and cell information from filenames
  group <- str_extract(object@meta.data$filenames, pattern = paste0('([^_]*_){', cell.index - 2, '}[^_]*'))
  cell <- str_split_i(object@meta.data$filenames, pattern = '_', cell.index)
  object@meta.data$cell <- cell

  if (cal_mean) {
    # Calculate mean background spectrum
    data <- get.nearest.dataset(object)
    data.mean <- hyperSpec::aggregate(data, by = list(group, cell), mean)
    object@datasets$mean_data <- as.matrix(data.mean[, -c(1, 2)])
    object@meta.data <- data.frame(
      group = str_split_i(data.mean[, 1], pattern = '_', 1),
      filenames = paste(data.mean[, 1], data.mean[, 2], sep = '_'),
      cell = data.mean[, 2]
    )
    group <- data.mean[,1]
    cell <- data.mean[,2]

  }

  # Perform background subtraction
  data.sub <- tapply(object, group, function(x) {
    data <- get.nearest.dataset(x)
    row.names(data) <- x@meta.data$filenames
    bg_spec <- data[str_detect(x@meta.data$cell, '[Bb][Gg]'), ]
    if (length(bg_spec) == 0) {
      return(data)
    } else if (is.null(nrow(bg_spec))) {
      return(t(t(data) - bg_spec))
    } else {
      return(t(t(data) - colMeans(bg_spec)))
    }
  }, simplify = FALSE, default = 0)
  data.sub <- do.call(rbind, list.map(data.sub, as.matrix(unlist(.))))

  # Update datasets and meta.data in Ramanome object
  object@datasets$sub_data <- data.sub
  filenames <- rownames(data.sub)
  object@meta.data <- data.frame(
    group = str_split_i(filenames, pattern = '_', 1),
    filenames = filenames
  )

  # Remove background samples from the Ramanome object
  return(object[!str_detect(str_split_i(filenames, pattern = '_', cell.index), '[Bb][Gg]')])
}
