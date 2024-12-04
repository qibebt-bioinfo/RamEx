#' Extract a single spectrum from Raman data
#'
#' This function extracts a single spectrum from Raman data based on a given x value.
#'
#' @param spec A data frame containing the Raman spectrum data.
#' @param x The x value to search for.
#'
#' @return A data frame containing the extracted spectrum.
single <- function(spec, x) {
  print(x)
  return(spec[paste(spec$V1, spec$V2) == x, 3:4])
}

#' Save Renishaw Raman data
#'
#' This function saves the Renishaw Raman data after processing.
#'
#' @param dir_path The directory path where the Raman data is located.
#' @export save.renishaw
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
    aa <- purrr::map(as.list(groups), single, spec = spec)
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
