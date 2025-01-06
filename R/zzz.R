#' @title Package Initialization Function
#' @description This function is called when the package is loaded. It loads the required data into the global environment.
#' @param libname The library name
#' @param RamEx The package name
#' @return None
#' @keywords internal
.onLoad <- function(libname, RamEx) {
  data_name <- "RamEx_data"
  data_file <- system.file("extdata", "data", package = "RamEx")
  if (file.exists(data_file)) {
    assign(data_name, read.spec.load(data_file,group.index = 2), envir = .GlobalEnv)
  } else {
    warning("Data file not found.")
  }

}
