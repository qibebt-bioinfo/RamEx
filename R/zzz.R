source("R/read-spec.R")
source("R/select_value.R")
source("R/rammanome-class.R")
.onLoad <- function(libname, RamEX) {
  data_name <- "RamEx_data"
  data_file <- system.file("extdata", "data", package = "RamEx")
  if (file.exists(data_file)) {
    assign(data_name, read.spec.load(data_file), envir = .GlobalEnv)
  } else {
    warning("Data file not found.")
  }

}
