#' Read and Process a Single Spectral File
#'
#' This function reads a spectral file from a specified path, checks if the path is valid,
#' and then uses spline interpolation to create a function that can be used to estimate
#' Feature.Reduction.Intensity values at a specified set of wavelengths (xout).
#'
#' @param cell_path The file path to the spectral data file. If the path is invalid or empty,
#' the function will return an empty result.
#' @param xout A vector of wavelengths at which to estimate the Feature.Reduction.Intensity values.
#' @return a dataframe
#' @importFrom data.table fread
#' @importFrom stats splinefun
#' @export read.single
#' @examples
#' # Create a temporary spectral data file
#' temp_file <- tempfile(fileext = ".txt")
#' wave <- seq(500, 3000, by = 10)
#' inten <- sin(wave/500) + rnorm(length(wave), 0, 0.1)
#' write.table(data.frame(wave, inten), temp_file, 
#'             sep = "\t", row.names = FALSE, col.names = FALSE)
#'
#' # Define wavelengths for interpolation
#' xout <- seq(600, 2900, by = 50)
#'
#' # Read and process the spectral file
#' result <- read.single(temp_file, xout)
#'
#' # Clean up
#' unlink(temp_file)
read.single <- function(cell_path,xout){
  if(!is.na(cell_path)&&cell_path!="")
    spec <- fread(cell_path, header = FALSE, sep = "\t")
  else
    return()
  wave <- spec$V1
  inten <- spec$V2
  fn  <- stats::splinefun(wave, inten, "natural")

  return(fn(xout))
}

#' Read and Process Spectral Data Files
#'
#' This function reads spectral data from text files in a specified directory, applies
#' filtering based on wavelength cutoffs, and optionally interpolates the data to
#' a uniform wavelength scale. It returns a Ramanome object containing the processed
#' data.
#'
#' @param data_path The path to the directory containing the spectral data files.
#' @param group.index The index of the group identifier in the filename. Defaults to 1.
#' @param group.levels A vector of levels for the group factor. If NULL, levels will be
#' determined from the data. Defaults to NULL.
#' @param cutoff A numeric vector of two values specifying the wavelength range to
#' include. Defaults to c(500, 3150).
#' @param interpolation A logical value indicating whether to interpolate the data to
#' a uniform wavelength scale. Defaults to FALSE.
#' @return A Ramanome object containing the processed spectral data.
#' @export read.spec
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom data.table fread
#' @importFrom stringr str_split
#' @examples
#' # Create a temporary directory for sample spectral files
#' temp_dir <- tempdir()
#'
#' # Create sample spectral files
#' groups <- c("Control", "Treatment")
#' for(i in 1:4) {
#'   for(g in groups) {
#'     fname <- file.path(temp_dir, paste0(g, "_sample_", i, ".txt"))
#'     wave <- seq(500, 3000, by = 10)
#'     inten <- sin(wave/500) + rnorm(length(wave), 0, 0.1)
#'     write.table(data.frame(wave, inten), fname,
#'                 sep = "\t", row.names = FALSE, col.names = FALSE)
#'   }
#' }
#'
read.spec <- function(data_path, group.index = 1, group.levels = NULL, cutoff = c(500, 3150), interpolation = FALSE) {
  filenames <- list.files(data_path, pattern = "*.txt",  full.names = TRUE,include.dirs = TRUE, recursive = TRUE)
  if(length(grep('Metadata.txt', filenames)) != 0)filenames <- filenames[-grep('Metadata.txt', filenames)]

  num_files <- length(filenames)
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)

  clusterEvalQ(cl, {
    suppressMessages(requireNamespace("data.table", quietly = TRUE))
    suppressMessages(requireNamespace("stringr", quietly = TRUE))
  })

  clusterExport(cl, varlist = c( "read.single", "cut.spec"))

  if (interpolation) {
    wavenumber <- seq(cutoff[1], cutoff[2])

    data_list <- parLapply(cl, filenames, function(x) {
      read.single(x, wavenumber)
    })
  } else {
    if(!is.na(filenames[1])&&filenames[1]!="")
      wavenumber <- cut.spec(fread(filenames[1], header = FALSE, sep = "\t"), cutoff)$V1
    else
      return()

    data_list <- parLapply(cl, filenames, function(x) {
      cut.spec(fread(x, header = FALSE, sep = "\t"), cutoff)$V2
    })
  }

  stopCluster(cl)

  data_mat <- do.call(rbind, data_list)
  colnames(data_mat) <- wavenumber

  group <- str_split(basename(filenames), pattern = "_", simplify = TRUE)[, group.index]
  if (is.null(group.levels)) {
    group.levels <- unique(group)
  }
  group <- factor(group, levels = group.levels)

  meta.data <- data.frame(group = group, filenames = filenames)

  Ramanome <- new("Ramanome", datasets = list(raw.data = data_mat), wavenumber = wavenumber, meta.data = meta.data)

  #show(Ramanome)
  return(Ramanome)
}


#' Read and Process Spectral Data Files
#'
#' This function reads spectral data from text files in a specified directory, applies
#' filtering based on wavelength cutoffs, and optionally interpolates the data to
#' a uniform wavelength scale. It returns a Ramanome object containing the processed
#' data.
#'
#' @param data_path The path to the directory containing the spectral data files.
#' @param group.index The index of the group identifier in the filename. Defaults to 1.
#' @param group.levels A vector of levels for the group factor. If NULL, levels will be
#' determined from the data. Defaults to NULL.
#' @param cutoff A numeric vector of two values specifying the wavelength range to
#' include. Defaults to c(500, 3150).
#' @param interpolation A logical value indicating whether to interpolate the data to
#' a uniform wavelength scale. Defaults to FALSE.
#' @return A Ramanome object containing the processed spectral data.
#' @export read.spec.load
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom data.table fread
#' @importFrom stringr str_split
#' @examples
#' # Create a temporary directory for sample spectral files
#' temp_dir <- tempdir()
#'
#' # Create sample spectral files without parallel processing
#' groups <- c("Control", "Treatment")
#' for(i in 1:4) {
#'   for(g in groups) {
#'     fname <- file.path(temp_dir, paste0(g, "_sample_", i, ".txt"))
#'     wave <- seq(500, 3000, by = 10)
#'     inten <- sin(wave/500) + rnorm(length(wave), 0, 0.1)
#'     write.table(data.frame(wave, inten), fname,
#'                 sep = "\t", row.names = FALSE, col.names = FALSE)
#'   }
#' }
read.spec.load <- function(data_path, group.index = 1, group.levels = NULL, cutoff = c(500, 3150), interpolation = FALSE) {
  filenames <- list.files(data_path, pattern = "*.txt",  full.names = TRUE,include.dirs = TRUE, recursive = TRUE)
  if(length(grep('Metadata.txt', filenames)) != 0)
    filenames <- filenames[-grep('Metadata.txt', filenames)]

  if (interpolation) {
    wavenumber <- seq(cutoff[1], cutoff[2])

    data_list <- lapply(filenames, function(x) {
      read.single(x, wavenumber)
    })
  } else {
    if(!is.na(filenames[1])&&filenames[1]!="")
      wavenumber <- cut.spec(fread(filenames[1], header = FALSE, sep = "\t"), cutoff)$V1
    else
      return()

    data_list <- lapply(filenames, function(x) {
      cut.spec(fread(x, header = FALSE, sep = "\t"), cutoff)$V2
    })
  }

  #stopCluster(cl)

  data_mat <- do.call(rbind, data_list)
  colnames(data_mat) <- wavenumber

  group <- str_split(basename(filenames), pattern = "_", simplify = TRUE)[, group.index]
  if (is.null(group.levels)) {
    group.levels <- unique(group)
  }
  group <- factor(group, levels = group.levels)

  meta.data <- data.frame(group = group, filenames = filenames)

  Ramanome <- new("Ramanome", datasets = list(raw.data = data_mat), wavenumber = wavenumber, meta.data = meta.data)

  #show(Ramanome)
  return(Ramanome)
}
