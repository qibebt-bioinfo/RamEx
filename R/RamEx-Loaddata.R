#' Read and Process a Single Spectral File
#'
#' Reads a spectral data file (e.g., .txt or .csv), automatically detects whether
#' samples are arranged by columns or rows, determines file structure type 
#' (single spectrum, mapping matrix, or coordinate-based data), and performs
#' optional spline interpolation within a specified wavenumber range.
#'
#'
#' @param filepath Character. Path to the spectral file to be read.
#' @param cutoff Numeric vector of length 2.
#'    Defines the wavenumber range (`[min, max]`) for trimming or interpolation.
#' @param interpolation Logical. If `TRUE`, applies spline interpolation across
#'   a uniform wavenumber grid defined by `cutoff`. Default: `FALSE`.
#' @return A list with components:  
#'   \item{wavenumber}{Numeric vector of wavenumber values.}  
#'   \item{intensity}{Matrix of intensity values. Each row corresponds to a sample.}  
#'   \item{orient}{Character. Indicates whether data were read as "row" or "col".}
#' @importFrom data.table fread
#' @importFrom stats splinefun
#' @noRd
read_spectral_file <- function(filepath, sep, cutoff = c(500, 3150), interpolation = FALSE) {
  if (is.na(filepath) || filepath == "" || !file.exists(filepath)) return(NULL)
  ext <- tools::file_ext(filepath)
  if(is.null(sep)) sep <- ifelse(ext == "csv", ",", "\t")
  spec <- fread(filepath, header = FALSE, sep = sep)
  spec <- as.matrix(spec)
  n_cols <- ncol(spec)
  n_rows <- nrow(spec)

  col_is_wave <- apply(spec,2, function(cl) is_wavenumber(as.numeric(cl)))
  if (any(col_is_wave) ){
    orient <- 'col'
  } else if(is_wavenumber(as.numeric(spec[1,])) ){
    spec <- t(spec)
    col_is_wave <- c(TRUE, rep(FALSE, ncol(spec)-1))
    orient <- 'row'
  } else {
    orient <- 'unknown'
    return(list( orient = orient))}

  if (sum(col_is_wave) == 1){
    wave_col <- which(col_is_wave)
    before_cols <- seq_len(wave_col - 1)
    after_cols  <- seq(from = wave_col + 1, to = ncol(spec))
    wave <- as.numeric(spec[,wave_col])

    if (length(after_cols) == 1 && length(before_cols) == 0) {
      # Type 1 : Single spectrum, col 1 represents wavenumber, col 2 represents intensity 
      inten_mat <- matrix(spec[,after_cols], nrow = 1)
    } else if (length(after_cols) > 1 && length(before_cols) == 0) {
      # Type 2 ：Mapping matrix. with col 1 represents wavenumber, others represent intensity of multipul cell
      inten_mat <- t(as.matrix(spec[, after_cols]))
    } else if (length(before_cols) > 0 && length(after_cols) == 1) {
      # Type 3: coordinate scan format
      coord_names <- paste0("coord", before_cols)
      spec <- as.data.frame(spec)
      colnames(spec) <- c(coord_names, "wave", "inten")
      
      split_vars <- as.data.frame(spec[, coord_names, drop = FALSE])
      key_strings <- do.call(paste, c(as.list(split_vars), sep = "_"))
      spec_list <- split(spec, key_strings, drop = TRUE)
      wave <- sort(unique(spec$wave))

      inten_mat <- do.call(rbind, lapply(spec_list, function(df) {
        df <- df[match(wave, df$wave), , drop = FALSE]
        df$inten
      }))
    } else {
      warning(sprintf("Unrecognized column arrangement in file: %s", basename(filepath)))
      return(list( orient = 'unknown'))
    }
  } else {
    warning(sprintf("No valid wavenumber column detected in file: %s", basename(filepath)))
    return(list( orient = 'unknown'))
  }

  if (interpolation) {
    xout <- seq(cutoff[1], cutoff[2])
    inten_mat <- t(apply(inten_mat, 1, function(inten) {
      fn <- stats::splinefun(wave, inten, "natural")
      fn(xout)
    }))
    wave <- xout
  } else {
    valid_idx <- which(wave >= cutoff[1] & wave <= cutoff[2])
    wave <- wave[valid_idx]
    inten_mat <- inten_mat[, valid_idx]
  }
  
  res <- list(
    wavenumber = wave,
    intensity = matrix(inten_mat, ncol=length(wave)),
    orient = orient
  )
  
  return(res)
}

#' Remove common characters in strings
#' @noRd
strip_common_prefix_suffix <- function(strings) {
  if (length(strings) <= 1) return(strings)
  
  prefix <- substring(strings[1], 1,
                      min(sapply(strings, function(x)
                        which(strsplit(x, "")[[1]] != strsplit(strings[1], "")[[1]])[1])) - 1)
  if (is.na(prefix)) prefix <- ""
  
  reversed <- sapply(strings, function(x) paste0(rev(strsplit(x, "")[[1]]), collapse = ""))
  suffix <- substring(reversed[1], 1,
                      min(sapply(reversed, function(x)
                        which(strsplit(x, "")[[1]] != strsplit(reversed[1], "")[[1]])[1])) - 1)
  if (is.na(suffix)) suffix <- ""
  
  suffix <- paste0(rev(strsplit(suffix, "")[[1]]), collapse = "")
  
  stripped <- gsub(paste0("^", stringr::fixed(prefix)), "", strings)
  stripped <- gsub(paste0(stringr::fixed(suffix), "$"), "", stripped)

  if (anyDuplicated(stripped)) {
    dup_counts <- ave(seq_along(stripped), stripped, FUN = seq_along)
    dup_suffix <- ifelse(dup_counts > 1, paste0("_", dup_counts), "")
    stripped <- paste0(stripped, dup_suffix)
  }
  
  stripped
}

#' Judge if a given vector is wavenumber information
#' @noRd
is_wavenumber <- function(x) {
  diffs <- diff(x)
  sign_diff <- sign(diffs[!is.na(diffs)])
  monotonic_ratio <- max(mean(sign_diff > 0), mean(sign_diff < 0))
  is_monotonic <- monotonic_ratio > 0.95
  rng <- range(x, na.rm = TRUE)
  is_range <- rng[1] > 50 && rng[2] < 6000
  is_monotonic & is_range
}

#' Print orient
#' @noRd
print_orient <- function(orient_list, filenames){
  if (all(orient_list == 'col') | all(orient_list == 'row')){
    message(sprintf("All spectra were read as %s-oriented.", unique(orient_list)))
  } else if(all(orient_list %in% c('col','row'))){
    message(sprintf("%d files read as column-oriented samples, %d as row-oriented.",
                    sum(orient_list == "col"), sum(orient_list == "row")))
  } else {
    message(sprintf("%d files have unknown orientation. Please verify file format.", length(which(orient_list == 'unknown'))))
    message("Affected files (first few shown): ",
            paste(head(filenames[orient_list == "unknown"]), collapse = ", "))
  }
}

#' Read Spectral Data Files from a Directory
#'
#' Reads multiple spectral data files from a target directory and automatically
#' identifies their structure — whether single spectra, mapping matrices, and coordinate‑based scans —
#' and merges them into a unified data matrix.
#' 
#' It optionally performs spline interpolation across the specified
#' wavenumber range and constructs a `Ramanome` object containing both
#' spectral data and associated metadata (such as file names and group labels).
#'
#' **Note:**  
#' If all files are reported as having `unknown` orientation or structure,
#' please verify your file format manually. Such cases may not be supported
#' by the automatic reader, and you can instead use
#' [build_ramanome_object()] to manually assemble your own `Ramanome` object
#' from pre‑processed data matrices.
#' 
#' @param data_path Character. 
#'   Path to the directory containing spectral data files.
#'   The function searches recursively for files matching the specified extensions (see `file_type`).
#'   Files named "Metadata" will be ignored.
#'
#' @param file_type Character vector.  
#'   Specifies which file formats to include.  
#'   Supported values include: "txt", "csv", "tsv"
#' 
#' @param group.index Integer. 
#'   Index specifying which split element is used as the group label.
#'   Splitting is performed according to the pattern set in `group_splits`.
#' 
#' @param group.levels Character vector. 
#'   Optional predefined levels for the grouping factors.  
#'   If `NULL` (default), unique group names will be derived automatically from file names.
#'
#' @param group_splits Character.  
#'   Regular expression defining how to split file paths for group extraction.  
#'   Default is `"/|_"`, meaning both directory `/` and underscore `_` are used as separators.
#'
#' @param cutoff Numeric vector of length 2.  
#'   Defines the wavenumber range (`[min, max]`) used for trimming or interpolation.  
#'   Default: `c(500, 3150)`.
#'
#' @param interpolation Logical.  
#'   If `TRUE`, performs natural spline interpolation across uniformly spaced wavenumbers
#'   within `cutoff`.  
#'   If `FALSE`, the spectra are directly truncated to the specified range without interpolation.
#'
#' @param n_cores Integer.  
#'   Number of CPU cores to use for parallel reading.  
#'   If `n_cores > 1`, the function runs in parallel;  
#'   if `1`, it falls back to single‑core serial processing.  
#'   Default: `1`.
#'
#' @return
#' A `Ramanome` object containing:
#' \itemize{
#'   \item \code{datasets$raw.data} — a numeric matrix of dimension (n_samples × n_wavenumbers)
#'   \item \code{wavenumber} — the common wavenumber vector
#'   \item \code{meta.data} — a data frame storing group labels and source filenames
#' }
#'
#' @details
#' File structure detection rules:
#' \itemize{
#'   \item **Type 1:** Two columns (wavenumber, intensity)
#'   \item **Type 2:** Mapping matrix — first column is wavenumber, remaining columns are multiple spectra
#'   \item **Type 3:** Coordinate scan — columns 1–2 are (x, y) positions, column 3 is wavenumber, column 4 is intensity
#' }
#'  
#' The function automatically determines the file type and extracts the spectra accordingly.
#'
#' @examples
#' \dontrun{
#' RamanData <- read.spec(
#'   data_path = "data/",
#'   group.index = 1,
#'   cutoff = c(600, 3000),
#'   interpolation = TRUE,
#'   n_cores = 4
#' )
#' }
#'
#' @export read.spec
#' @importFrom parallel makeCluster
#' @importFrom parallel detectCores
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
read.spec <- function(
    data_path, 
    file_type = 'txt', 
    sep = NULL,
    group.index = 1, 
    group.names = 'Group', 
    group.levels = NULL, 
    group_splits = '/|_',  
    cutoff = c(500, 3150), 
    interpolation = FALSE, 
    n_cores = 1) {

  pattern <- paste0("\\.(", paste(file_type, collapse = "|"), ")$")
  filenames <- list.files(data_path, pattern = pattern,  full.names = FALSE,include.dirs = TRUE, recursive = TRUE)
  if(length(grep('Metadata', filenames)) != 0)filenames <- filenames[-grep('Metadata', filenames)]

  full_files <- file.path(data_path, filenames)
  if (n_cores > 1){
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {
      suppressMessages(requireNamespace("data.table", quietly = TRUE))
      suppressMessages(requireNamespace("stringr", quietly = TRUE))
    })
    clusterExport(cl, varlist = c( "read_spectral_file"))
    data_list <- parLapply(cl, full_files, function(x) {
      read_spectral_file(x, sep, cutoff, interpolation)
    })
    stopCluster(cl)
  } else{
    data_list <- lapply( full_files, function(x) {read_spectral_file(x, sep, cutoff, interpolation) })
  }

  # Report reading formats within the directory
  orient_list = unlist(lapply(data_list, function(x) x$orient)) 
  print_orient(orient_list, filenames)
  valid_sample <- orient_list != 'unknown'
  data_list <- data_list[valid_sample]
  filenames <- filenames[valid_sample]

  # Identify batches by wavenumber
  wave_list = lapply(data_list, function(x) x$wavenumber)
  wave_len <- sapply(wave_list, length)
  unique_lens <- unique(wave_len)
  if (length(unique_lens) > 1) {
    warning(sprintf('Different wavenumber lengths detected. Returning %d Ramanome objects.', length(unique_lens)) ) 
  }

  inten_list <- lapply(data_list, function(x) x$intensity)
  Ram_list <- list()
  for (len_value in unique_lens){
    idx <- which(wave_len == len_value)
    sub_waves <- wave_list[idx]
    sub_files <- filenames[idx]
    sub_spec <- inten_list[idx]
    data_mat <- do.call(rbind, sub_spec)
    rep_spec <- unlist(sapply(sub_spec, nrow))
    filenames_rep <- rep(sub_files, times = rep_spec)
    
    waves_rep <- rep(sub_waves, times = rep_spec)
    
    wave_keys <- sapply(waves_rep, function(wv) paste(round(wv, 1), collapse = "_"))
    unique_keys <- unique(wave_keys)
    if(length(unique_keys) > 1){
      tab_key <- table(wave_keys)
      most_common_key <- names(tab_key)[which.max(tab_key)]
      wavenumber <- waves_rep[[which(wave_keys == most_common_key)[1]]]
      batch_idx <- match(wave_keys, unique_keys)
      batch_idx <- paste('Batch',as.numeric(as.factor(batch_idx)),sep='_')
      warning(sprintf('Different batches detected, they were combined into the same matrix by row'),  length(unique_keys)) 
    } else wavenumber <- sub_waves[[1]]

    Ram_object <- build_ramanome_object(data_mat, wavenumber,  file_infor = filenames_rep,
                                      meta.data = if (length(unique_keys) > 1) data.frame(batch = batch_idx) else NULL, 
                                      group.index = group.index, group.names = group.names, group.levels = group.levels, group_splits = group_splits)
    Ram_list[[paste0('Wave_points_',len_value)]] <- Ram_object
  }
  
  return(if(length(unique_lens) > 1) Ram_list else Ram_list[[1]])
}

#' Build a Ramanome Object from Parsed Spectral Data
#'
#' @description
#' Constructs a standardized **Ramanome** object from a spectral intensity matrix 
#' and its corresponding wavenumber vector.  
#'
#' This function enables manually assemble of a `Ramanome` object 
#' from pre‑processed spectral data, meta.data.
#'
#' @param data_matrix A numeric matrix of Raman spectral intensities.  
#'        Rows correspond to spectra (samples), and columns correspond to wavenumbers.
#' @param wavenumber A numeric vector giving the wavenumber values associated with each spectral column.  
#' @param file_infor Optional character vector of source file informations.  
#'        If provided, these are used to derive sample grouping and for naming rows.
#' @param meta.data Optional data frame of metadata.  
#'        If supplied, it will override automatically generated grouping information.
#' @param group.index Integer vector specifying which elements from filename splits to use as grouping variables.
#' @param group.names Character vector giving names to group columns.
#' @param group.levels Optional character vector of predefined factor levels.
#' @param group_splits Character pattern (regular expression) used for splitting filenames (default: `"/|_"`).
#'
#' @return A single `Ramanome` object.
#' @examples
#' \dontrun{
#' # Build a Ramanome object from a supplied spectral matrix
#' raw_mat <- matrix(rnorm(10 * 2000), nrow = 10)   # 10 spectra × 2000 wavenumbers
#' wnum <- seq(500, 3150, length.out = 2000)
#' files <- paste0("Sample_", 1:10)
#' Raman <- build_ramanome_object(raw_mat, wnum, filenames = files)
#' }
#' @export build_ramanome_object

build_ramanome_object <- function(
    data_matrix,
    wavenumber,
    file_infor = NULL,
    meta.data = NULL,
    group.index = 1,
    group.names = "Group",
    group.levels = NULL,
    group_splits = "/|_"
) {
  meta.data_defined <- meta.data
  order_index <- order(wavenumber, decreasing = F)
  wavenumber <- wavenumber[order_index]
  data_matrix <- data_matrix[,order_index]
  len_group <- length(group.index)
  colnames(data_matrix) <- wavenumber

  if (!is.null(file_infor)){
    split_list <- strsplit(file_infor, split = group_splits)
    split_mat <- do.call(rbind, split_list)
    group_df <- as.data.frame(split_mat[, group.index])

    if (length(group.names) != len_group && len_group != 1) {
      warning(sprintf("Number of group indices (%d) and names (%d) differ; auto-renaming applied.", len_group, length(group.names)))
    }
    group_name <- ifelse(length(group.names) == len_group, group.names, paste(group.names[1], 1:len_group, sep = '_'))
    colnames(group_df) <- group_name
    group_df[] <- lapply(group_df, factor)

    group <- do.call(paste, c(group_df,sep='_'))

    if (!all(unique(group) %in% group.levels) && !is.null(group.levels)){
      undefined <- setdiff(unique(group), group.levels)
      warning(sprintf("Undefined groups detected: %s. 'group.levels' ignored.", undefined))
      group.levels <- NULL
    }
    if (is.null(group.levels)) group.levels <- unique(group)
    group <- factor(group, levels = group.levels) 

    meta.data <- data.frame(group_df, group = group, filenames = file_infor)
    rownames(data_matrix) <- strip_common_prefix_suffix(meta.data$filenames)
  }
  if ((!is.null(meta.data_defined)) & (!is.null(file_infor))){
    meta.data <- data.frame(meta.data_defined, meta.data)
  }
  
  Ramanome <- new("Ramanome", datasets = list(raw.data = data_mat), wavenumber = wavenumber, meta.data = meta.data)
  show(Ramanome)
  return(Ramanome)
}
