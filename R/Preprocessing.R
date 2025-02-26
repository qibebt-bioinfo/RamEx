#' Preprocess baseline
#'
#' This function performs baseline pre-processing on the input object.
#'
#' @param object The input object
#' @param order The order of the polynomial used in the baseline calculation (default: 1)
#'
#' @return The modified object with preprocessed baseline data
#' @importFrom hyperSpec spc.fit.poly.below
#' @importFrom hyperSpec spc.fit.poly
#' @export Preprocessing.Baseline.Polyfit
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
Preprocessing.Baseline.Polyfit <- function(object, order = 1) {
  wavenumber <- object@wavenumber
  spc2hs <- new(
    "hyperSpec",
    spc = get.nearest.dataset(object),
    wavelength = wavenumber
  )
  
  data_hyperSpec <- spc2hs
  wave_max <- max(data_hyperSpec@wavelength)
  wave_min <- min(data_hyperSpec@wavelength)
  hyperspec <- data_hyperSpec
  if (wave_max >= 3050 & wave_min <= 600) {
    hyperspec_1 <- hyperspec[, , floor(wave_min) ~ 1790] -
      spc.fit.poly.below(
        hyperspec[, , floor(wave_min) ~ 1790],
        hyperspec[, , floor(wave_min) ~ 1790],
        poly.order = order
      )
    hyperspec_2 <- hyperspec[, , 1790 ~ floor(wave_max)] -
      spc.fit.poly(
        hyperspec[, , c(1790 ~ 2065, 2300 ~ 2633, 2783, floor(wave_max))],
        hyperspec[, , 1790 ~ floor(wave_max)],
        poly.order = 6
      )
    hyperspec_baseline <- cbind(hyperspec_1, hyperspec_2)
    print(paste0("The Ramanome contains ", length(hyperspec_baseline), " spectra"))
  }
  else if (wave_max > 1750 & wave_min < 600) {
    hyperspec_baseline <- hyperspec[, , 550 ~ 1750] -
      spc.fit.poly.below(
        hyperspec[, , 550 ~ 1750],
        hyperspec[, , 550 ~ 1750],
        poly.order = order
      )
    print(paste0("The Ramanome contains ", length(hyperspec_baseline), "spectra"))
  }
  else if (wave_max > 3050 &
           wave_min < 1800 &
           wave_min >
           600) {
    hyperspec_baseline <- hyperspec[, , 1800 ~ 3050] -
      spc.fit.poly(
        hyperspec[, , c(1800 ~ 2065, 2300 ~ 2633, 2783, 3050)],
        hyperspec[, , 1800 ~ 3050],
        poly.order = order
      )
    print(paste0("The Ramanome contains ", length(hyperspec_baseline), "spectra"))
  }
  else if (wave_max < 3050 & wave_min > 600) {
    hyperspec_baseline <- hyperspec[, , floor(wave_min) ~ floor(wave_max)] -
      spc.fit.poly(
        hyperspec[, , c(floor(wave_min) ~ floor(wave_min) + 10, 1800~floor(wave_max))],
        hyperspec[, ,floor(wave_min) ~ floor(wave_max) ],
        poly.order = order
      )}
  else if (wave_max < 1750 & wave_min > 600) {
    print("the spc is too small to baseline")
  }
  else {
    print("something is wrong in your spc data")
  }
  
  data_hyperSpec_baseline <- hyperspec_baseline
  
  data <- data_hyperSpec_baseline$spc[, !duplicated(colnames(data_hyperSpec_baseline$spc))]
  object@datasets$baseline.data <- data
  object@wavenumber <- as.numeric(colnames(data))
  
  return(object)
}


#' Remove baseline by bubble
#'
#' This function performs baseline pre-processing on the input object.
#'
#' @param object Ramanome object
#'
#' @return The baseline corrected spectra with bubble method
#' @export Preprocessing.Baseline.Bubble
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
Preprocessing.Baseline.Bubble <- function(object){
  band_1 <- which(object@wavenumber < 1800)
  band_2 <- which(object@wavenumber > 1800 & object@wavenumber < 2700)
  band_3 <- which(object@wavenumber > 2700)
  data_matrix <- get.nearest.dataset(object)
  data_matrix[,band_1] <- t(apply(data_matrix[,band_1], 1, function(x)bubblefill(x)$raman))
  data_matrix[,band_2] <- polyfit(data_matrix[,band_2], degree = 1, tol = 0.001, rep = 100)$corrected
  data_matrix[,band_3] <- polyfit(data_matrix[,band_3], degree = 1, tol = 0.001, rep = 100)$corrected
  data_matrix[, band_2] <- data_matrix[, band_2] - (data_matrix[,band_2[1]]-data_matrix[,band_2[1]-1])
  data_matrix[, band_3] <- data_matrix[, band_3] - (data_matrix[,band_3[1]]-data_matrix[,band_3[1]-1])
  object@datasets$baseline.data <- data_matrix
  return(object)
}


#' Pre-process Raman data by cutting the wavenumber range
#'
#' This function cuts the wavenumber range of Raman data and updates the Raman object.
#'
#' @param object A Ramanome object.
#' @param from The lower bound of the wavenumber range to cut.
#' @param to The upper bound of the wavenumber range to cut.
#'
#' @return The updated Ramanome object with the wavenumber range cut and the corresponding dataset.
#'
#' @export Preprocessing.Cutoff
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocessing.Normalize(data_baseline, "ch")
#' Preprocessing.Cutoff(data_normalized,550, 1800)
Preprocessing.Cutoff <- function(object, from, to) {
  waves <- object@wavenumber
  pred.data <- get.nearest.dataset(object)
  bands <- waves > from & waves < to
  pred.data <- pred.data[, bands]
  object@wavenumber <- waves[bands]
  object@datasets$cut.data <- pred.data
  return(object)
}

#' Pre-process Raman data by ch normalization method
#'
#' This function applies normalization to the Raman data based on the ch method/max method/specific method/area method and updates the Raman object.
#'
#' @param object A Ramanome object.
#' @param method The normalize method
#' @param wave The interested wavenumber
#' @return The updated Ramanome object with normalized Raman data.
#' @export Preprocessing.Normalize
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(data_smoothed)
#' data_normalized_ch <- Preprocessing.Normalize(data_baseline, "ch")
#' #data_normalized_specific <- Preprocessing.Normalize(data_baseline, "specific")
#' data_normalized_max <- Preprocessing.Normalize(data_baseline, "max")
#' data_normalized_area <- Preprocessing.Normalize(data_baseline, "area")
Preprocessing.Normalize <- function(object, method, wave = NULL){
  
  dataset <- get.nearest.dataset(object)
  if(method=="area"){
    value <- rowSums(dataset)
    object@datasets$normalized.data <- dataset / value
    return(object)
  }
  if(method=="specific"){
    if (is.null(wave)) {
      stop('Error! Please input the interested wavenumber!')
    }
    loc <- which.min(abs(object@wavenumber - wave))
    value <- dataset[, loc]
    object@datasets$normalized.data <- dataset / value
    return(object)
  }
  if(method=="max"){
    value <- apply(dataset, 1, max)
    object@datasets$normalized.data <- dataset / value
    return(object)
  }
  if(method=="ch"){
    range <- object@wavenumber > 2850 & object@wavenumber < 3000
    value <- apply(dataset[, range], 1, max)
    object@datasets$normalized.data <- dataset / value
    return(object)
  }
}

#' Smooth Spectral Data Using Savitzky-Golay Filter
#'
#' This function applies the Savitzky-Golay filter to the spectral data in a Ramanome
#' object to reduce noise. The filter is applied to the central portion of the data,
#' preserving the first and last few columns as specified by the `w` parameter.
#'
#' @param object A Ramanome object containing the spectral data.
#' @param m The order of the polynomial used in the Savitzky-Golay filter. Defaults to 0.
#' @param p The number of points to use in the Savitzky-Golay filter. Defaults to 5.
#' @param w The window size used in the Savitzky-Golay filter. Defaults to 11.
#' @param delta.wav The change in wavelength between each data point. Defaults to 2.
#' @return The updated Ramanome object with the smoothed spectral data.
#' @importFrom prospectr savitzkyGolay
#' @export Preprocessing.Smooth.Sg
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
Preprocessing.Smooth.Sg <- function(object, m = 0, p = 5, w = 11, delta.wav = 2) {
  pred.data <- get.nearest.dataset(object)
  pred.data <- cbind(
    pred.data[, 1:((w - 1) / 2)],
    prospectr::savitzkyGolay(pred.data, m = m, p = p, w = w, delta.wav = delta.wav),
    pred.data[, (ncol(pred.data) - (w - 1) / 2 + 1):ncol(pred.data)]
  )
  object@datasets$smooth.data <- pred.data
  return(object)
}



#' Apply Standard Normal Variate (SNV) Transformation
#'
#' This function applies the Standard Normal Variate (SNV) transformation to the raw
#' spectral data in a Ramanome object. The SNV transformation is a type of preprocessing
#' technique used to normalize the data and reduce the impact of multiplicative scatter.
#' The function also reports the number of spectra and the wavelength range after
#' the transformation.
#'
#' @param object A Ramanome object containing the raw spectral data.
#' @return The updated Ramanome object with the SNV-transformed data.
#' @importFrom prospectr standardNormalVariate
#' @export Preprocessing.Smooth.Snv
#' @examples
#' data(RamEx_data)
#' data_smoothed_snv <- Preprocessing.Smooth.Snv(RamEx_data)

Preprocessing.Smooth.Snv <- function(object) {
  #pred.data <- get.nearest.dataset(object)
  
  object@datasets$raw.data <- prospectr::standardNormalVariate(as.data.frame(object@datasets$raw.data))
  
  writeLines(paste("after snv there is ", length(object), " spec data", sep = ""))
  min_wave <- min(object@wavenumber)
  max_wave <- max(object@wavenumber)
  writeLines(paste("wave range is ", min_wave, "-", max_wave, sep = ""))
  return(object)
}


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
#' @export Preprocessing.Background.Remove

Preprocessing.Background.Remove <- function(object, cell.index, cal_mean = FALSE) {
  
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


#' Convert Image Matrix to Vector Form
#'
#' This function converts a 2D image matrix into a vectorized form by sliding a kernel window
#' over the image and extracting blocks. Each block is reshaped into a column vector.
#'
#' @param image A numeric matrix representing the input image
#' @param kernel_height Height of the kernel window
#' @param kernel_width Width of the kernel window
#' @param device Device to use for computation - either "CPU" or "GPU"
#'
#' @return A matrix where each column is a vectorized block from the input image.
#'         For CPU device, returns a regular matrix.
#'         For GPU device, returns a gpuMatrix.
matrix2vector <- function(image, kernel_height, kernel_width, device) {
  image_height <- nrow(image)
  image_width <- ncol(image)
  output_height <- image_height - kernel_height + 1
  output_width <- image_width - kernel_width + 1
  if(device == "CPU"){
    cols <- matrix(0, nrow = kernel_height * kernel_width, ncol = output_height * output_width)
  }
  if(device == "GPU"){
    cols <- gpuMatrix(0, nrow = kernel_height * kernel_width, ncol = output_height * output_width)
  }
  idx <- 1
  for (i in 1:output_height) {
    for (j in 1:output_width) {
      
      block <- image[(i):((i + kernel_height - 1)), (j):((j + kernel_width - 1))]
      cols[, idx] <- as.vector(block)
      idx <- idx + 1
    }
  }
  return(cols)
}






#' Pre spike the spec. with Rcpp
#'
#' This function performs Pre spike function
#'
#' @param all_data The rammannome class
#' @param device Parallel device, default is "CPU", can be set to "GPU"
#' @return A rammanome object


gpupre_spike_matrix <- function(all_data, device = "CPU") {
  all_spc <- as.matrix(all_data[,-1])
  
  if (nrow(all_spc) < 3 || ncol(all_spc) < 11) {
    stop("all_spc must be at least 3 rows and 11 columns")
  }
  
  kernel_ori <- matrix(c(-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,33, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1),nrow = 3, byrow = TRUE)
  
  
  kernel_matrix <- as.vector(t(kernel_ori))
  
  
  if(device == "GPU"){
    vector_result <- matrix2vector(all_spc, 3, 11, device)
    gpu_kernel_matrix <- as.gpuVector(kernel_matrix)
    conv_result <-  gpu_kernel_matrix %*% vector_result
  }
  if(device == "CPU"){
    vector_result <- matrix2vector(all_spc, 3, 11, device)
    conv_result <-  kernel_matrix %*% vector_result
  }
  
  
  all_spc_1 <- matrix(conv_result, nrow = nrow(all_spc) - 2, ncol = ncol(all_spc) - 10, byrow = TRUE)
  
  slope_inds_2 <- which(all_spc_1 > 10 * apply(all_spc, 1, max), arr.ind = TRUE)
  
  slope_inds <- slope_inds_2
  
  spc_new <- all_data
  wavenumber <- as.numeric(as.character(colnames(spc_new[,-1])))
  
  for (ind in unique(slope_inds[,1])) {
    spike_pos <- slope_inds[which(slope_inds[,1] == ind), 2]
    if (ind <= 3) {
      inds_new <- c((ind + 1):(ind + 3))
    } else if (ind > 3 & ind <= (nrow(all_data) - 3)) {
      inds_new <- c((ind - 3):(ind - 1), (ind + 1):(ind + 3))
    } else if (ind > (nrow(all_data) - 3)) {
      inds_new <- c((ind - 3):(ind - 1))
    }
    
    for (i in spike_pos) {
      if (i == 1) {
        spc_new[ind, 2:4] <- colMeans(spc_new[inds_new, 2:4])
      } else if (i >= length(wavenumber)) {
        spc_new[ind, (i - 2):i] <- colMeans(spc_new[inds_new, (i - 2):i])
      } else {
        spc_new[ind, i:(i + 2)] <- colMeans(spc_new[inds_new, i:(i + 2)])
      }
    }
  }
  
  return(spc_new)
}

#' Pre spike the spec. with Rcpp
#'
#' This function performs Pre spike function
#'
#' @param object The rammanome class
#' @param device Parallel device, default is null for CPU, can be set to "GPU"
#' @return A rammanome object
#'
#' @export Preprocessing.Background.Spike
#' @examples
#' data(RamEx_data)
#' data_spike <- Preprocessing.Background.Spike(RamEx_data,"CPU")
Preprocessing.Background.Spike <- function(object, device){
  data <- get.nearest.dataset(object)
  pre.data <- gpupre_spike_matrix(data, device)
  object@datasets$spike.data <- pre.data
  return(object)
}