#' One step spectral preprocessing
#'
#' This function performs one-step spectral preprocessing on the input Ramanome, including smoothing(SG), baseline correction(Polyfit), normalization(ch), and quality control(ICOD).
#' If specific preprocessing methods are needed, please define/run them step by step.
#' 
#' @param object Ramanome object
#'
#' @return The Ramanome object with preprocessed spectral data
#' @export Preprocessing.OneStep
#' @examples
#' data(RamEx_data)
#' data_preprocessed <- Preprocessing.OneStep(RamEx_data) 
Preprocessing.OneStep <- function(object) {
     data_normalized <- Preprocessing.Smooth.Sg(object) %>% Preprocessing.Baseline.Polyfit %>% Preprocessing.Normalize(., "ch")
     qc_icod <- Qualitycontrol.ICOD(data_normalized)
     data_cleaned <- data_normalized[qc_icod$quality,]
     return(data_cleaned)}


#' Remove baseline (low frequency noises) by polynomial fitting
#'
#' Partitioned polynomial fitting was performed using a two-region fitting approach. A first-order (alterable) polynomial was fitted to the 500–1800 cm⁻¹ region, while a sixth-order (alterable) polynomial was applied to the 1800–3050 cm⁻¹ region.
#'
#' @param object Ramanome object
#' @param order_1 The order of the polynomial used in the baseline calculation for the 500–1800 cm⁻¹ region (default: 1)
#' @param order_2 The order of the polynomial used in the baseline calculation for the 1800–3050 cm⁻¹ region (default: 6)
#'
#' @return The Ramanome object with baseline corrected spectra, named as "baseline.data" in the 'datasets' slot
#' @export Preprocessing.Baseline.Polyfit
#' @examples
#' data(RamEx_data)
#' data_baseline <- Preprocessing.Baseline.Polyfit(RamEx_data)

Preprocessing.Baseline.Polyfit <- function(object, order_1 = 1, order_2 = 6) {
  wavenumber <- object@wavenumber
  spc_data <- get.nearest.dataset(object)
  
  wave_max <- max(wavenumber)
  wave_min <- min(wavenumber)
  
  if (wave_max >= 2600 & wave_min <= 1100) {
    idx_1 <- which(wavenumber >= floor(wave_min) & wavenumber <= 1790)
    baseline1 <- spc.fit.poly.below(spc_data[, idx_1, drop = FALSE], spc_data[, idx_1, drop = FALSE], poly.order = order_1)
    corrected1 <- spc_data[, idx_1, drop = FALSE] - baseline1
    
    idx_2 <- which(wavenumber >= 1790 & wavenumber <= wave_max)
    support_idx <- which(wavenumber >= 1790 & wavenumber <= 2065 | 
                           wavenumber >= 2300 & wavenumber <= 2633 | 
                           wavenumber == wavenumber[which.min(abs(wavenumber - 2783))] | 
                           wavenumber == wave_max)
    baseline2 <- spc.fit.poly(spc_data[, support_idx, drop = FALSE], spc_data[, idx_2, drop = FALSE], poly.order = order_2)
    corrected2 <- spc_data[, idx_2, drop = FALSE] - baseline2
    corrected <- cbind(corrected1, corrected2)
    message(paste0("The Ramanome contains ", nrow(corrected), " spectra"))
    
  } else if (wave_min <= 1100) {
    baseline1 <- spc.fit.poly.below(spc_data, spc_data, poly.order = order_1)
    corrected <- spc_data - baseline
  } else if (wave_min > 1100 & wave_max > 2600) {
    support_idx <- which(wavenumber == wave_min |
                          wavenumber >= 1790 & wavenumber <= 2065 | 
                          wavenumber >= 2300 & wavenumber <= 2633 | 
                          wavenumber == wavenumber[which.min(abs(wavenumber - 2783))] | 
                          wavenumber == wave_max)
    baseline <- spc.fit.poly(spc_data[, support_idx, drop = FALSE], spc_data, poly.order = order_2)
    corrected <- spc_data - baseline
    message(paste0("The Ramanome contains ", nrow(corrected), " spectra"))
    
  } else {
    baseline1 <- spc.fit.poly.below(spc_data, spc_data, poly.order = order_1)
    corrected <- spc_data - baseline
  }
  
  corrected <- corrected[, !duplicated(colnames(corrected), fromLast = FALSE), drop = FALSE]
  object@datasets$baseline.data <- corrected
  object@wavenumber <- as.numeric(colnames(corrected))
  
  return(object)
}


vanderMonde <- function(x, poly.order) {
  sapply(0:poly.order, function(i) x^i)
}

spc.fit.poly <- function(fit.to, apply.to = NULL, poly.order = 1, offset.wl = !is.null(apply.to)) {
  x <- as.numeric(colnames(fit.to))
  if (offset.wl) {
    minx <- min(x)
    x_adj <- x - minx
  } else {
    minx <- 0
    x_adj <- x
  }
  
  vdm <- vanderMonde(x_adj, poly.order)
  
  p <- t(apply(fit.to, 1, function(y, vdm) {
    valid <- !is.na(y)
    if (sum(valid) < (poly.order + 1)) {
      rep(NA, poly.order + 1)
    } else {
      qr.solve(vdm[valid, , drop = FALSE], y[valid])
    }
  }, vdm))
  
  if (is.null(apply.to)) {
    colnames(p) <- paste0("(x-minx)^", 0:poly.order)
    return(list(spc = p, wavelength = 0:poly.order, min.x = minx))
  } else {
    wl_new <- as.numeric(colnames(apply.to)) - minx
    vdm_new <- vanderMonde(wl_new, poly.order)
    baseline <- t(apply(p, 1, function(coefs) as.numeric(vdm_new %*% coefs)))
    return(baseline)
  }
}

spc.fit.poly.below <- function(fit.to, apply.to = fit.to, poly.order = 1,
                               npts.min = max(round(ncol(fit.to) * 0.05), 3 * (poly.order + 1)),
                               noise = 0, offset.wl = FALSE, max.iter = ncol(fit.to),
                               stop.on.increase = FALSE, debuglevel = 0) {
  x <- as.numeric(colnames(fit.to))
  if (offset.wl) {
    minx <- min(x)
    x <- x - minx
  } else {
    minx <- 0
  }
  vdm <- vanderMonde(x, poly.order)
  n_spectra <- nrow(fit.to)
  pcoef <- matrix(NA, nrow = n_spectra, ncol = poly.order + 1)
  
  for(i in 1:n_spectra) {
    y <- fit.to[i, ]
    valid <- !is.na(y)
    use <- valid
    use_old <- rep(FALSE, length(y))
    iter <- 0
    repeat {
      iter <- iter + 1
      if(sum(use) < (poly.order + 1)) {
        warning("Spectrum ", i, ": not enough valid points for fitting.")
        break
      }
      coef_i <- qr.solve(vdm[use, , drop = FALSE], y[use])
      pcoef[i,] <- coef_i
      bl <- as.vector(vdm %*% coef_i)
      use_old <- use
      use <- (y < (bl + noise)) & valid
      if(sum(use, na.rm = TRUE) < npts.min || all(use == use_old)) break
      if(stop.on.increase && sum(use, na.rm = TRUE) > sum(use_old, na.rm = TRUE)) break
      if(iter >= max.iter) break
    }
  }
  
  wl <- as.numeric(colnames(apply.to)) - minx
  vdm_new <- vanderMonde(wl, poly.order)
  baseline <- t(apply(pcoef, 1, function(coefs) as.numeric(vdm_new %*% coefs)))
  return(baseline)
  
}


#' Remove baseline (low frequency noises) by bubble
#'
#' Bubble method is more suitable for spectra with complex background and uneven baseline. Bubble filling was performed on the 500–1800 cm⁻¹ region, while polynomial fitting was applied to the 1800–3050 cm⁻¹ region.
#'
#' @param object Ramanome object
#' @param min_bubble_widths The minimum width of the bubbles.
#'
#' @return The Ramanome object with baseline corrected spectra, named as "baseline.data" in the 'datasets' slot
#' @export Preprocessing.Baseline.Bubble
#' @examples
#' data(RamEx_data)
#' data_baseline_bubble <- Preprocessing.Baseline.Bubble(RamEx_data)
Preprocessing.Baseline.Bubble <- function(object,min_bubble_widths = 50){
  band_1 <- which(object@wavenumber < 1800)
  band_2 <- which(object@wavenumber > 1800 & object@wavenumber < 2700)
  band_3 <- which(object@wavenumber > 2700)
  data_matrix <- get.nearest.dataset(object)
  data_matrix[,band_1] <- t(apply(data_matrix[,band_1], 1, function(x)bubblefill(x,min_bubble_widths)$raman))
  data_matrix[,band_2] <- polyfit(data_matrix[,band_2], degree = 1, tol = 0.001, rep = 100)$corrected
  data_matrix[,band_3] <- polyfit(data_matrix[,band_3], degree = 1, tol = 0.001, rep = 100)$corrected
  data_matrix[, band_2] <- data_matrix[, band_2] - (data_matrix[,band_2[1]]-data_matrix[,band_2[1]-1])
  data_matrix[, band_3] <- data_matrix[, band_3] - (data_matrix[,band_3[1]]-data_matrix[,band_3[1]-1])
  object@datasets$baseline.data <- data_matrix
  return(object)
}


#' Spectral truncation
#' 
#' Cutoff the spectra to the specific wavenumber range.
#'
#' @param object A Ramanome object.
#' @param from The lower bound of the wavenumber range to cut.
#' @param to The upper bound of the wavenumber range to cut.
#'
#' @return A Ramanome object with truncated spectra, named as "cut.data" in the 'datasets' slot
#'
#' @export Preprocessing.Cutoff
#' @examples
#' data(RamEx_data)
#' Preprocessing.Cutoff(RamEx_data,550, 1800)
Preprocessing.Cutoff <- function(object, from, to) {
  waves <- object@wavenumber
  pred.data <- get.nearest.dataset(object)
  bands <- waves > from & waves < to
  pred.data <- pred.data[, bands]
  object@wavenumber <- waves[bands]
  object@datasets$cut.data <- pred.data
  return(object)
}

#' Normalize the spectra matrix
#'
#' @param object A Ramanome object.
#' @param method The normalize method. "area" for area normalization, "specific" for specific wavenumber normalization, "max" for maximum normalization, "ch" for continuum normalization.
#' @param wave The interested wavenumber for specific normalization.
#' 
#' @return A Ramanome object with normalized spectra, named as "normalized.data" in the 'datasets' slot
#' @export Preprocessing.Normalize
#' 
#' @examples
#' data(RamEx_data)
#' data_normalized_ch <- Preprocessing.Normalize(RamEx_data, "ch")
#' #data_normalized_specific <- Preprocessing.Normalize(RamEx_data, "specific", 1650)
#' data_normalized_max <- Preprocessing.Normalize(RamEx_data, "max")
#' data_normalized_area <- Preprocessing.Normalize(RamEx_data, "area")
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

#' Smooth Spectral Data Using Savitzky-Golay Filter to reduce high-frequency noises
#'
#' The Savitzky-Golay filter is a digital smoothing technique that preserves the shape and height of spectral peaks by fitting a polynomial to a moving window of data points.
#'
#' @param object A Ramanome object.
#' @param m Whether perform differentiation. Defaults to 0.
#' @param p The order of the polynomial used in the Savitzky-Golay filter. Defaults to 5.
#' @param w The window size used for polynomial fitting. Defaults to 11.
#' @return A Ramanome object with smoothed spectra, named as "smooth.data" in the 'datasets' slot
#' @importFrom prospectr savitzkyGolay
#' @export Preprocessing.Smooth.Sg
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocessing.Smooth.Sg(RamEx_data)
Preprocessing.Smooth.Sg <- function(object, m = 0, p = 5, w = 11) {
  pred.data <- get.nearest.dataset(object)
  pred.data <- cbind(
    pred.data[, 1:((w - 1) / 2)],
    prospectr::savitzkyGolay(pred.data, m = m, p = p, w = w, delta.wav = 2),
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
#' @param object A Ramanome object.
#' @return A Ramanome object with SNV-transformed spectra, named as "smooth.data" in the 'datasets' slot
#' @importFrom prospectr standardNormalVariate
#' @export Preprocessing.Smooth.Snv
#' @examples
#' data(RamEx_data)
#' data_smoothed_snv <- Preprocessing.Smooth.Snv(RamEx_data)

Preprocessing.Smooth.Snv <- function(object) {
  pred.data <- get.nearest.dataset(object)
  
  object@datasets$smooth.data <- prospectr::standardNormalVariate(as.data.frame(pred.data))
  
  writeLines(paste("after snv there is ", length(object), " spec data", sep = ""))
  min_wave <- min(object@wavenumber)
  max_wave <- max(object@wavenumber)
  writeLines(paste("wave range is ", min_wave, "-", max_wave, sep = ""))
  return(object)
}


#' Subtract background spectra from the cell spectra
#'
#' This function performs background subtraction on a Ramanome object by subtracting the background spectra from the dataset.
#'
#' @param object A Ramanome object.
#' @param cell.index The index of the cell component in the filenames, where the background spectra are recorded as 'bg' or 'BG'.
#' @param cal_mean Logical value indicating whether to calculate the mean spectra with in a cell when a cell contains multiple spectra.
#'
#' @return A Ramanome object with background-subtracted spectra, named as "sub_data" in the 'datasets' slot
#' @importFrom stringr str_extract
#' @importFrom stringr str_split_i
#' @importFrom stringr str_detect
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
    data.mean <- aggregate(data, by = list(group, cell), mean)
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
#' @noRd 
matrix2vector <- function(image, kernel_height, kernel_width, device) {
  image_height <- nrow(image)
  image_width <- ncol(image)
  padding_width <- floor(kernel_width/2 )
  padding_height <- floor(kernel_height/2 )
  
  image <- rbind(image[rep(1,padding_height),],image, image[rep(image_height,padding_height),])
  image <- cbind(image[,rep(1,padding_width)],image, image[,rep(image_width, padding_width)])
  if(device == "CPU"){
    cols <- matrix(0, nrow = kernel_height * kernel_width, ncol = image_height * image_width)
  }
  if(device == "GPU"){
    cols <- gpuMatrix(0, nrow = kernel_height * kernel_width, ncol = image_height * image_width)
  }
  idx <- 1
  for (i in 1:image_height) {
    for (j in 1:image_width) {
      
      block <- image[(i):((i + kernel_height - 1)), (j):((j + kernel_width - 1))]
      cols[, idx] <- as.vector(block)
      idx <- idx + 1
    }
  }
  return(cols)
}






#' Remove cosmic rays from the spectra
#'
#' recognize the cosmic rays by convolution with a kernel and replace them with the average of the neighboring pixels
#'
#' @param all_data The rammannome class
#' @param device Parallel device, default is null for CPU, can be set to "GPU"
#' @return A rammanome object
#' @noRd


gpupre_spike_matrix <- function(all_data, device = NULL, sharpen=33) {
  kernel_ori <- matrix(c(-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,sharpen, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1),nrow = 3, byrow = TRUE)
  
  kernel_matrix <- as.vector(t(kernel_ori))
  
  if(device == "GPU"){
    vector_result <- matrix2vector(all_data, 3, 11, device)
    gpu_kernel_matrix <- as.gpuVector(kernel_matrix)
    conv_result <-  gpu_kernel_matrix %*% vector_result
  }
  if(device == "CPU"){
    vector_result <- matrix2vector(all_data, 3, 11, device)
    conv_result <-  kernel_matrix %*% vector_result
  }
  
  all_spc_1 <- matrix(conv_result, nrow = nrow(all_data) , ncol = ncol(all_data) , byrow = TRUE)
  
  slope_inds <- which(all_spc_1 > 10 * apply(all_data, 1, max), arr.ind = TRUE)
  
  spc_new <- all_data
  
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
        spc_new[ind, 1:3] <- colMeans(spc_new[inds_new, 1:3])
      } else if (i >= ncol(all_data)) {
        spc_new[ind, c(i - 2,i-1,i)] <- colMeans(spc_new[inds_new, c(i - 2,i-1,i)])
      } else {
        spc_new[ind, c(i-1,i,i+1)] <- colMeans(spc_new[inds_new, c(i-1,i,i+1)])
      }
    }
  }
  
  return(spc_new)
}

#' Remove cosmic rays from the spectra
#'
#' recognize the cosmic rays by a sharp convolution kernel and replace them with the average of the neighboring pixels
#'
#' @param object The rammanome class
#' @param device Parallel device, default is null for CPU, can be set to "GPU"  
#' @param sharpen The sharpening factor of the convolution kernel, default is 33
#' 
#' @return The Ramanome object with cosmic rays removed spectra, named as "spike.data" in the 'datasets' slot
#'
#' @export Preprocessing.Spike
#' @examples
#' data(RamEx_data)
#' data_spike <- Preprocessing.Spike(RamEx_data,"CPU")
#' @note For Windows users: GPU acceleration requires proper installation of OpenCL SDK and GPU drivers.
#' If you encounter issues with gpuR installation, please refer to the package documentation
#' or use CPU mode instead.
Preprocessing.Spike <- function(object, device='CPU', sharpen=33){
  data <- get.nearest.dataset(object)
  if (device == "GPU") {
    if (!requireNamespace("gpuR", quietly = TRUE)) {
        stop("gpuR package is required for GPU acceleration. Please install it first.")
    }
    if (.Platform$OS.type == "windows") {
        warning("GPU acceleration on Windows may require additional setup. Please ensure OpenCL is properly configured.")
    }
}
  pre.data <- gpupre_spike_matrix(data, device, sharpen)
  object@datasets$spike.data <- pre.data
  return(object)
}