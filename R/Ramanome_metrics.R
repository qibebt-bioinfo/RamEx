# Definition of Ramanome properties
library(stats)
library(MASS)
library(entropy)
get_file_time <- function(file_list){
  file_creation_times <- apply(file_list, 1,function(file) {  
      file_info <- file.info(file)  
      creation_time <- as.Date(file_info$mtime)  
      return(as.character(creation_time))  })
  return(file_creation_times)}

find_variable_features <- function(matrix, num_features = 500) {
  spec_means <- colMeans(matrix)
  spec_vars <- apply(matrix, 2, var)
  spec_disp <- spec_vars / spec_means
  
  spec_means_z <- scale(spec_means)
  spec_disp_z <- scale(spec_disp)
  
  spec_stats <- data.frame(mean = spec_means, disp = spec_disp, mean_z = spec_means_z, disp_z = spec_disp_z)
  
  spec_stats <- spec_stats[order(spec_stats$disp_z, decreasing = TRUE), ]
  variable_wave <- rownames(spec_stats)[1:min(num_features,ncol(matrix))]
  
  return(variable_wave)
}

#' @export 
Raman_attribute <- function(object, min.width=100, max.range=min.width/2, coverage=0.9,n_bins=20){
  matrix <- get.nearest.dataset(object)
  waves <- object@wavenumber
  resolution <- (max(waves)-min(waves))/length(waves)
  n <- nrow(matrix)
  group <- object@meta.data$group
  properties <- list()
  
  tryCatch({  
    if (is.null(group) | length(unique(group)) == 1) {  
      properties$batch_index <- properties$interpretability <- data.frame(score = 0, p = 0)
      warning('Lack of group labels, so the batch index or interpretability could not been estimated! ')  
    }
  }, warning = function(w) {  
    message("Warining: ", w$message)  
  })
  

  # batch_score
  batch <- object@meta.data$batch
  tryCatch({  
    if (is.null(batch)) {  
      batch <- get_file_time(as.matrix(object@meta.data$filenames, ncol=1)) 
      warning('Lack of batch labels, so they are set as the acquisition date')  
    }
  }, warning = function(w) {  
    message("Warining: ", w$message)  
  })
  peaks_locs <- apply(matrix,1, function(x)bubblefill(x,min_bubble_widths = floor(min.width/resolution))$peaks)
  if(is.list(peaks_locs))total_peaks <- do.call(rbind, peaks_locs, quote = F) else total_peaks <- t(peaks_locs)
  
  total_peaks <- sort(unique(as.vector(total_peaks)))
  gaps <- c(1,which(total_peaks[-1]-total_peaks[-length(total_peaks)] >= 2),length(total_peaks))
  common_peaks <- sapply(2:length(gaps), function(i)total_peaks[gaps[i-1]:gaps[i]])
  num_peaks <- list2DF(lapply(common_peaks, function(x)sum(list2DF(lapply(peaks_locs,function(y)any(y %in% x)))) ))
  common_peaks <- common_peaks[which(num_peaks/n > coverage)]
  peaks_wave <- lapply(common_peaks,function(peak_range)apply(matrix,1,function(y)min(peak_range)-1+which.max(y[min(peak_range):max(peak_range)])))
  peaks_wave <- do.call(cbind, peaks_wave)
  peaks_inten <- lapply(1:nrow(peaks_wave), function(x)matrix[x,peaks_wave[x,]])
  peaks_inten <- do.call(rbind, peaks_inten)
  peaks_info <- cbind(peaks_inten,peaks_wave)
  if (!is.null(group) & length(unique(group)) != 1){
    mano_peaks <- summary(manova(peaks_info ~ group + batch))
    properties$batch_index <- data.frame(score=mano_peaks$stats[1,2]/sum(mano_peaks$stats[1:2,2]), p = mano_peaks$stats[2,6])
  }
  
  # SNR Level
  signal_band <- waves < 1800
  noise_band <- waves > 1800 & waves < 2000 | waves > 2350 & waves < 2700
  S <- apply(matrix[,signal_band],1,max)
  N_mean <- apply(matrix[,noise_band],1,mean)
  N_sd <- apply(matrix[,noise_band],1,sd)
  SNR_level <- 1-N_mean/(N_sd*S)
  SNR_level <- fitdistr(SNR_level, densfun = "normal")
  properties$SNR_level <- data.frame(score=SNR_level$estimate[1], p=SNR_level$sd[1])
  
  
  # Variation
  mean_spec <- colMeans(matrix[,!noise_band])
  diff_matrix <- t(t(matrix[,!noise_band]) - mean_spec) 
  Var_signal <- apply(matrix[,!noise_band],2,sd)
  variation <- fitdistr((Var_signal-mean(N_sd))/Var_signal, densfun = 'normal')
  properties$variation <- data.frame(score = variation$estimate[1], p = variation$sd[1])
  
  # Raman Entropy
  compute_bins <- function(x, n_bins = 20) {  
    breaks <- seq(min(x), max(x), length.out = n_bins + 1)  
    cut_data <- cut(x, breaks = breaks, include.lowest = TRUE)  
    freq_table <- table(cut_data)  
    return(freq_table)  
  } 
  entropy_wave <- apply(peaks_wave,2,function(x)entropy(table(x))/log(n))
  entropy_inten <- apply(peaks_inten,2, function(x) entropy(compute_bins(x, n_bins)))/log(n)
  peaks_entropy_wave <- fitdistr( entropy_wave, densfun = "t")
  peaks_entropy_inten <- fitdistr(as.numeric(entropy_inten), densfun = "normal")  
  properties$peaks_entropy <- data.frame(score = mean(peaks_entropy_wave$estimate[1], peaks_entropy_inten$estimate[1]), 
                                        p = mean(peaks_entropy_wave$sd[1], peaks_entropy_inten$sd[1]))
  
  # interpretability
  if (!is.null(group) & length(unique(group)) != 1){
    variable_wave <- find_variable_features(matrix)
    manova_results <- manova(matrix[,waves %in% variable_wave] ~ group)
    properties$interpretability <- data.frame(score = median(1-abs(manova_results[["residuals"]]/matrix[,waves %in% variable_wave] ))) #, p = summary(manova_results)$stats[1,6]
  }
    
  return(properties)
  }