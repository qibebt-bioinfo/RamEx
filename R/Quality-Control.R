#' Polynomial fitting for spectral baseline correction
#'
#' This function performs polynomial fitting to correct the baseline of spectral data.
#' It iteratively fits a polynomial to the data, ensuring that the fitted values
#' do not exceed the original data values.
#'
#' @param spectra A matrix of spectral data where each row represents a spectrum.
#' @param t A numeric vector of x-values for the polynomial fitting. If missing or FALSE,
#' a sequence from 1 to the number of columns in `spectra` is used.
#' @param degree The degree of the polynomial to fit.
#' @param tol The tolerance for convergence of the iterative fitting process.
#' @param rep The maximum number of iterations to perform.
#' @return A list containing the baseline and the corrected spectra.
polyfit <- function (spectra, t, degree = 4, tol = 0.001, rep = 100)
{
  dimnames(spectra) <- NULL
  np <- dim(spectra)
  baseline <- matrix(0, np[1], np[2])
  if (missing(t) || (t == FALSE))
    t <- 1:np[2]
  polx <- cbind(1/sqrt(np[2]), stats::poly(t, degree = degree))
  for (i in 1:np[1]) {
    ywork <- yold <- yorig <- spectra[i, ]
    nrep <- 0
    repeat {
      nrep <- nrep + 1
      ypred <- polx %*% crossprod(polx, yold)
      ywork <- pmin(yorig, ypred)
      crit <- sum(abs((ywork - yold)/yold), na.rm = TRUE)
      if (crit < tol || nrep > rep)
        break
      yold <- ywork
    }
    baseline[i, ] <- ypred
  }
  list(baseline = baseline, corrected = spectra - baseline)
}


#' Grow a Bubble Chart
#'
#' This function generates a bubble chart representation of a spectrum by expanding a bubble
#' until it touches the spectrum. The bubble can be aligned to the left, right, or centered.
#'
#' @param spectrum A numeric vector representing the spectrum.
#' @param alignment A character string specifying the alignment of the bubble, which can be "left", "right", or "center".
#' @return A list containing the bubble values and the relative touching point.

grow_bubble <- function(spectrum, alignment = "center") {
  xaxis <- 0:(length(spectrum)-1)
  # Adjusting bubble parameter based on alignment
  if (alignment == "left") {
    # half bubble right
    width <- 2 * length(spectrum)
    middle <- 1
  } else if (alignment == "right") {
    # half bubble left
    width <- 2 * length(spectrum)
    middle <- length(spectrum)
  } else {
    # Centered bubble
    width <- length(spectrum)
    middle <- length(spectrum) / 2
  }

  squared_arc <- (width / 2) ^ 2 - (xaxis - middle) ^ 2  # squared half circle
  # squared_arc[squared_arc < 0] <- 0
  bubble <- sqrt(squared_arc) - width
  # find new intersection
  touching_point <- which.min(spectrum - bubble)

  # grow bubble until touching
  bubble <- bubble + min(spectrum - bubble)

  return(list(bubble=bubble, relative_touching_point=touching_point))
}


#' Keep the largest values between baseline and bubble
#'
#' This function compares the baseline and bubble values and keeps the larger ones.
#'
#' @param baseline A vector of baseline values.
#' @param bubble A vector of bubble values.
#' @return A vector of the larger values

keep_largest <- function(baseline, bubble) {
  for (i in 1:length(baseline)) {
    if (baseline[i] < bubble[i]) {
      baseline[i] <- bubble[i]
    }
  }
  return(baseline)
}


#' Process a spectrum with bubble loop
#'
#' This function iteratively applies the bubble growing algorithm to a spectrum, adjusting the baseline
#' and identifying the bounds of each bubble.
#'
#' @param spectrum A numeric vector representing the spectrum.
#' @param baseline A numeric vector representing the baseline of the spectrum.
#' @param min_bubble_widths A numeric value or vector specifying the minimum width of the bubbles.
#' @return A list containing the adjusted baseline and the bounds of each bubble.
bubbleloop <- function(spectrum, baseline, min_bubble_widths) {
  range_cue <- list(c(1, length(spectrum)))
  bounds <- list()
  i <- 1
  while (i <= length(range_cue)) {
    left_bound <- range_cue[[i]][1]
    right_bound <- range_cue[[i]][2]

    i <- i + 1

    if (right_bound == left_bound) {
      next
    }

    if (is.numeric(min_bubble_widths)) {
      min_bubble_width <- min_bubble_widths
    } else {
      min_bubble_width <- min_bubble_widths[(left_bound + right_bound) %/% 2]
    }
    if (left_bound == 1 & right_bound != length(spectrum)) {
      alignment <- "left"
    } else if (left_bound != 1 & right_bound == length(spectrum)) {
      alignment <- "right"
    } else {
      if ((right_bound - left_bound) < min_bubble_width) {
        next
      }
      alignment <- "center"
    }
    bubble <- grow_bubble(spectrum[left_bound:right_bound], alignment)$bubble
    relative_touching_point <- grow_bubble(spectrum[(left_bound:right_bound)], alignment)$relative_touching_point
    touching_point <- relative_touching_point + left_bound -1
    bounds <- append(bounds, c(left_bound,right_bound))
    baseline[left_bound:right_bound] <- keep_largest(baseline[left_bound:right_bound], bubble)
    if (touching_point == left_bound) {
      range_cue <- c(range_cue, list(c(touching_point+1, right_bound)))
    } else if (touching_point == right_bound) {
      range_cue <- c(range_cue, list(c(left_bound, touching_point-1)))
    } else {
      range_cue <- c(range_cue, list(c(left_bound, touching_point - 1)), list(c(touching_point, right_bound)))
    }
  }
  return(list(baseline=baseline, bounds=range_cue))
}


#' Fill bubbles in a spectrum to correct baseline
#'
#' This function applies a bubble-filling algorithm to a spectrum to correct the baseline.
#' It uses the `bubbleloop` function to identify and adjust the baseline, and then
#' identifies the peaks in the spectrum after baseline correction.
#'
#' @param spectrum A numeric vector representing the spectrum.
#' @param min_bubble_widths A numeric value or vector specifying the minimum width of the bubbles.
#'   If a vector, it should be of the same length as the spectrum, specifying the minimum width for each point.
#' @return A list containing the corrected spectrum (raman), the peak locations (peaks), and the band boundaries (bands).
#' @importFrom prospectr savitzkyGolay
bubblefill <- function(spectrum, min_bubble_widths = 50) {
  if(is.null(names(spectrum))) xaxis <- as.numeric(colnames(spectrum)) else xaxis <- as.numeric(names(spectrum))
  spectrum_ <- spectrum
  # smin <- min(spectrum_)
  # spectrum_ <- spectrum_ - smin
  # scale <- max(spectrum_) / length(spectrum)
  # spectrum_ <- spectrum_ / scale
  baseline <- rep(0, length(spectrum))
  baseline <- bubbleloop(spectrum_, baseline, min_bubble_widths)
  # baseline$baseline <- baseline$baseline * scale + poly_fit + smin
  if (!is.numeric(min_bubble_widths)) {
    filter_width <- max(pmin(min_bubble_widths))
  } else {
    filter_width <- max(min_bubble_widths)
  }
  # baseline <- prospectr::savitzkyGolay(baseline, m = 0, p = 3, w = round(2 * (filter_width %/% 4) + 3), delta.wav = 0)
  raman <- spectrum - baseline$baseline
  bands <- list2DF(baseline$bounds)
  bands <- bands[,which(bands[2,]-bands[1,] < min_bubble_widths)]
  bands <- bands[,which(bands[2,]-bands[1,] > 3)]
  # peaks <- unique(as.vector(as.matrix(bands[1,])))
  peaks <- apply(bands,2,function(x)x[1]+which.max(spectrum[x[1]:x[2]])-1)
  return(list(raman=raman, peaks=peaks, bands=bands))
}


#' Custom Convolution Function
#'
#' This function performs a convolution operation on a numeric vector using a specified kernel.
#' The convolution is performed using the filter type, which applies the kernel to the vector.
#'
#' @param vec A numeric vector to be convolved.
#' @param kernel A numeric vector representing the convolution kernel.
#' @return A numeric vector resulting from the convolution operation.
#' @importFrom stats convolve


convolve_custom <- function(vec, kernel) {
  result <- convolve(as.numeric(vec), rev(kernel), type = "filter")
  return(result)
}


#' Estimate outliers using the Minimum Covariance Determinant (MCD) method
#'
#' This function implements the Minimum Covariance Determinant method to estimate
#' outliers in a dataset. It calculates the MCD estimate of scatter and location,
#' and identifies observations that are likely to be outliers based on their distance
#' from the centroid.
#'
#' @param x The input dataset.
#' @param index_good A vector indicating the indices of good observations to use for the estimation.
#' @param h The fraction of the data to use for computing the MCD (between 0 and 1). Defaults to 0.75.
#' @param alpha The significance level for the MCD. Defaults to 0.01.
#' @param na.rm A logical value indicating whether to remove NA values. Defaults to TRUE.
#' @return A list containing the following components:
#'   - MaxDist: The maximum distance threshold for outliers.
#'   - center: The estimated centroid of the data.
#'   - outliers_pos: The positions of the outliers.
#'   - outliers_val: The values of the outliers.
#'   - dist_from_center: The distances of each observation from the centroid.
#' @importFrom stats qchisq
#' @importFrom stats mahalanobis


outliers_mcdEst <- function(x,index_good,
                            h = .75, # fraction of data we wanna keep
                            # to compute the MCD (between 0 and 1)
                            alpha = .01,
                            na.rm = TRUE){

  if (na.rm == TRUE ) {
    data <- na.omit(x)
  } else {data <- x}

  for (i in seq_len(ncol(data))){
    if(inherits(data[,i],c("numeric","integer")) == FALSE)
      stop("Data are neither numeric nor integer")
  }

  #Creating covariance matrix for Minimum Covariance Determinant
  # by default, use the "best" method = exhaustive method
  print(nrow(data))
  # output <- fastMCD(data, h=0)
  output <- covFastMCD(data[index_good,], alpha = 0.6, m=10, l=1, delta = 0.05)
  print('FastMCD finished!')
  # output <- cov.mcd(data,cor = FALSE,quantile.used = nrow(data)*h)

  cutoff <- (qchisq(p = 1-alpha, df = ncol(data)))
  # cor = FALSE to avoid useless output(correlation matrix)

  #Distances from centroid for each matrix
  cat('Calculating the mahalanobis distance \n')
  dist <- mahalanobis(data,output$center,output$cov, tol=-1) # distance

  #Detecting outliers
  names_outliers <- which(dist > cutoff)
  coordinates <- data.frame(
    matrix(NA, nrow = length(names_outliers), ncol = ncol(data))
  )
  for (k in seq_len(ncol(data))){
    coordinates[,k] <- data[,k][dist > cutoff]
  }


  # print results
  meth <- "Minimum Covariance Determinant estimator"

  # Return results in list()
  invisible(
    list(MaxDist = cutoff,
         center = output$center,
         outliers_pos = names_outliers,
         outliers_val=coordinates,
         dist_from_center = dist)
  )
}


#' Fast Minimum Covariance Determinant Estimation
#'
#' This function implements the Fast Minimum Covariance Determinant (MCD) method to estimate
#' the covariance matrix and identify potential outliers in the data. The MCD method is
#' robust to outliers and can be used for multivariate data.
#'
#' @param x The input dataset, a matrix or data frame.
#' @param alpha The significance level for the MCD estimation.
#' @param m The number of initial subsets to consider.
#' @param l The number of final subsets to consider.
#' @param delta The quantile for the reweighting step.
#' @return A list containing the raw center, raw covariance matrix, best subset indices,
#' weights for reweighting, and the reweighted center and covariance matrix.
#' @importFrom stats qchisq
#' @importFrom stats pgamma
#' @importFrom stats mahalanobis
#' @importFrom robustbase h.alpha.n
#' @importFrom robustbase .MCDcnp2
covFastMCD <- function(x, alpha, m, l, delta) {
  # Convert input to data frame
  x <- data.frame(as.matrix(x, nrow = nrow(x), ncol = ncol(x)))
  p <- ncol(x)  # Number of variables
  n <- nrow(x)  # Number of observations
  h <- h.alpha.n(alpha, n, p)  # Compute subset size h based on alpha, n, and p

  determinant_vec <- rep(0, m)  # Initialize determinant vector
  q_alpha <- qchisq(alpha, df = p)
  q_delta <- qchisq(1 - delta, df = p)
  c_alpha <- alpha / pgamma(q_alpha / 2, shape = p / 2 + 1, scale = 1)
  c_delta <- (1 - delta) / pgamma(q_delta / 2, shape = p / 2 + 1, scale = 1)
  cnp <- .MCDcnp2(p, n, alpha)

  # Check conditions for input data
  if (n <= p) {
    stop("Error: n <= p. The sample size is too small and MCD cannot be performed!")
  }
  if (n == p + 1) {
    stop("Error: n == p + 1. The sample size is too small for MCD.")
  }
  if (n < 2 * p) {
    warning("Warning: n < 2 * p. The sample size might be too small.")
  }
  if (h > n) {
    stop("Error: h value is too large: h > n. Choose a lower alpha.")
  }
  if (h == n) {
    # When h equals n, the MCD location estimator is the mean of the entire dataset
    center_mcd <- apply(x, 2, mean)
    scatter_mcd <- (h - 1) / h * var(x)
    return(list(centerh = center_mcd, scatterh = scatter_mcd))
  }

  # Initialize lists to store subsets and their statistics
  center_0 <- list()
  scatter_0 <- list()
  subset_0 <- list()

  center_1 <- list()
  scatter_1 <- list()
  subset_1 <- list()

  center_2 <- list()
  scatter_2 <- list()
  subset_2 <- list()

  center_3 <- list()
  scatter_3 <- list()
  subset_3 <- list()

  final_sets <- list()
  final_scatter <- list()

  list_sets <- list()
  S_final <- list()
  sigma <- list()
  center <- list()

  # Iterate over the number of initial subsets
  for (i in 1:m) {
    # Randomly select p + 1 observations to form the initial subset
    subset_0[[i]] <- x[sample(1:n, p + 1), ]
    center_0[[i]] <- apply(subset_0[[i]], 2, mean)
    scatter_0[[i]] <- (p / (p + 1)) * var(subset_0[[i]])

    # Ensure the covariance matrix is non-singular
    while (det(scatter_0[[i]]) == 0) {
      for (k in 1:(n - p - 1)) {
        subset_new <- x[sample(n, (p + 1) + k), ]
      }
      subset_0[[i]] <- subset_new
      center_0[[i]] <- apply(subset_new, 2, mean)
      scatter_0[[i]] <- ((p + k) / (p + 1 + k)) * var(subset_new)
    }

    # Compute Mahalanobis distances for the current subset
    mahal_dist0 <- mahalanobis(x, center_0[[i]], scatter_0[[i]], tol = -1)
    d_sqr0 <- sort(mahal_dist0, decreasing = FALSE, index.return = TRUE)
    idx_h1 <- d_sqr0$ix[1:h]
    subset_1[[i]] <- x[idx_h1, ]
    center_1[[i]] <- apply(subset_1[[i]], 2, mean)
    scatter_1[[i]] <- ((h - 1) / h) * var(subset_1[[i]])

    mahal_dist1 <- mahalanobis(x, center_1[[i]], scatter_1[[i]], tol = -1)
    d_sqr1 <- sort(mahal_dist1, decreasing = FALSE, index.return = TRUE)
    idx_h2 <- d_sqr1$ix[1:h]
    subset_2[[i]] <- x[idx_h2, ]
    center_2[[i]] <- apply(subset_2[[i]], 2, mean)
    scatter_2[[i]] <- ((h - 1) / h) * var(subset_2[[i]])

    mahal_dist2 <- mahalanobis(x, center_2[[i]], scatter_2[[i]], tol = -1)
    d_sqr2 <- sort(mahal_dist2, decreasing = FALSE, index.return = TRUE)
    idx_h3 <- d_sqr2$ix[1:h]
    subset_3[[i]] <- x[idx_h3, ]
    center_3[[i]] <- apply(subset_3[[i]], 2, mean)
    scatter_3[[i]] <- ((h - 1) / h) * var(subset_3[[i]])
    determinant_vec[i] <- det(scatter_3[[i]])
  }

  # Select the l subsets with the smallest determinants
  D_sorted <- sort(determinant_vec, decreasing = FALSE, index.return = TRUE)
  idx_top_l <- D_sorted$ix[1:l]

  # Store the selected subsets and their covariance matrices
  for (r in 1:l) {
    final_sets[[r]] <- subset_3[[idx_top_l[r]]]
    final_scatter[[r]] <- scatter_3[[idx_top_l[r]]]
  }

  object_mcd <- list(
    call = match.call(),
    H0sets = subset_0,
    H1sets = subset_1,
    H2sets = subset_2,
    H3sets = subset_3,
    final_l_set = final_sets,
    finalS = final_scatter,
    index_l = idx_top_l
  )

  # Perform C-steps for convergence on the selected subsets
  best_centers <- list()
  best_scatters <- list()
  best_subsets <- list()

  for (k in 1:l) {
    subset_current <- data.frame(as.matrix(object_mcd$final_l_set[[k]], nrow = h, ncol = p))
    det_new <- 1
    det_old <- 2

    while (det_new != det_old && det_new != 0) {
      center_old <- apply(subset_current, 2, mean)
      scatter_old <- ((h - 1) / h) * var(subset_current)
      det_old <- det(scatter_old)

      mahal_dist_current <- mahalanobis(x, center_old, scatter_old, tol = -1)
      d_sqr_current <- sort(mahal_dist_current, decreasing = FALSE, index.return = TRUE)
      idx_new <- d_sqr_current$ix[1:h]

      # Update subset with new indices
      subset_old <- subset_current
      subset_current <- x[idx_new, ]
      center_new <- apply(subset_current, 2, mean)
      scatter_new <- ((h - 1) / h) * var(subset_current)
      det_new <- det(scatter_new)
    }

    best_centers[[k]] <- center_new
    best_scatters[[k]] <- scatter_new
    best_subsets[[k]] <- subset_current
  }

  # Select the subset with the smallest determinant from the final subsets
  final_dets <- sapply(best_scatters, det)
  final_sorted <- sort(final_dets, decreasing = FALSE, index.return = TRUE)
  best_idx <- final_sorted$ix[1]

  center_mcd <- best_centers[[best_idx]]
  scatter_mcd <- c_alpha * cnp * best_scatters[[best_idx]]
  best_subset <- best_subsets[[best_idx]]
  idx_best_subset <- as.integer(rownames(best_subset))

  # Raw results
  result_raw <- list(
    raw.center = center_mcd,
    raw.cov = scatter_mcd,
    best = idx_best_subset
  )

  # Reweighting step
  weights_vector <- rep(0, n)
  dist_raw <- mahalanobis(x, result_raw$raw.center, result_raw$raw.cov, tol = -1)
  cutoff <- qchisq(1 - delta, df = p)
  weights_vector[dist_raw <= cutoff] <- 1
  weights <- t(weights_vector)

  # Reweighted MCD location estimator
  x_matrix <- as.matrix(x, nrow = n, ncol = p)
  center_rwgt <- as.numeric((weights %*% x_matrix) / sum(weights))

  # Reweighted MCD covariance estimator
  S_accum <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
    diff <- matrix(x_matrix[i, ] - center_rwgt, ncol = 1)
    S_accum <- S_accum + weights[i] * (diff %*% t(diff))
  }
  scatter_rwgt <- (S_accum / sum(weights)) * c_delta * cnp

  # Final output combining raw and reweighted estimates
  result <- list(
    raw.center = result_raw$raw.center,
    raw.cov = result_raw$raw.cov,
    best = result_raw$best,
    weights = weights,
    center = center_rwgt,
    cov = scatter_rwgt
  )

  return(result)
}


#' Quality Control by Identifying Jumps in a Matrix
#'
#' This function performs quality control on a matrix by identifying jumps or outliers.
#' It uses a combination of bubble filling, peak detection, and clustering to identify
#' regions of the matrix that are significantly different from the rest.
#'
#' @param matrix A numeric matrix containing the data to be analyzed.
#' @param var_tol A numeric value specifying the tolerance for variation in peak detection.
#' Defaults to 0.5.
#' @param max_iterations The maximum number of iterations to perform in the clustering process.
#' Defaults to 100.
#' @param kernel A numeric vector specifying the kernel to use in the convolution step.
#' Defaults to c(1,1,1).
#' @return A list containing the time gaps, time clusters, index of good data points,
#' and the iterations matrix.
#' @export Qualitycontrol.ICOD
#' @importFrom stats kmeans
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)

Qualitycontrol.ICOD <- function(matrix, var_tol=0.5, max_iterations=100, kernel = c(1,1,1)){
  resolution <- (max(as.numeric(colnames(matrix))) - min(as.numeric(colnames(matrix))))/ncol(matrix)
  group <- rep(1, nrow(matrix))
  group[is.na(rowSums(matrix))] <- 0
  time_gap <- c()
  time_cluster <- c()
  interations <- group
  # cl <- makeCluster(getOption("cl.cores", detectCores() ))
  index_matrix <- group == 1
  for (inter in 1:max_iterations){
    print(inter)
    new_group <- as.numeric(group)
    index_good <- group == 1
    start_time <- Sys.time()
    peaks <- bubblefill(colMeans(matrix[index_good,]),min_bubble_widths = floor(100/resolution))$peaks
    # peaks <- 1:5
    # peaks <- phenofit::findpeaks(colMeans(matrix[index_good,]), nups = 3, minpeakdistance = 0)$X$pos
    q <- apply(matrix[index_good,peaks],2,quantile)
    IQR <- q[4,]-q[2,]
    gap <- 1.5
    top_range <- q[4,] + IQR * gap
    bottom_range <- q[2,] - IQR * gap
    group_matrix <- t(t(matrix[index_matrix,peaks]) > bottom_range & t(matrix[index_matrix,peaks]) < top_range)
    end_time <- Sys.time()
    time_gap <- c(time_gap, end_time - start_time)
    if(min(group_matrix)==1){
      return(list(time_gap=time_gap,time_cluster=time_cluster,index_good=index_good,interations=as.matrix(interations)))
    } else if(max(group_matrix)==0){
      stop('The gap level is too small, please change a larger one!')
    }
    start_time <- Sys.time()
    group_matrix <- t(apply(!group_matrix, 1, convolve_custom, kernel = kernel))
    kmeans_group <- kmeans(group_matrix, 2, 1)$cluster
    asso_group <- names(which.max(table(kmeans_group[index_matrix[index_good]])))
    new_group[index_matrix][kmeans_group!=asso_group] <- 0
    end_time <- Sys.time()
    time_cluster <- c(time_cluster, end_time - start_time)
    interations <- cbind(interations, new_group)
    if(all.equal(group, new_group) == TRUE ){
      print(max(rowSums(group_matrix[group[index_matrix]!=0,] > 1.1)))
      print(length(peaks))
      if(max(rowSums(group_matrix[group[index_matrix]!=0,] > 1.1)) < ceiling(length(peaks)*var_tol))
        break
      else
        index_matrix <- group==1
    } else group <- new_group
  }
  return(list(time_gap=time_gap,time_cluster=time_cluster,index_good=index_good,interations=as.matrix(interations)))
}


#' Detect outliers in a dataset using Hotelling's T2 method
#'
#' This function identifies outliers in a dataset by applying Hotelling's T2 statistical method.
#' It uses parallel processing to speed up the computation.
#'
#' @param x A matrix of spectral data where each row represents a single spectrum.
#' @return A data frame containing two columns: 'out' indicating whether each spectrum is an outlier,
#' and 'dis' containing the corresponding Hotelling's T2 distance values.
#' @export Qualitycontrol.T2
#' @importFrom disprofas get_hotellings
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' #qc_t2 <- Qualitycontrol.T2(data_normalized@datasets$normalized.data)
Qualitycontrol.T2 <- function(x){
  out_na <- which(is.na(rowSums(x)))
  pred_outliers <- rep(TRUE, nrow(x))
  if(length(out_na)!=0)pred_outliers[out_na,] <- FALSE
  n_out <- 0
  temp_outliers <- pred_outliers
  i <- 1
  while(TRUE){
    print(i)
    mean_spec <- as.matrix(apply(x[temp_outliers,],2, median))
    temp_outliers <- pred_outliers
    hotel_p <- apply(as.matrix(x),1,function(x,spec=mean_spec){
      return(disprofas::get_hotellings(as.matrix(x), spec, 0.05)$Parameters['p.F'])
    })
    temp_outliers[hotel_p < 0.05] <- FALSE
    if(n_out==length(temp_outliers[!temp_outliers]))break
    if(i > 30)break
    i <- i +1
    n_out <- length(pred_outliers[!temp_outliers])
  }
  return(data.frame(quality=temp_outliers, p.F=hotel_p))
}

#' Detect outliers in a dataset based on distance from the mean spectrum
#'
#' This function identifies outliers in a dataset by calculating the distance of each spectrum
#' from the mean spectrum of non-outlier spectra. The distance is calculated as the square
#' root of the sum of the squared differences between each spectrum and the mean spectrum.
#' Spectra with a distance greater than a specified threshold are marked as outliers.
#'
#' @param x A matrix of spectral data where each row represents a single spectrum.
#' @param min.dis The minimum distance threshold for marking a spectrum as an outlier.
#'                 Defaults to 1.
#' @return A data frame containing two columns: 'out' indicating whether each spectrum is an
#'         outlier, and 'dis' containing the corresponding distance values.
#' @export Qualitycontrol.Dis
#' @importFrom stats quantile
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' qc_dis <- Qualitycontrol.Dis(data_normalized@datasets$normalized.data)
Qualitycontrol.Dis <- function(x, min.dis=1){
  pred_outliers <- rep(TRUE, nrow(x))
  out_na <- which(is.na(rowSums(x)))
  if(length(out_na)!=0)pred_outliers[out_na,] <- FALSE
  n_out <- length(pred_outliers[!pred_outliers])
  mean_spec <- colMeans(x[pred_outliers,])
  dis_matrix <- apply(x, 1, function(row) sqrt(sum((row - mean_spec)^2)))
  temp.dis <- quantile(dis_matrix, probs = 0.9)
  while(TRUE){
    pred_outliers[dis_matrix>max(temp.dis, min.dis)] <- FALSE
    mean_spec <- colMeans(x[pred_outliers,])
    dis_matrix <- apply(x, 1, function(row) sqrt(sum((row - mean_spec)^2)))
    if(n_out==length(pred_outliers[!pred_outliers]))  {temp.dis <- temp.dis*0.95;pred_outliers[dis_matrix>max(temp.dis, min.dis)] <- FALSE}
    else{temp.dis <- quantile(dis_matrix, probs = 0.9)}
    n_out <- length(pred_outliers[!pred_outliers])
    if(all(dis_matrix[pred_outliers]<min.dis))
      break
  }
  return(data.frame(out=pred_outliers, dis=dis_matrix))
}

#' Detect outliers using the Minimum Covariance Determinant (MCD) method
#'
#' This function extends the outliers_mcdEst function to detect outliers in a dataset
#' using the Minimum Covariance Determinant method. It modifies the output to include
#' additional information such as the distance and center of the outliers.
#'
#' @param x The input dataset.
#' @param index_good A vector indicating the indices of good observations.
#' @param h The fraction of the data to use for computing the MCD.
#' @param alpha The significance level for the MCD.
#' @param na.rm A logical value indicating whether to remove NA values.
#' @return An object of class "Qualitycontrol.Mcd" containing the following components:
#'   - distance: The maximum distance threshold for outliers.
#'   - center: The center of the outliers.
#'   - outliers_pos: The positions of the outliers.
#'   - nb: A named number giving the total number of outliers detected.
#' @export Qualitycontrol.Mcd
#' @importFrom stats qchisq
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' #qc_mcd <- Qualitycontrol.Mcd(data_normalized@datasets$normalized.data)

Qualitycontrol.Mcd <- function(x,index_good,h = .5,alpha = .01, na.rm = TRUE){
  x <- prcomp_irlba(x, n = 20, center = TRUE, scale. = TRUE)$x[,1:5]
  out <- outliers_mcdEst(x,index_good,h,alpha,na.rm)
  out$distance <- out$MaxDist
  out$center <- out$center
  out$call <- match.call()
  out$nb <- c(total = length(out$outliers_pos))

  class(out) <- "Qualitycontrol.Mcd"
  return(out)
}

#' Calculate the Signal-to-Noise Ratio (SNR) for a dataset
#'
#' This function calculates the SNR for a dataset by comparing the signal and noise levels.
#'
#' @param data A matrix of spectral data where each row represents a single spectrum.
#' @param level The level of strictness for the SNR calculation. Defaults to "medium".
#' @return A vector indicating whether each spectrum passes the SNR test.
#' @export Qualitycontrol.Snr
#' @importFrom stats sd
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' qc_snr <- Qualitycontrol.Snr(data_normalized@datasets$normalized.data)

Qualitycontrol.Snr <- function(data, level = "medium")
{
  data <- as.matrix(data)
  waves <- as.numeric(colnames(data))
  band_noise <- c(which((waves > 1800) & (waves < 1900)))
  band_signal <- c(which(waves > 600 & waves < 3050))
  if (level == 'hard')
  {
    QC1 <- 0.05
    QC2 <- 0.03
    QC3 <- 10
  }
  if (level == 'medium')
  {
    QC1 <- 0.03
    QC2 <- 0.03
    QC3 <- 5
  }
  if (level == 'easy')
  {
    QC1 <- 0.3
    QC2 <- 0.3
    QC3 <- 3
  }
  IBG <- apply(data[,band_noise],1,mean)
  IRAM <- apply(data[,band_signal],1,max)-IBG
  sdBG <- apply(data[,band_signal],1,sd)
  SNR <- IRAM/sdBG


  index_good <- apply(abs(data[,band_noise]),1,function(x) {return((sd(x)<QC2)&(mean(x)<QC1))}) &
    SNR > QC3 &     !duplicated(data)
  index_good[is.na(index_good)] <- 'FALSE'
  return(index_good)
}

#' Detect outliers in a dataset using various quality control methods
#'
#' This function applies multiple quality control methods to detect outliers in a dataset.
#' It uses Hotelling's T2, SNR, MCD, and distance methods to identify and mark outliers.
#'
#' @param matrix A matrix of spectral data where each row represents a single spectrum.
#' @param var_tol Tolerance for variation in the data. Defaults to 0.5.
#' @return A list containing the outliers and the run time for each method.
#' @export Qualitycontrol.All
#' @importFrom disprofas get_hotellings
#' @importFrom stats quantile
#' @importFrom magrittr %<>%
#' @importFrom stats mahalanobis
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' #qc_all <- Qualitycontrol.All(data_normalized@datasets$normalized.data)
Qualitycontrol.All <- function(matrix, var_tol=0.5){
  data <- new('Ramanome', datasets = list(data=matrix), wavenumber=as.numeric(colnames(matrix)))
  data %<>% Preprocesssing.Background.Spike(.,"CPU") %>% Preprocesssing.Smooth.Sg(.,m = 0, p = 5, w = 7, delta.wav = 2) %>% Preprocesssing.Baseline.Bubble %>% Preprocesssing.Normalize(.,'ch')
  matrix <- data@datasets$normalized.data
  outliers_all_ <- outliers_all <- data.frame(matrix('Bad',ncol=5, nrow=nrow(matrix)))
  colnames(outliers_all_) <- colnames(outliers_all) <- c('T2','SNR','mcd','dis','jump') #
  cal_time <- list()
  outlier_na <- which(is.na(rowSums(matrix)))
  if(length(outlier_na)!=0) {matrix=matrix[-outlier_na,];outliers_all=outliers_all_[-outlier_na,]}
  print('Start Hotelling T2!')
  start_time <- Sys.time()
  outliers_all$T2 <- Qualitycontrol.T2(matrix)$out %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$T2 <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start SNR!')
  start_time <- Sys.time()
  outliers_all$SNR <- Qualitycontrol.Snr(matrix, level = 'medium') %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$SNR <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start PC-MCD!')
  start_time <- Sys.time()
  #data.red <- prcomp_irlba(matrix, n = 20, center = TRUE, scale. = TRUE)$x[,1:5]
  outliers_all$mcd <- rep('Good', nrow(outliers_all))
  outliers_all$mcd[Qualitycontrol.Mcd(matrix,index_good=outliers_all$mcd=='Good', alpha = 0.001)$outliers_pos] <- 'Bad'
  end_time <- Sys.time()
  cal_time$mcd <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start Distance!')
  start_time <- Sys.time()
  outliers_all$dis <- Qualitycontrol.Dis(matrix, min.dis=ncol(matrix)*0.01)$out %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$dis <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start Jump!')
  start_time <- Sys.time()
  outliers_all$jump <- Qualitycontrol.ICOD(matrix, var_tol=var_tol, max_iterations=100)$index_good %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$jump <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if(length(outlier_na)!=0) outliers_all_[-outlier_na,]=outliers_all else outliers_all_=outliers_all
  return(list(outliers=outliers_all_, run_time=cal_time))
}
