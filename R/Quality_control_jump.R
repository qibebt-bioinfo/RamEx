#' Polynomial fitting for spectral baseline correction
#'
#' This function performs polynomial fitting to correct the baseline of spectral data.
#' It iteratively fits a polynomial to the data, ensuring that the fitted values
#' do not exceed the original data values.
#'
#' @param spectra A matrix of spectral data where each row represents a spectrum.
#' @param t A numeric vector of x-values for the polynomial fitting. If missing or FALSE,
#'   a sequence from 1 to the number of columns in `spectra` is used.
#' @param degree The degree of the polynomial to fit.
#' @param tol The tolerance for convergence of the iterative fitting process.
#' @param rep The maximum number of iterations to perform.
#' @return A list containing the baseline and the corrected spectra.
#' @export polyfit
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
#' @export grow_bubble
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
#' @return A vector of the larger values.
#' @export keep_largest
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
#' @export bubbleloop
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
#' @export bubblefill
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

#' Quality Control by Identifying Jumps in a Matrix
#'
#' This function performs quality control on a matrix by identifying jumps or outliers.
#' It uses a combination of bubble filling, peak detection, and clustering to identify
#' regions of the matrix that are significantly different from the rest.
#'
#' @param matrix A numeric matrix containing the data to be analyzed.
#' @param var_tol A numeric value specifying the tolerance for variation in peak detection.
#'   Defaults to 0.5.
#' @param max_iterations The maximum number of iterations to perform in the clustering process.
#'   Defaults to 100.
#' @param kernel A numeric vector specifying the kernel to use in the convolution step.
#'   Defaults to c(1,1,1).
#' @return A list containing the time gaps, time clusters, index of good data points,
#'   and the iterations matrix.
#' @export qc_jump
#' @importFrom stats kmeans
qc_jump <- function(matrix, var_tol=0.5, max_iterations=100, kernel = c(1,1,1)){
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


#' Custom Convolution Function
#'
#' This function performs a convolution operation on a numeric vector using a specified kernel.
#' The convolution is performed using the filter type, which applies the kernel to the vector.
#'
#' @param vec A numeric vector to be convolved.
#' @param kernel A numeric vector representing the convolution kernel.
#' @return A numeric vector resulting from the convolution operation.
#' @export convolve_custom
#' @importFrom stats convolve
convolve_custom <- function(vec, kernel) {
  result <- convolve(as.numeric(vec), rev(kernel), type = "filter")
  return(result)
}


