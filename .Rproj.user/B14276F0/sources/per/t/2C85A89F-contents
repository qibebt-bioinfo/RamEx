#' find spikes
#'
#' @param mat = matrix
#' @return a spike_list(sample_i_ind and wavenumber_ind)
#' @examples
#' pre_spike_find(hyperspec, width = 5)
#'
convolve_matrix <- function(mat, kernel) {


  rows <- sapply(1:nrow(mat), function(i) {
    row <- mat[i, ]
    stats::convolve(row, kernel, type = "open")
  })

  cols <- sapply(1:ncol(rows), function(i) {
    col <- rows[, i]
    convolve(col, kernel, type = "open")
  })

  result <- t(cols[(length(kernel)+1):(length(kernel)+ncol(mat)),])
  return(result)
}


#' @import parallel
convolve_matrix <- function(mat, kernel) {
  mat <- as.matrix(mat)

  num_cores <- detectCores() - 1


  cl <- makeCluster(num_cores)


  clusterExport(cl, c("mat", "kernel"))


  rows <- parSapply(cl, 1:nrow(mat), function(i) {
    row <- mat[i, ]
    stats::convolve(row, kernel, type = "open")
  })


  stopCluster(cl)


  rows <- t(rows)


  cl <- makeCluster(num_cores)


  clusterExport(cl, c("rows", "kernel"))


  cols <- parSapply(cl, 1:ncol(rows), function(i) {
    col <- rows[, i]
    stats::convolve(col, kernel, type = "open")
  })


  stopCluster(cl)


  result <- t(cols[(length(kernel) + 1):(length(kernel) + ncol(mat)), ])

  return(result)
}


#' @import Rcpp
cppFunction('
NumericMatrix convolve_matrix_rcpp(NumericMatrix mat, NumericVector kernel) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int ksize = kernel.size();

  NumericMatrix result(nrow, ncol);

  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      double sum = 0.0;
      for (int k = 0; k < ksize; ++k) {
        int idx = j - k + ksize / 2;
        if (idx >= 0 && idx < ncol) {
          sum += mat(i, idx) * kernel[k];
        }
      }
      result(i, j) = sum;
    }
  }

  return result;
}
')


pre_spike_matrix <- function(all_data) {

  all_spc <- as.matrix(all_data[,-1])
  all_spc <- matrix(all_spc, nrow = nrow(all_spc), ncol = ncol(all_spc), byrow = FALSE)
  # slope_inds_1 <- which(all_spc > 1,arr.ind = TRUE)
  kernel_ori <- matrix(c(-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,33, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1),nrow = 3)

  all_spc_1 <- convolve_matrix_rcpp(all_spc, kernel_ori)
  print('=======aaa=======')
  slope_inds_2 <- which(all_spc_1 > 10*apply(all_spc,1,max),arr.ind = TRUE)

  # slope_inds <- unique(rbind(slope_inds_1, slope_inds_2))
  slope_inds <- slope_inds_2

  spc_new <- all_data
  wavenumber <- as.numeric(as.character(colnames(spc_new[,-1])))

  print(length(unique(slope_inds[,1])))

  for (ind in unique(slope_inds[,1])) {
    spike_pos <- slope_inds[which(slope_inds[,1] == ind),2]
    if(ind <= 3 ) { inds_new <- c((ind+1):(ind+3))}
    if(ind > 3 & ind <= (nrow(all_data)-3)) { inds_new <- c((ind-3):(ind-1),(ind+1):(ind+3))}
    if(ind > (nrow(all_data)-3) ) { inds_new <- c((ind-3):(ind-1))}

    for(i in spike_pos){
      if (i == 1) {
        spc_new[ind, 2:4] <- colMeans(spc_new[inds_new,2:4])
      } else if (i >= (length(wavenumber))) {
        spc_new[ind, (i-2):i] <- colMeans(spc_new[inds_new,(i-2):i])
      } else {
        spc_new[ind, i:(i+2)] <- colMeans(spc_new[inds_new,i:(i+2)])
      }
    }
  }

  return(spc_new)
}



matrix2vector <- function(image, kernel_height, kernel_width) {
  image_height <- nrow(image)
  image_width <- ncol(image)
  output_height <- image_height - kernel_height + 1
  output_width <- image_width - kernel_width + 1

  cols <- matrix(0, nrow = kernel_height * kernel_width, ncol = output_height * output_width)
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

gpupre_spike_matrix <- function(all_data) {
  all_spc <- as.matrix(all_data[,-1])

  if (nrow(all_spc) < 3 || ncol(all_spc) < 11) {
    stop("all_spc must be at least 3 rows and 11 columns")
  }

  kernel_ori <- matrix(c(-1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,33, -1,-1,-1, -1,-1,
                         -1,-1,-1, -1,-1,-1, -1,-1,-1, -1,-1),nrow = 3, byrow = TRUE)
  vector_result <- matrix2vector(all_spc, 3, 11)

  kernel_matrix <- as.vector(t(kernel_ori))

  conv_result <-  kernel_matrix %*% vector_result

  all_spc_1 <- matrix(conv_result, nrow = nrow(all_spc) - 2, ncol = ncol(all_spc) - 10, byrow = TRUE)

  slope_inds_2 <- which(all_spc_1 > 10 * apply(all_spc, 1, max), arr.ind = TRUE)
  slope_inds_2 <- which(all_spc_1 > 10 * apply(all_spc, 1, max), arr.ind = TRUE)

  slope_inds <- slope_inds_2

  spc_new <- all_data
  wavenumber <- as.numeric(as.character(colnames(spc_new[,-1])))

  print(length(unique(slope_inds[,1])))

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
#'
#' @return A rammanome object
#' @export pre.spike
pre.spike <- function(object){
  data <- get.nearest.dataset(object)
  pre.data <- gpupre_spike_matrix(data)
  object@datasets$spike.data <- pre.data
  return(object)
}


