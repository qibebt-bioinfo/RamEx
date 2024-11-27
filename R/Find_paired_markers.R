#' Compute Correlations with a Target Variable
#'
#' This function calculates the correlations between each column of a dataset X and a target
#' variable y. It also computes the correlations for combinations of two and three
#' variables, after applying certain transformations.
#'
#' @param X A data frame containing the dataset.
#' @param y A vector containing the target variable.
#' @param min.cor The minimum correlation threshold. Defaults to 0.8.
#' @return A list containing the high correlation elements, individual correlations,
#' correlations for two-variable combinations, and correlations for three-variable combinations.
#' @export compute_correlations
#' @importFrom dplyr filter
#' @importFrom stats cor
#' @importFrom utils combn
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom doParallel registerDoParallel
compute_correlations <- function(X, y,min.cor=0.8) {
  X <- as.data.frame(X)


  correlations <- sapply(X, function(col) cor(col, y, use = "complete.obs"))

  high_cor_elements <- names(which(correlations > min.cor))


  combinations_two <- expand.grid(col1 = colnames(X), col2 = colnames(X)) %>%
    filter(col1 != col2)

  combinations_three <- t(combn(colnames(X), 3))

  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  clusterExport(cl, varlist = c("X", "y"), envir = environment())


  print('=========== 2 ============')
  combination_correlations_two <- do.call(rbind, parLapply(cl, 1:nrow(combinations_two), function(i) {
    row <- combinations_two[i, ]
    col1 <- row$col1
    col2 <- row$col2
    cor_div <- try(cor(X[[col1]] / X[[col2]], y, use = "complete.obs"), silent = TRUE)
    if (!inherits(cor_div, "try-error") && abs(cor_div) > min.cor) {
      data.frame(col1, col2, cor_div)
    } else {
      NULL
    }
  }))


  print('=========== 3 ============')
  combination_correlations_three <- do.call(rbind, parLapply(cl, 1:nrow(combinations_three), function(i) {
    row <- combinations_three[i, ]
    col1 <- row[1]
    col2 <- row[2]
    col3 <- row[3]
    comb_value <- (X[[col1]] - X[[col2]]) / X[[col3]]
    cor_comb <- try(cor(comb_value, y, use = "complete.obs"), silent = TRUE)
    if (!inherits(cor_comb, "try-error") && abs(cor_comb) > min.cor) {
      data.frame(col1, col2, col3, cor_comb)
    } else {
      NULL
    }
  }))

  # 停止并行计算
  stopCluster(cl)

  list(
    high_cor_elements = high_cor_elements,
    correlations = correlations,
    combination_correlations_two = combination_correlations_two,
    combination_correlations_three = combination_correlations_three
  )
}
#' Find Markers Using ROC Analysis
#'
#' This function performs Receiver Operating Characteristic (ROC) analysis to identify
#' markers in a dataset that are associated with different groups. It calculates the
#' Area Under the Curve (AUC) for single markers and paired markers.
#'
#' @param matrix A matrix of spectral data where each row represents a sample and each
#' column represents a wavelength.
#' @param group A vector of group identifiers corresponding to each row in the matrix.
#' @param threshold The minimum AUC threshold for considering a marker significant.
#' Defaults to 0.75.
#' @return A list containing two data frames: 'markers_single' for single markers and
#' 'markers_paired' for paired markers, each with their corresponding AUC scores.
#' @export find_markers_roc
#' @importFrom utils combn
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom MLmetrics AUC
find_markers_roc <- function(matrix,group, threshold=0.75){
  wave <- as.numeric(colnames(matrix))
  u_group <- unique(group)

  combinations_two <- as.data.frame(t(combn(colnames(matrix), 2)))
  colnames(combinations_two) <- c('col1','col2')

  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  clusterEvalQ(cl, library(MLmetrics))
  clusterExport(cl, varlist = c("group","combinations_two","wave","matrix", 'threshold'), envir = environment())


  # Single
  cat('Finding single markers ... \n')
  raman_markers <- lapply(u_group,
                          function(x)parApply(cl,matrix,2,function(i)MLmetrics::AUC(i,group==x)))
  names(raman_markers) <- u_group
  markers_single <- lapply(names(raman_markers), function(name) {
    data.frame(group = rep(name, length(raman_markers[[name]])), wave=names(raman_markers[[name]]),auc = unlist(raman_markers[[name]]))
  })

  markers_single <- do.call(rbind, markers_single)

  # Paired
  cat('Finding paired markers ... \n')

  markers_paired <- lapply(u_group, function(x) parLapply(cl, 1:nrow(combinations_two), function(i) {
    row <- combinations_two[i, ]
    col1 <- row$col1
    col2 <- row$col2
    auc_score <- MLmetrics::AUC(matrix[,wave==col1] / matrix[,wave==col2], group==x)
    print(auc_score)
    if ( auc_score > threshold) {
      return(data.frame(col1, col2, auc_score))
    }
  }))

  stopCluster(cl)

  names(markers_paired) <- u_group
  markers_paired <- do.call(rbind, lapply(names(markers_paired), function(name) {
    df <- markers_paired[[name]]
    df <- do.call(rbind,Filter(Negate(is.null), df))
    df$source <- name
    return(df[,c(4,1:3)])
  }))

  return(list(markers_single = markers_single, markers_paired = markers_paired))
}
