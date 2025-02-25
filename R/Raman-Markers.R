#' Calculate mislabel scores for classification results
#'
#' This function calculates mislabel scores based on class probabilities for classification results.
#'
#' @param y A vector containing the true labels.
#' @param y.prob A matrix containing class probabilities.
#'
#' @return A matrix containing mislabel scores for each class.
#' @importFrom stats model.matrix


get.mislabel.scores <- function(y, y.prob) {
  result <- matrix(0, nrow = length(y), ncol = 3)
  mm <- stats::model.matrix(~0 + y)
  y.prob.other.max <- apply(y.prob * (1 - mm), 1, max)
  y.prob.alleged <- apply(y.prob * mm, 1, max)
  result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
  rownames(result) <- rownames(y.prob)
  colnames(result) <- c('P(alleged label)', 'P(second best)', 'P(alleged label)-P(second best)')
  return(result)
}

#' Generate balanced folds for cross-validation
#'
#' This function generates balanced folds for cross-validation, ensuring an equal number of samples from each class in each fold.
#'
#' @param y A vector containing the class labels.
#' @param nfolds The number of folds to generate.
#'
#' @return A vector containing the fold assignments for each sample.


balanced.folds <- function(y, nfolds = 10) {
  folds <- rep(0, length(y))
  classes <- levels(y)
  Nk <- table(y)
  if (nfolds == -1 || nfolds == length(y)) {
    invisible(seq_along(y))
  } else {
    nfolds <- min(nfolds, max(Nk))
    for (k in seq_along(classes)) {
      ixs <- which(y == classes[k])
      folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
      folds_k <- folds_k[seq_along(ixs)]
      folds_k <- sample(folds_k)
      folds[ixs] <- folds_k
    }
    invisible(folds)
  }
}

#' Perform cross-validation for random forest classification
#'
#' This function performs cross-validation for random forest classification on Raman data.
#'
#' @param x A data frame or matrix containing the features.
#' @param y A vector containing the class labels.
#' @param nfolds The number of folds for cross-validation. Default is 10.
#' @param verbose A logical value indicating whether to display progress messages. Default is FALSE.
#' @param ... Additional arguments to be passed to the randomForest function.
#'
#' @return A list containing the cross-validation results including predicted labels, class probabilities, feature importances, and error rates.
#' @importFrom randomForest randomForest

rf.cross.validation <- function(x, y, nfolds = 10, verbose = FALSE, ...) {
  if (nfolds == -1) nfolds <- length(y)
  folds <- balanced.folds(y, nfolds = nfolds)
  result <- list()
  result$y <- as.factor(y)
  result$predicted <- result$y
  result$probabilities <- matrix(0, nrow = length(result$y), ncol = length(levels(result$y)))
  rownames(result$probabilities) <- rownames(x)
  colnames(result$probabilities) <- levels(result$y)
  result$importances <- matrix(0, nrow = ncol(x), ncol = nfolds)
  result$errs <- numeric(length(unique(folds)))
  
  # K-fold cross-validation
  for (fold in sort(unique(folds))) {
    if (verbose) cat(sprintf('Fold %d...\n', fold))
    foldix <- which(folds == fold)
    model <- randomForest::randomForest(
      x[-foldix,],
      factor(result$y[-foldix]),
      importance = TRUE,
      do.trace = verbose, ...
    )
    newx <- x[foldix,]
    if (length(foldix) == 1) newx <- matrix(newx, nrow = 1)
    result$predicted[foldix] <- stats::predict(model, newx)
    probs <- stats::predict(model, newx, type = 'prob')
    result$probabilities[foldix, colnames(probs)] <- probs
    result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
    result$importances[, fold] <- model$importance[, 'MeanDecreaseAccuracy']
  }
  
  result$nfolds <- nfolds
  result$params <- list(...)
  result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y == level])))
  return(result)
}


#' Perform out-of-bag analysis for random forest classification
#'
#' This function performs out-of-bag analysis for random forest classification on Raman data.
#'
#' @param x A matrix or data frame containing the predictors.
#' @param y A vector containing the response variable.
#' @param verbose A logical value indicating whether to print the progress of the model.
#' @param ntree An integer specifying the number of trees to grow.
#' @param ... Additional arguments to be passed to the randomForest function.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{rf.model}{The random forest model object.}
#'   \item{probabilities}{The out-of-bag class probabilities for each observation.}
#'   \item{y}{The true response variable.}
#'   \item{predicted}{The predicted response variable based on the random forest model.}
#'   \item{confusion.matrix}{A confusion matrix showing the count of true and predicted levels.}
#'   \item{params}{A list of parameter values used for the random forest model.}
#'   \item{errs}{A vector of error rates for each observation.}
#'   \item{importances}{A vector of feature importances based on mean decrease accuracy.}
#' }
#' @importFrom randomForest randomForest


rf.out.of.bag <- function(x, y, verbose = FALSE, ntree = 500, ...) {
  rf.model <- randomForest::randomForest(
    x,
    y,
    keep.inbag = TRUE,
    importance = TRUE,
    do.trace = verbose,
    ...
  )
  result <- list()
  result$rf.model <- rf.model
  result$probabilities <- get.oob.probability.from.forest(rf.model, x)
  result$y <- y
  result$predicted <- rf.model$predicted
  result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y == level])))
  result$params <- list(ntree = ntree)
  result$errs <- as.numeric(result$predicted != result$y)
  result$importances <- rf.model$importance[, 'MeanDecreaseAccuracy']
  return(result)
}


#' Get out-of-bag class probabilities from a random forest model
#'
#' This function calculates the out-of-bag class probabilities for each sample in the dataset, based on the aggregated class votes from the out-of-bag trees.
#'
#' @param model The random forest model object.
#' @param x A matrix or data frame containing the predictors.
#'
#' @return A matrix with the predicted class probabilities for each sample.

get.oob.probability.from.forest <- function(model, x) {
  # get aggregated class votes for each sample using only OOB trees
  votes <- get.oob.votes.from.forest(model, x)
  # convert to probs
  probs <- sweep(votes, 1, rowSums(votes), '/')
  rownames(probs) <- rownames(x)
  colnames(probs) <- model$classes
  
  return(probs)
}
#' Get Out-of-Bag class votes from a random forest model
#'
#' This function extracts the aggregated class votes for each sample using only the Out-of-Bag trees from a random forest model.
#'
#' @param model The random forest model.
#' @param x A matrix or data frame containing the predictor variables.
#'
#' @return A matrix of class votes for each sample, with rows representing samples and columns representing classes.
#' @importFrom stats predict


get.oob.votes.from.forest <- function(model, x) {
  # Get aggregated class votes for each sample using only OOB trees
  votes <- matrix(0, nrow = nrow(x), ncol = length(model$classes))
  
  rf.pred <- stats::predict(model, x, type = "vote", predict.all = TRUE)
  for (i in 1:nrow(x)) {
    # Find which trees are not in-bag for this sample
    outofbag <- model$inbag[i, ] == 0
    # Get OOB predictions for this sample
    votes[i, ] <- table(factor(rf.pred$individual[i, ][outofbag], levels = model$classes))
  }
  rownames(votes) <- rownames(x)
  colnames(votes) <- model$classes
  
  return(invisible(votes))
}


#' Save random forest classification results
#'
#' This function saves the various results of random forest classification to the specified output directory.
#'
#' @param result The result object returned by the random forest classification.
#' @param rf.opts The options used for random forest classification.
#' @param feature.ids The ids of the features used for classification.
#' @param outdir The directory where the results will be saved.
#'
#' @return None. The function saves multiple files to the specified output directory:
#'   - A summary file with error metrics and model parameters
#'   - A file containing classification probabilities for each sample
#'   - A file with potential mislabeling information
#'   - A file with feature importance scores
#'   - A confusion matrix file showing classification performance
#'
#' @importFrom utils write.table

save.rf.results <- function(result, rf.opts, feature.ids, outdir) {
  save.rf.results.summary(result, rf.opts, outdir = outdir)
  save.rf.results.probabilities(result, outdir = outdir)
  save.rf.results.mislabeling(result, outdir = outdir)
  save.rf.results.importances(result, feature.ids = feature.ids, outdir = outdir)
  save.rf.results.confusion.matrix(result, outdir = outdir)
}

#' Save summary of random forest classification results
#'
#' This function saves a summary of the random forest classification results to a file.
#'
#' @param result The result object from random forest classification
#' @param rf.opts The options used for random forest classification
#' @param filename The filename for the summary file (default is "summary.xls")
#' @param outdir The output directory to save the summary file
#'
#' @return None

save.rf.results.summary <- function(result, rf.opts, filename = 'summary.xls', outdir) {
  err <- hyperSpec::mean(result$errs)
  err.sd <- stats::sd(result$errs)
  baseline.err <- 1 - max(table(y)) / length(y)
  filepath <- sprintf('%s/%s', outdir, filename)
  sink(filepath)
  cat(sprintf('Model\tRandom Forest\n'))
  cat(sprintf('Error type\t%s\n', result$error.type))
  if (rf.opts$errortype == 'oob' || rf.opts$errortype == 'cvloo') {
    cat(sprintf('Estimated error\t%.5f\n', err))
  } else {
    cat(sprintf('Estimated error (mean +/- s.d.)\t%.5f +/- %.5f\n', err, err.sd))
  }
  cat(sprintf('Baseline error (for random guessing)\t%.5f\n', baseline.err))
  cat(sprintf('Ratio baseline error to observed error\t%.5f\n', baseline.err / err))
  cat(sprintf('Number of trees\t%d\n', result$params$ntree))
  sink(NULL)
}


#' Save probabilities of random forest classification results
#'
#' This function saves the probabilities of random forest classification results to a file.
#'
#' @param result The result object from random forest classification
#' @param filename The filename for the probabilities file (default is "cv_probabilities.xls")
#' @param outdir The output directory to save the probabilities file
#' @return NULL. The function saves a tab-delimited text file containing class probabilities
#'         for each sample, with rows for samples and columns for class probabilities.

save.rf.results.probabilities <- function(result, filename = 'cv_probabilities.xls', outdir) {
  filepath <- sprintf('%s/%s', outdir, filename)
  sink(filepath)
  cat('SampleID\t')
  write.table(result$probabilities, sep = '\t', quote = FALSE)
  sink(NULL)
}

#' Save mislabeling scores from random forest classification
#'
#' @param result The result object from random forest classification
#' @param filename The filename for the probabilities file (default is "mislabeling.xls")
#' @param outdir The output directory to save the probabilities file
#' @return NULL. The function saves a tab-delimited text file containing mislabeling scores
#'         for each sample, with columns for sample ID and corresponding mislabeling score.

save.rf.results.mislabeling <- function(result, filename = 'mislabeling.xls', outdir) {
  filepath <- sprintf('%s/%s', outdir, filename)
  sink(filepath)
  cat('SampleID\t')
  write.table(get.mislabel.scores(result$y, result$probabilities), sep = '\t', quote = FALSE)
  sink(NULL)
}


#' Save feature importance scores of random forest classification results
#'
#' This function saves the feature importance scores of random forest classification results to a file.
#'
#' @param result The result object obtained from random forest classification
#' @param feature.ids A vector of feature identifiers
#' @param filename The name of the output file
#' @param outdir The path to the output directory
#'
#' @return None


save.rf.results.importances <- function(result, feature.ids, filename = 'feature_importance_scores.xls', outdir) {
  filepath <- sprintf('%s/%s', outdir, filename)
  
  if (is.null(dim(result$importances))) {
    imp <- result$importances
    imp.sd <- rep(NA, length(imp))
  } else {
    imp <- rowMeans(result$importances)
    imp.sd <- apply(result$importances, 1, sd)
  }
  
  output.table <- cbind(imp, imp.sd)
  rownames(output.table) <- feature.ids
  output.table <- output.table[sort(imp, decreasing = TRUE, index = TRUE)$ix,]
  colnames(output.table) <- c('Mean_decrease_in_accuracy', 'Standard_deviation')
  
  sink(filepath)
  cat('Feature_id\t')
  write.table(output.table, sep = '\t', quote = FALSE)
  sink(NULL)
}


#' Save confusion matrix of random forest classification results
#'
#' This function saves the confusion matrix of random forest classification results to a file.
#'
#' @param result The result object obtained from random forest classification
#' @param filename The name of the output file
#' @param outdir The path to the output directory
#'
#' @return None


save.rf.results.confusion.matrix <- function(result, filename = 'confusion_matrix.xls', outdir) {
  filepath <- sprintf('%s/%s', outdir, filename)
  
  # add class error column to each row
  x <- result$confusion.matrix
  class.errors <- rowSums(x * (1 - diag(nrow(x)))) / rowSums(x)
  output <- cbind(result$confusion.matrix, class.errors)
  colnames(output)[ncol(output)] <- "Class error"
  
  sink(filepath)
  cat('True\\Predicted\t')
  write.table(output, sep = '\t', quote = FALSE)
  sink(NULL)
}


#' Perform between group statistical tests on Raman data
#'
#' This function performs between group statistical tests on Raman data,
#' including t-tests, Wilcoxon tests, variances tests, Oneway tests, and Kruskal-Wallis tests.
#'
#' @param data A numeric matrix or data frame containing the Raman data.
#' @param group A factor indicating the grouping variable.
#' @param p.adj.method The method used for p-value adjustment. Default is "bonferroni".
#' @param paired Logical indicating whether paired tests should be performed for t-tests and Wilcoxon tests. Default is FALSE.
#'
#' @return A data frame with the results of the statistical tests.
#' @importFrom stats var.test
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' @importFrom stats oneway.test
#' @importFrom stats kruskal.test
#' @importFrom stats bartlett.test
#' @importFrom stats p.adjust
BetweenGroup.test <- function(data, group, p.adj.method = "bonferroni", paired = FALSE) {
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  
  n_group <- nlevels(group)
  if (!is.numeric(n_group) | n_group == 1)
    stop("group must be a numeric and up to two levels\n")
  if (n_group == 2) {
    output1 <- matrix(NA, ncol = 9, nrow = ncol(data))
    rownames(output1) <- colnames(data)
    colnames(output1) <- c(
      paste0("mean_", levels(group)[1]),
      paste0("mean_", levels(group)[2]),
      paste0("sd_", levels(group)[1]),
      paste0("sd_", levels(group)[2]),
      "Var.test",
      "T-test",
      "Wilcoxon.test",
      paste0("T-test_", p.adj.method),
      paste0("Wilcoxon.test_", p.adj.method)
    )
    
    # 计算样本均值和样本标准差
    means <- sapply(levels(group), function(g) colMeans(data[, group == g]))
    sds <- sapply(levels(group), function(g) apply(data[, group == g], 2, sd))
    
    for (i in seq_len(ncol(data)))
    {
      output1[i, 1] <- means[1, i]
      output1[i, 2] <- means[2, i]
      output1[i, 3] <- sds[1, i]
      output1[i, 4] <- sds[2, i]
      output1[i, 5] <- var.test(data[, i] ~ group)$p.value
      
      if (output1[i, 5] < 0.01)
        output1[i, 6] <- t.test(data[, i] ~ group, paired = paired)$p.value
      else
        output1[i, 6] <- t.test(data[, i] ~ group, var.equal = TRUE, paired = paired)$p.value
      
      output1[i, 7] <- wilcox.test(data[, i] ~ group, paired = paired, conf.int = TRUE, exact = FALSE, correct = FALSE)$p.value
    }
    
    output1[, 8] <- p.adjust(output1[, 6], method = p.adj.method, n = ncol(data))
    output1[, 9] <- p.adjust(output1[, 7], method = p.adj.method, n = ncol(data))
    
    return(data.frame(output1))
  }else {
    output2 <- matrix(NA, ncol = n_group + 5, nrow = ncol(data))
    rownames(output2) <- colnames(data)
    colnames.output2 <- array(NA)
    for (j in seq_len(ncol(output2))) {
      if (j <= n_group) {
        colnames.output2[j] <- paste0("mean_", levels(group)[j])
      }else {
        colnames.output2[(n_group + 1):(n_group + 5)] <- c(
          "Var.test", "Oneway-test", "Kruskal.test",
          paste0("Oneway-test_", p.adj.method),
          paste0("Kruskal.test_", p.adj.method)
        )
      }
    }
    colnames(output2) <- colnames.output2
    
    means <- sapply(levels(group), function(g) colMeans(data[, group == g]))
    
    for (i in seq_len(ncol(data)))
    {
      for (j in 1:n_group)
      {
        output2[i, j] <- means[j, i]
      }
      output2[i, (n_group + 1)] <- bartlett.test(data[, i] ~ group)$p.value
      
      if (output2[i, (n_group + 1)] < 0.01)
        output2[i, (n_group + 2)] <- oneway.test(data[, i] ~ group)$p.value
      else
        output2[i, (n_group + 2)] <- oneway.test(data[, i] ~ group, var.equal = TRUE)$p.value
      
      output2[i, (n_group + 3)] <- kruskal.test(data[, i] ~ group)$p.value
    }
    
    output2[, (n_group + 4)] <- p.adjust(output2[, (n_group + 2)], method = p.adj.method, n = ncol(data))
    output2[, (n_group + 5)] <- p.adjust(output2[, (n_group + 3)], method = p.adj.method, n = ncol(data))
    
    return(data.frame(output2))
  }
}

#' Perform random forest classification and feature importance analysis
#'
#' This function performs random forest classification on Raman data and generates feature importances.
#'
#' @param object A Ramanome object.
#' @param outfolder The output folder to save the results.
#' @param draw Logical value indicating whether to draw a heatmap.
#' @param threshold The threshold value for feature importances.
#' @param ntree The number of trees in the random forest model.
#' @param errortype The error type used in random forest model evaluation.
#' @param verbose Logical value indicating whether to display progress messages.
#' @param nfolds The number of folds used in cross-validation.
#' @param rf.cv Logical value indicating whether to perform random forest with cross-validation.
#'
#' @return A data frame containing the selected feature importances.
#' @importFrom pheatmap pheatmap
#' @importFrom hyperSpec aggregate
#' @importFrom dplyr %>%
#' @importFrom magrittr %<>%
#' @importFrom grDevices pdf
#' @export Raman.Markers.Rbcs
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' RBCS.markers <- Raman.Markers.Rbcs(data_cleaned, threshold = 0.003, draw = FALSE)


Raman.Markers.Rbcs <- function(
    object,
    outfolder = "RF_imps",
    draw = TRUE,
    threshold = 0.002,
    ntree = 1000,
    errortype = 'oob',
    verbose = FALSE,
    nfolds = 3,
    rf.cv = FALSE
) {
  x <- get.nearest.dataset(object)
  y <- object@meta.data$group
  
  # Out-of-bag classification in training data and prediction in test data
  y <- as.factor(y)
  oob.result <- rf.out.of.bag(x, y, verbose = verbose, ntree = ntree)
  
  # RF importances of features
  dir.create(outfolder)
  
  imps <- oob.result$importances
  imps <- imps[order(imps, decreasing = TRUE)]
  imps.cutoff <- imps[which(imps > threshold)]
  
  sink(paste0(outfolder, "/RBCS_imps", "_All_sorted.xls"))
  cat("\t")
  write.table(imps, quote = FALSE, sep = "\t")
  sink()
  
  sink(paste0(outfolder, "/RBCS_imps", "_cutoff_sorted.xls"))
  cat("\t")
  write.table(imps.cutoff, quote = FALSE, sep = "\t")
  sink()
  
  data.select <- x[, names(imps.cutoff)]
  means <- aggregate(data.select, by = list(y), mean)
  rownames(means) <- means[, 1]
  means %<>% .[, -1]
  colnames(means) %<>% as.numeric() %>% round(., 0)
  if (draw) {
    pheatmap::pheatmap(
      means[, -1],
      show_colnames = TRUE,
      angle_col = 45,
      scale = "none",
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      filename = 'RBCS.heatmap.png',
      width = length(imps.cutoff) / 4,
      height = length(unique(y)) / 2
    )
  }
  return(imps.cutoff)
}



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
#' @export Raman.Markers.Correlations
#' @importFrom dplyr filter
#' @importFrom stats cor
#' @importFrom utils combn
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom doParallel registerDoParallel
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' options(mc.cores = 2)
#' #cor_markers <- Raman.Markers.Correlations(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],as.numeric(data_cleaned@meta.data$group), min.cor = 0.6)
Raman.Markers.Correlations <- function(X, y,min.cor=0.8) {
  X <- as.data.frame(X)
  waves <- colnames(X)
  
  correlations <- sapply(X, function(col) cor(col, y, use = "complete.obs"))
  
  cat('Finding singular markers ... \n')
  high_cor_elements <- which(abs(correlations) > min.cor)
  
  combinations_two <- expand.grid(col1 = waves, col2 = waves) %>% dplyr::filter(col1 != col2)
  
  numCores <- detectCores() - 1
  cl <- makeCluster(numCores)
  clusterExport(cl, varlist = c("X", "y"), envir = environment())
  
  
  cat('Finding paired markers ... \n')
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
  # 停止并行计算
  stopCluster(cl)
  
  list(
    correlations_singular = if(length(high_cor_elements) == 0) NA else data.frame(wave=waves[high_cor_elements], corrlations = correlations[high_cor_elements]),
    combination_correlations_two = if(length(combination_correlations_two) == 0) NA else combination_correlations_two
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
#' @param threshold The minimum AUC threshold for considering a marker significant.Defaults to 0.75.
#' @param paired Marker is paired or no. Default is FALSE
#' @param batch_size Set the batch size. Default is 1,000
#' @return A list containing two data frames: 'markers_single' for single markers and
#' 'markers_paired' for paired markers, each with their corresponding AUC scores.
#' @export Raman.Markers.Roc
#' @importFrom utils combn
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @importFrom MLmetrics AUC
#' @importFrom Rcpp sourceCpp
#' @import RcppParallel
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' ROC_markers <- Raman.Markers.Roc(data_cleaned@datasets$normalized.data[,sample(1:1000, 50)],data_cleaned@meta.data$group)
Raman.Markers.Roc <- function(matrix, group, threshold = 0.75, paired = FALSE, batch_size = 1000){
  os_type <- .Platform$OS.type
  if (os_type == "unix") {
    file_path <- system.file("libs/RamEx.so", package = "RamEx")
  } else if (os_type == "windows") {
    file_path <- system.file("libs\\x64\\RamEx.dll", package = "RamEx")
  } else {
    stop("Unsupported operating system.")
  }
  dyn.load(file_path)
  wave <- as.numeric(colnames(matrix))
  group <- as.factor(group)
  u_group <- levels(group)
  
  cat('Finding singular markers ... \n')
  raman_markers <- data.frame(calculateAUCParallel(matrix, group))
  names(raman_markers) <- u_group
  raman_markers$wave <- wave
  raman_markers <- reshape2::melt(raman_markers, id.var='wave', variable.name = 'group', value.name = 'AUC')
  
  markers_all <- list(markers_single = subset(raman_markers, AUC > threshold) )
  
  if(paired){
    cat('Finding paired markers ... \n')
    raman_markers <- calculatePairedMarkersAUC(matrix, group, threshold = threshold, batch_size = batch_size)
    raman_markers$col1 <- wave[raman_markers$col1]
    raman_markers$col2 <- wave[raman_markers$col2]
    raman_markers$group <- u_group[raman_markers$group]
    markers_all$markers_paired <- subset(raman_markers, AUC > threshold)
  }
  
  return(markers_all)
}