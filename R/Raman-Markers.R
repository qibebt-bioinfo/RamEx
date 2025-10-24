#' Calculate mislabel scores for classification results
#'
#' This function calculates mislabel scores based on class probabilities for classification results.
#'
#' @param y A vector containing the true labels.
#' @param y.prob A matrix containing class probabilities.
#'
#' @return A matrix containing mislabel scores for each class.
#' @importFrom stats model.matrix
#' @noRd


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
#' @noRd


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
#' @noRd

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
#' @noRd

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
#' @noRd
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
#' @noRd

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
#' @noRd
#' 
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
#' @noRd

save.rf.results.summary <- function(result, rf.opts, filename = 'summary.xls', outdir) {
  err <- mean(result$errs)
  err.sd <- sd(result$errs)
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
#' @noRd
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
#' @noRd

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
#' @noRd

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
#' @noRd

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
#' @noRd
#' 
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
#' @param threshold The threshold value for RBCS feature importances.
#' @param ntree The number of trees to build the random forest model.
#' @param errortype The error type used in random forest model evaluation.
#' @param verbose Logical value indicating whether to display progress messages.
#' @param nfolds The number of folds used in cross-validation.
#' @param rf.cv Logical value indicating whether to perform random forest with cross-validation.
#' @param show Logical value indicating whether to show the heatmap.
#' @param save Logical value indicating whether to draw the heatmap.
#'
#' @return A data frame containing the selected feature importances.
#' @importFrom pheatmap pheatmap
#' @importFrom dplyr %>%
#' @importFrom grDevices pdf
#' @export Raman.Markers.Rbcs
#' 
#' @details
#' more details about RBCS could be found in the reference below.
#' Teng L., Wang X.,  Wang X.,  Gou H.,  Ren L., & Wang T., 2016. Label-free, rapid and quantitative phenotyping of stress response in e. coli via ramanome. *Scientific Reports* 
#' 
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' Markers_RBCS <- Raman.Markers.Rbcs(data_processed, threshold = 0.003, show = TRUE)


Raman.Markers.Rbcs <- function(
    object,
    outfolder = "RF_imps",
    threshold = 0.002,
    ntree = 1000,
    errortype = 'oob',
    verbose = FALSE,
    nfolds = 3,
    rf.cv = FALSE,
    show = TRUE,
    save = FALSE
) {
  x <- get.nearest.dataset(object)
  y <- object@meta.data$group
  
  # Out-of-bag classification in training data and prediction in test data
  y <- as.factor(y)
  oob.result <- rf.out.of.bag(x, y, verbose = verbose, ntree = ntree)
  
  # RF importances of features  
  imps <- oob.result$importances
  imps <- imps[order(imps, decreasing = TRUE)]
  imps.cutoff <- imps[which(imps > threshold)]
  
  
  
  data.select <- x[, names(imps.cutoff)]
  means <- aggregate(data.select, by = list(y), mean)
  rownames(means) <- means[, 1]
  means <- means[, -1]
  colnames(means) <-  as.numeric(colnames(means)) %>% round(., 0)
  if (save) {
    cat('Saving RBCS importances of Raman features and heatmap to the current working directory: ', getwd(), '\n')
    dir.create(outfolder)
    sink(paste0(outfolder, "/RBCS_imps", "_All_sorted.xls"))
    cat("\t")
    write.table(imps, quote = FALSE, sep = "\t")
    sink()
  
    sink(paste0(outfolder, "/RBCS_imps", "_cutoff_sorted.xls"))
    cat("\t")
    write.table(imps.cutoff, quote = FALSE, sep = "\t")
    sink()
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
  if (show) {
    pheatmap::pheatmap(
      means[, -1],
      show_colnames = TRUE,
      angle_col = 45,
      scale = "none",
      cluster_cols = TRUE,
      cluster_rows = TRUE
    )
  }
  return(imps.cutoff)
}



#' Compute Correlations with a Target Variable
#'
#' Calculates the correlations between each column of a dataset X and a target variable y by traversing all possible single and paired markers.
#'
#' @param object A Ramanome object.
#' @param group A factor vector (should be the same length as the cell number in the object), as the target variable. default is object$group.
#' @param paired Whether to consider paired markers, e.g. wave1 / wave2. Default is FALSE
#' @param min.cor The minimum correlation threshold. Defaults to 0.8.
#' @param by_average Whether to average the spectra of each group. Default is FALSE. if TRUE, the spectra will be execute baseline correction in peak regions to remove the baseline caused by the local, adjacent peak signal.
#' @param min.range The minimum wavenumber range of the bands when bubble method considering them as a peak (see 'by_average'). Default is 30.
#' @param extract_num Whether to extract numeric values from the group column. Default is TRUE. If else, the correlations will be calculated based on the group factor level.
#' @return A list containing the high correlation elements and their corresponding correlations with the target variable, containing two data frames: 'correlations_singular' for single markers (group ~ wave) and 'combination_correlations_two' for paired markers (group ~ wave1 / wave2).
#' 
#' @export Raman.Markers.Correlations
#' @importFrom dplyr filter
#' @importFrom stats cor
#' @importFrom utils combn
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel parLapply
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' options(mc.cores = 2)
#' Markers_corr <- Raman.Markers.Correlations(preprocessing.cutoff(data_processed, 800, 1300), min.cor = 0.6)
Raman.Markers.Correlations <- function(object,group=object$group, paired = FALSE,
min.cor=0.8, by_average = FALSE, min.range = 30,  extract_num = TRUE) {
  if(by_average){
    aver_spec <- average_spectra(get.nearest.dataset(object), group, min.range)
    X <- as.data.frame(aver_spec$means)
    group <- aver_spec$group
  }else{
    X <- as.data.frame(get.nearest.dataset(object))
  }
  waves <- object@wavenumber
  y <- if(extract_num) as.numeric(gsub("[^0-9]", "", group)) else as.numeric(group)
  
  correlations <- sapply(X, function(col) cor(col, y, use = "complete.obs"))
  
  cat('Finding singular markers ... \n')
  high_cor_elements <- which(abs(correlations) > min.cor)
  
  if(paired){
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
  stopCluster(cl)
  }else{
    combination_correlations_two <- NULL
  }
  
  list(
    markers_singular = if(length(high_cor_elements) == 0) NA else data.frame(wave=object@wavenumber[high_cor_elements], corrlations = correlations[high_cor_elements]),
    markers_paired = if(length(combination_correlations_two) == 0) NA else combination_correlations_two
  )
}


#' Average spectra of each group
#'
#' This function averages the spectra of each group in a Ramanome object.
#'
#' @param matrix A matrix of Raman spectra.
#' @param group A factor indicating the grouping variable.
#' @param min.range The minimum wavenumber range of the bands when bubble method considering them as a peak. 
#' @noRd
average_spectra <- function(matrix, group, min.range ){
  means <- aggregate(matrix, by = list(group), mean)
  group <- means[, 1]
  means <- as.matrix(means[, -1])
  bands <- apply(means, 1, function(x)bubblefill(x,min_bubble_widths = min.range)$bands)
  bands <- bands[[which.max(sapply(bands,length))]]
  bands[[1]][colnames(means)[2]] <- 1
  bands[[length(bands)]][colnames(means)[ncol(means)]] <- ncol(means)
  corre_means <- lapply(bands,function(x)polyfit(means[,c(x[1]:x[2])], degree = 2 )$corrected)
  corre_means <- do.call(cbind, corre_means)
  colnames(corre_means) <- colnames(means)
  return(list(means = corre_means, group = group))
}

#' Find Markers Using receiver operator characteristic (ROC) curve
#'
#' Performs receiver operator characteristic curve (ROC) analysis to identify significant markers by traversing all possible single and paired markers, and calculating the Area Under the Curve (AUC) for them, thus identify Raman features that best discriminate between different groups.
#' Each group/sample will get their own marker list.
#' 
#' @param object A Ramanome object.
#' @param threshold The minimum AUC threshold for considering a marker significant. Defaults to 0.75.
#' @param paired Whether to consider paired markers, e.g. wave1 / wave2. Default is FALSE
#' @param batch_size Sample size for each batch when traversing all possible paired markers. Default is 1,000
#' @return A list containing:
#' \describe{
#'   \item{markers_single}{A data frame of single wavelength markers with their AUC scores}
#'   \item{markers_paired}{A data frame of paired wavelength markers with their AUC scores (only if paired=TRUE)}
#' }
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
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' Markers_ROC <- Raman.Markers.Roc(preprocessing.cutoff(data_processed, 800, 1300), threshold = 0.75, paired = FALSE, batch_size = 1000)
Raman.Markers.Roc <- function(object, threshold = 0.75, paired = FALSE, batch_size = 1000, n_threads = NULL){
  matrix <- get.nearest.dataset(object)
  group <- object@meta.data$group
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
  
  markers_all <- list(markers_singular = subset(raman_markers, AUC > threshold) )
  
  if(paired){
    cat('Finding paired markers ... \n')
    cat('This may take a while, please be patient ... \n')
    print(ifelse(is.null(n_threads),0,n_threads))
    raman_markers <- calculatePairedMarkersAUC(matrix, group, threshold = threshold, batch_size = batch_size, n_threads = ifelse(is.null(n_threads),0,n_threads))
    raman_markers$col1 <- wave[raman_markers$col1]
    raman_markers$col2 <- wave[raman_markers$col2]
    raman_markers$group <- u_group[raman_markers$group]
    markers_all$markers_paired <- subset(raman_markers, AUC > threshold)
  }
  
  return(markers_all)
}

#' Find Markers Using Discriminant Analysis of Principal Components (DAPC)
#'
#' DAPC is a statistical method for classifying groups and identifying features that discriminate between them.
#' Principal Component Analysis (PCA) is firstly used to reduce data dimensionality and remove correlations among variables.
#' Then for each group, a one-vs-rest Linear Discriminant Analysis (LDA) is performed on the selected PCs to find the linear combinations that best separate the predefined groups.
#' And the loadings are back-projected to the original feature space to trace variable importance. The top N features most important for discriminating each group are returned.
#' 
#' @param object A Ramanome object.
#' @param n_pc  Integer. The number of principal components to retain for LDA. Default is 20.
#' @param top_n Integer. For each group, how many top variables to return according to absolute importance. Default is 100.
#' @return A data frame with the following columns:
#' \describe{
#'   \item{group}{The group label being discriminated}
#'   \item{wave}{The wavelength}
#'   \item{importance}{The calculated importance (projected LDA coefficient)}
#' }
#' Features of each group are sorted descendingly by absolute importance, and only the top N are reported.
#' @details
#' For each unique group label, a one-vs-rest LDA is trained on the PCA-reduced data, 
#' then the LDA coefficients are reconstructed back into the original feature space to rank feature importances.
#' @export Raman.Markers.DAPC
#' @importFrom irlba prcomp_irlba
#' @importFrom MASS lda
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' Markers_DAPC <- Raman.Markers.DAPC(preprocessing.cutoff(data_processed, 800, 1300), n_pc = 15, top_n = 40)
Raman.Markers.DAPC <- function(object, n_pc = 20, top_n = 100) {
  dat <- get.nearest.dataset(object)
  group <- object@meta.data$group
  # Step 1. PCA
  pca_res <- prcomp_irlba(dat, n = n_pc, center = TRUE, scale. = TRUE)
  X_pca <- pca_res$x[, 1:n_pc, drop=FALSE]
  out_list <- list()
  grp_all <- unique(group)
  
  for (g in grp_all) {
    # Step 2. One-vs-rest LDA
    y_bin <- factor(ifelse(group == g, as.character(g), "Other"))
    lda_res <- MASS::lda(X_pca, grouping = y_bin)
    # Step 3. Trace back the contribution of the original variable
    # lda_res$scaling: n_pc x 1
    coeffs_raw <- pca_res$rotation[, 1:n_pc, drop=FALSE] %*% as.matrix(lda_res$scaling)
    rownames(coeffs_raw) <- colnames(dat)
    # Step 4. Sort & Extract the top variable
    idx <- order(abs(coeffs_raw[,1]), decreasing=TRUE)[1:top_n]
    out_list[[as.character(g)]] <- data.frame(
      group = rep(as.character(g), top_n),
      wave = rownames(coeffs_raw)[idx],
      importance = coeffs_raw[idx, 1],
      stringsAsFactors = FALSE
    )
  }
  out_all <- do.call(rbind, out_list)
  rownames(out_all) <- NULL
  return(out_all)
}
