#' Perform stratified partitioning of data indices
#'
#' This function performs stratified sampling to split data indices into training and test sets
#' while maintaining class proportions.
#'
#' @param labels A vector of class labels
#' @param p The proportion of data to include in the training set (default: 0.7)
#' @return A vector of indices for the training set


stratified_partition <- function(labels, p = 0.7) {
  unique_labels <- unique(labels)
  train_indices <- unlist(lapply(unique_labels, function(label) {
    label_indices <- which(labels == label)
    sample(label_indices, size = round(length(label_indices) * p), replace = FALSE)
  }))
  return(train_indices)
}


#' Perform classification use LDA model
#'
#' This function performs classification using PC-LDA
#'
#' @param train The training data object
#' @param test The test data object (optional)
#' @param show Whether user want to show the results
#' @param save Wether user want to save the results
#' @param seed The random seed
#' @return A list containing:
#'   \item{pred.train}{The confusion matrix plot for training data}
#'   \item{pred.test}{The confusion matrix plot for test data}
#' @importFrom MASS lda
#' @importFrom ggplot2 ggsave
#' @export Classification.Lda
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' Classification.Lda(data_cleaned)
Classification.Lda <- function(train, test = NULL, show=TRUE, save=FALSE, seed=42) {
  if (is.null(test)) {
    data_set <- get.nearest.dataset(train)
    labels <- train@meta.data$group
    index <- stratified_partition(labels, p = 0.7)
    data_train <- data_set[index,]
    label_train <- labels[index]
    data_val <- data_set[-index,]
    label_val <- labels[-index]
  } else {
    data_train <- get.nearest.dataset(train)
    label_train <- train@meta.data$group
    data_val <- get.nearest.dataset(test)
    label_val <- test@meta.data$group
  }
  
  data.pca <- prcomp(data_train, scale = TRUE, retx = TRUE)
  data_20 <- scale(data_train, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame
  
  model.lda <- lda(label_train ~ ., data = data_20)
  data_pre <- predict(model.lda, data_20)
  cat('Training accuracy of PC-LDA: ')
  pred.train <- confusion.plot(label_train, data_pre$class)
  if(save){ggsave('Classification_PC-LDA_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}
  
  test_20 <- scale(data_val, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame
  
  data_pre <- predict(model.lda, test_20)
  cat('Test accuracy of PC-LDA: ')
  pred.test <- confusion.plot(label_val, data_pre$class)
  if(save){ggsave('Classification_PC-LDA_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  if(show){print(pred.train)
    print(pred.test)}
  if (is.null(test)) return(list(model=model.lda)) else return(list(model=model.lda, pred_test = data_pre))
}


#' Perform classification use svm model
#'
#' This function performs classification using svm
#'
#' @param train The training data object
#' @param test The test data object (optional)
#' @param show Whether user want to show the results
#' @param save Wether user want to save the results
#' @param seed The random seed
#' @return A list containing:
#'   \item{pred.train}{The confusion matrix plot for training data}
#'   \item{pred.test}{The confusion matrix plot for test data}
#' @importFrom e1071 svm
#' @importFrom ggplot2 ggsave
#' @export Classification.Svm
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' Classification.Svm(data_cleaned)
Classification.Svm <- function(train, test = NULL, show=TRUE, save=FALSE, seed=42) {
  if (is.null(test)) {
    data_set <- get.nearest.dataset(train)
    labels <- train@meta.data$group
    index <- stratified_partition(labels, p = 0.7)
    data_train <- data_set[index,]
    label_train <- labels[index]
    data_val <- data_set[-index,]
    label_val <- labels[-index]
  } else {
    data_train <- get.nearest.dataset(train)
    label_train <- train@meta.data$group
    data_val <- get.nearest.dataset(test)
    label_val <- test@meta.data$group
  }
  model.svm <- svm(x = data_train, y = as.factor(label_train), type = "C-classification", scale = TRUE, kernel = 'linear')
  cat('Training accuracy of SVM: ')
  pred.train <- confusion.plot(label_train, predict(model.svm, data_train))
  if(save){ggsave('Classification_SVM_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}
  
  cat('Test accuracy of SVM: ')
  data_pre <- predict(model.svm, data_val)
  pred.test <- confusion.plot(label_val, data_pre)
  if(save){ggsave('Classification_SVM_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  if(show){print(pred.train)
    print(pred.test)}
  if (is.null(test)) return(list(model=model.svm)) else return(list(model=model.svm, pred_test = data_pre))
}


#' Perform classification use RandomForest model
#'
#' This function performs classification using RF
#'
#' @param train The training data object
#' @param test The test data object (optional)
#' @param show Whether user want to show the results
#' @param save Wether user want to save the results
#' @param seed The random seed
#' @return A list containing:
#'   \item{pred.train}{The confusion matrix plot for training data}
#'   \item{pred.test}{The confusion matrix plot for test data}
#' @importFrom ggplot2 ggsave
#' @importFrom randomForest randomForest
#' @export Classification.Rf
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' Classification.Rf(data_cleaned)
Classification.Rf <- function(train, test = NULL, show=TRUE, save=FALSE, seed=42) {
  set.seed(seed)
  if (is.null(test)) {
    data_set <- get.nearest.dataset(train)
    labels <- train@meta.data$group
    index <- stratified_partition(labels, p = 0.7)
    data_train <- data_set[index,]
    label_train <- labels[index]
    data_val <- data_set[-index,]
    label_val <- labels[-index]
  } else {
    data_train <- get.nearest.dataset(train)
    label_train <- train@meta.data$group
    data_val <- get.nearest.dataset(test)
    label_val <- test@meta.data$group
  }
  model.rf <- randomForest(data_train, as.factor(label_train), ntree = 100, mtry = 2, replace = TRUE)
  cat('Training accuracy of Random forest: ')
  pred.train <- confusion.plot(label_train, predict(model.rf, data_train))
  if(save){ggsave('Classification_RF_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}
  
  data_pre <- predict(model.rf, data_val)
  cat('Test accuracy of Random forest: ')
  pred.test <- confusion.plot(label_val, data_pre)
  if(save){ggsave('Classification_RF_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  if(show){print(pred.train)
    print(pred.test)}
  if (is.null(test)) return(list(model=model.rf)) else return(list(model=model.rf, pred_test = data_pre))
}

#' Gaussian Mixture Model (GMM)
#'
#' A probabilistic model that assumes data is generated from
#' a mixture of Gaussian distributions and assigns probabilities
#'  to each class
#'
#' @param train The training data object
#' @param test The test data object (optional)
#' @return The GMM model.
#' @export Classification.Gmm
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' Classification.Gmm(data_cleaned)
Classification.Gmm <- function(train, test = NULL) {
  if (is.null(test)) {
    data_set <- get.nearest.dataset(train)
    labels <- train@meta.data$group
    index <- stratified_partition(labels, p = 0.7)
    data_train <- data_set[index,]
    label_train <- as.numeric(labels[index])
    data_val <- data_set[-index,]
    label_val <- labels[-index]
  } else {
    data_train <- get.nearest.dataset(train)
    label_train <- train@meta.data$group
    data_val <- get.nearest.dataset(test)
    label_val <- test@meta.data$group}
  data.pca <- prcomp(data_train, scale = TRUE, retx = TRUE)
  data_train_20 <- scale(data_train, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame
  gmm_model <- Mclust(data_train_20)
  data_test_20 <- scale(data_val, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame
  data_pre <- predict(gmm_model, data_test_20)
  if (is.null(test)) return(list(model=gmm_model)) else return(list(model=gmm_model, pred_test = data_pre$classification))}