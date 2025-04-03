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


#' Perform classification use linear discriminant analysis LDA model
#'
#' This function performs classification using PC-LDA
#'
#' @param train The training data object
#' @param test The test data object (optional). If not provided, the function will perform a stratified cross-validation (Training : Test = 7 : 3).
#' @param show Whether user want to show the confusion matrix plot of the results
#' @param save Wether user want to save the confusion matrix plot of the results (default path : getwd())
#' @param seed The random seed
#' @param n_pc The number of principal components to use
#' 
#' @return A list containing:
#' \describe{
#'   \item{model}{The LDA model}
#'   \item{pred_test}{The prediction for test data if test is provided}
#' }
#' @importFrom MASS lda
#' @importFrom ggplot2 ggsave
#' @export Classification.Lda
#' 
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' model.lda <- Classification.Lda(data_processed)
Classification.Lda <- function(train, test = NULL, show=TRUE, save=FALSE, seed=42, n_pc = 20) {
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
  
  data.pca <- prcomp(data_train, scale = TRUE, retx = TRUE)
  data_20 <- scale(data_train, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:n_pc] %>% as.data.frame
  
  model.lda <- lda(label_train ~ ., data = data_20)
  data_pre <- predict(model.lda, data_20)
  cat('Training accuracy of PC-LDA: ')
  pred.train <- confusion.plot(label_train, data_pre$class)
  if(save){
    cat('Saving plot to the current working directory: ', getwd(), '\n')
    ggsave('Classification_PC-LDA_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}
  
  test_20 <- scale(data_val, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:n_pc] %>% as.data.frame
  
  data_pre <- predict(model.lda, test_20)
  pca_params <- list(
    center = data.pca$center,
    scale = data.pca$scale,
    rotation = data.pca$rotation[, 1:n_pc]
  )
  cat('Test accuracy of PC-LDA: ')
  pred.test <- confusion.plot(label_val, data_pre$class)
  if(save){
    cat('Saving plot to the current working directory: ', getwd(), '\n')
    ggsave('Classification_PC-LDA_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  if(show){print(pred.train)
    print(pred.test)}
  if (is.null(test)) return(list(model=model.lda, pca_params = pca_params)) else return(list(model=model.lda, pred_test = data_pre, pca_params = pca_params))
}


#' Perform classification use support vector machine (SVM) model
#'
#' This function performs classification using SVM
#'
#' @param train The training data object
#' @param test The test data object (optional). If not provided, the function will perform a stratified cross-validation (Training : Test = 7 : 3).
#' @param show Whether user want to show the confusion matrix plot of the results
#' @param save Wether user want to save the confusion matrix plot of the results (default path : getwd())
#' @param seed The random seed
#' 
#' @return A list containing:
#' \describe{
#'   \item{model}{The SVM model}
#'   \item{pred_test}{The prediction for test data if test is provided}
#' }
#' @importFrom e1071 svm
#' @importFrom ggplot2 ggsave
#' @export Classification.Svm
#' 
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' model.svm <- Classification.Svm(data_processed)
Classification.Svm <- function(train, test = NULL, show=TRUE, save=FALSE, seed=42) {
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
  model.svm <- svm(x = data_train, y = as.factor(label_train), type = "C-classification", scale = TRUE, kernel = 'linear')
  cat('Training accuracy of SVM: ')
  pred.train <- confusion.plot(label_train, predict(model.svm, data_train))
  if(save){ggsave('Classification_SVM_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}
  
  cat('Test accuracy of SVM: ')
  data_pre <- predict(model.svm, data_val)
  pred.test <- confusion.plot(label_val, data_pre)
  if(save){
    cat('Saving plot to the current working directory: ', getwd(), '\n')
    ggsave('Classification_SVM_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  if(show){print(pred.train)
    print(pred.test)}
  if (is.null(test)) return(list(model=model.svm)) else return(list(model=model.svm, pred_test = data_pre))
}


#' Perform classification use RandomForest (RF) model
#'
#' This function performs classification using RF
#'
#' @param train The training data object
#' @param test The test data object (optional). If not provided, the function will perform a stratified cross-validation (Training : Test = 7 : 3).
#' @param ntree The number of trees in the forest
#' @param mtry The number of variables randomly sampled as candidates for splitting at each node
#' @param show Whether user want to show the confusion matrix plot of the results
#' @param save Wether user want to save the confusion matrix plot of the results (default path : getwd())
#' @param seed The random seed
#' 
#' @return A list containing:
#' \describe{ 
#'   \item{model}{The RF model}
#'   \item{pred_test}{The prediction for test data if test is provided}
#' }
#' @importFrom ggplot2 ggsave
#' @importFrom randomForest randomForest
#' @export Classification.Rf
#' 
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' model.rf <- Classification.Rf(data_processed)
Classification.Rf <- function(train, test = NULL, ntree = 100, mtry = 2, show=TRUE, save=FALSE, seed=42) {
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
  model.rf <- randomForest(data_train, as.factor(label_train), ntree = ntree, mtry = mtry, replace = TRUE)
  cat('Training accuracy of Random forest: ')
  pred.train <- confusion.plot(label_train, predict(model.rf, data_train))
  if(save){ggsave('Classification_RF_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}
  
  data_pre <- predict(model.rf, data_val)
  cat('Test accuracy of Random forest: ')
  pred.test <- confusion.plot(label_val, data_pre)
  if(save){
    cat('Saving plot to the current working directory: ', getwd(), '\n')
    ggsave('Classification_RF_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
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
#' @param test The test data object (optional). If not provided, the function will perform a stratified cross-validation (Training : Test = 7 : 3).
#' @param n_pc The number of principal components to use
#' 
#' @return A list containing:
#' \describe{ 
#'   \item{model}{The GMM model}
#'   \item{pred_test}{The prediction for test data if test is provided}
#' }
#' @export Classification.Gmm
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC

#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' model.gmm <- Classification.Gmm(data_processed)
Classification.Gmm <- function(train, test = NULL, n_pc = 20, show=TRUE, save=FALSE, seed=42) {
  set.seed(seed)
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
  data_train_20 <- scale(data_train, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:n_pc] %>% as.data.frame
  gmm_model <- Mclust(data_train_20)
  data_test_20 <- scale(data_val, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:n_pc] %>% as.data.frame
  data_pre <- predict(gmm_model, data_test_20)
  cat('Training accuracy of GMM: ')
  pred.train <- confusion.plot(label_train, predict(gmm_model, data_train_20)$class)
  pca_params <- list(
    center = data.pca$center,
    scale = data.pca$scale,
    rotation = data.pca$rotation[, 1:n_pc]
  )
  cat('Test accuracy of GMM: ')
  pred.test <- confusion.plot(label_val, data_pre$class)
  if(save){
    cat('Saving plot to the current working directory: ', getwd(), '\n')
    ggsave('Classification_GMM_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  if(show){print(pred.train)
    print(pred.test)}

  if (is.null(test)) return(list(model=gmm_model, pca_params = pca_params)) else return(list(model=gmm_model, pred_test = data_pre$class, pca_params = pca_params))}


#' Predict using a trained classification model
#'
#' This function uses a saved classification model to predict labels for new data.
#'
#' @param model The saved classification model (from Classification.Lda, Classification.Svm, Classification.Rf, or Classification.Gmm)
#' @param new_data The new Ramanome object to predict
#' @param show Whether to show the confusion matrix plot of the results
#' @param save Whether to save the confusion matrix plot of the results (default path: getwd())
#'
#' @return A list containing:
#' \describe{
#'   \item{predictions}{The predicted labels for the new data}
#'   \item{probabilities}{The prediction probabilities for each class (if available)}
#' }
#'
#' @export predict_classification
#'
#' @examples
#' data(RamEx_data)
#' data_processed <- Preprocessing.OneStep(RamEx_data)
#' # Train a model
#' model <- Classification.Rf(data_processed)
#' # Use the model to predict new data
#' predictions <- predict_classification(model, data_processed)
predict_classification <- function(model, new_data, show = TRUE, save = FALSE) {
  # Get the new data
  new_data_matrix <- get.nearest.dataset(new_data)
  
  # Determine model type and make predictions
  if (inherits(model$model, "lda")) {
    # For LDA model
    new_data_20 <- scale(new_data_matrix, center = model$pca_params$center, scale = model$pca_params$scale) %*% model$pca_params$rotation %>% as.data.frame
    pred <- predict(model$model, new_data_20)
    predictions <- pred$class
    probabilities <- pred$posterior
    
  } else if (inherits(model$model, "svm")) {
    # For SVM model
    pred <- predict(model$model, new_data_matrix, probability = TRUE)
    predictions <- pred
    probabilities <- attr(pred, "probabilities")
    
  } else if (inherits(model$model, "randomForest")) {
    # For Random Forest model
    pred <- predict(model$model, new_data_matrix, type = "response")
    predictions <- pred
    probabilities <- predict(model$model, new_data_matrix, type = "prob")
    
  } else if (inherits(model$model, "Mclust")) {
    # For GMM model
    new_data_20 <- scale(new_data_matrix, center = model$pca_params$center, scale = model$pca_params$scale) %*% model$pca_params$rotation %>% as.data.frame
    pred <- predict(model$model, new_data_20)
    predictions <- pred$classification
    probabilities <- pred$z
    
  } else {
    stop("Unsupported model type")
  }
  
  # If true labels are available, show confusion matrix
  if (!is.null(new_data@meta.data$group)) {
    true_labels <- new_data@meta.data$group
    if (show) {
      pred.plot <- confusion.plot(true_labels, predictions)
      print(pred.plot)
    }
    if (save) {
      pred.plot <- confusion.plot(true_labels, predictions)
      cat('Saving plot to the current working directory: ', getwd(), '\n')
      ggsave('Classification_Prediction.png', pred.plot, 
             width = length(unique(true_labels)) + 1, 
             height = length(unique(true_labels)))
    }
  }
  
  return(list(predictions = predictions, probabilities = probabilities))
}