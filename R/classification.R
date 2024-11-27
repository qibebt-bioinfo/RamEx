#' Perform classification
#'
#' This function performs classification using various methods such as PC-LDA, SVM, and Random Forest.
#'
#' @param train The training data object
#' @param test The test data object (optional)
#' @param method The classification method to use (default: 'PC-LDA')
#' @param show Whether user want to show the results
#' @param save Wether user want to save the results
#' @param seed The random seed
#' @importFrom MASS lda
#' @importFrom e1071 svm
#' @importFrom randomForest randomForest
#' @importFrom ggplot2 ggsave
#'
#' @export classification
classification <- function(train, test = NULL, method = 'PC-LDA', show=TRUE, save=FALSE, seed=42) {
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

  if (method == 'PC-LDA') {
    data.pca <- prcomp(data_train, scale = TRUE, retx = T)
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
  }

  if (method == 'SVM') {
    model.svm <- svm(x = data_train, y = as.factor(label_train), type = "C-classification", scale = TRUE, kernel = 'linear')
    cat('Training accuracy of SVM: ')
    pred.train <- confusion.plot(label_train, predict(model.svm, data_train))
    if(save){ggsave('Classification_SVM_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}

    cat('Test accuracy of SVM: ')
    pred.test <- confusion.plot(label_val, predict(model.svm, data_val))
    if(save){ggsave('Classification_SVM_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  }

  if (method == 'Random forest') {
    model.rf <- randomForest(data_train, as.factor(label_train), ntree = 100, mtry = 2, replace = TRUE)
    cat('Training accuracy of Random forest: ')
    pred.train <- confusion.plot(label_train, predict(model.rf, data_train))
    if(save){ggsave('Classification_RF_Train.png', pred.train, width = length(unique(label_train)) + 1, height = length(unique(label_train)))}

    data_pre <- predict(model.rf, data_val)
    cat('Test accuracy of Random forest: ')
    pred.test <- confusion.plot(label_val, data_pre)
    if(save){ggsave('Classification_RF_Test.png', pred.test, width = length(unique(label_val)) + 1, height = length(unique(label_train)))}
  }
  if(show){print(pred.train)
    print(pred.test)}
}

stratified_partition <- function(labels, p = 0.7) {
  unique_labels <- unique(labels)
  train_indices <- unlist(lapply(unique_labels, function(label) {
    label_indices <- which(labels == label)
    sample(label_indices, size = round(length(label_indices) * p), replace = FALSE)
  }))
  return(train_indices)
}
