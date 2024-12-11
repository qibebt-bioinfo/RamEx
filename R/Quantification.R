#' Partial Least Squares (PLS)
#'
#' A regression method that reduces predictors to a smaller set of latent
#' variables while maximizing the covariance between predictors and response
#'  variables, ideal for situations where predictors are highly collinear
#'  or when the number of predictors exceeds the number of observations
#'
#' @param train The training data object
#' @param test The test data object (optional)
#' @return The PLS model.
#' @export Quantification.Pls
#' @importFrom mdatools pls
Quantification.Pls <- function(train, test = NULL) {
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
    label_val <- test@meta.data$group
  }
  pls_model <- pls(data_train,label_train)
  summary(pls_model)
  pre_result <- predict(pls_model, data_val)
  return(pre_result)
}



#' Multiple linear regression (MLR)
#' Simple and interpretable regression method, but assumes
#' no multicollinearity among predictors and a linear
#' relationship between predictors and response
#' @param train The training data object
#' @param test The test data object (optional)
#' @return The MLR prediction result.
#' @export Quantification.Mlr
#' @importFrom stats lm
Quantification.Mlr <- function(train, test = NULL) {
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
  mlr_model <- lm(data_train ~ label_train)
  summary(mlr_model)
  pre_result <- predict(mlr_model, as.data.frame(data_val))
  return(pre_result)
}


#' Generalized linear model (GLM)
#' SExtends linear regression by allowing the dependent
#' variable to follow distributions other than
#' normal (e.g., binomial, Poisson) and uses a link function
#'  to relate predictors to the response
#'
#' @param train The training data object
#' @param test The test data object (optional)
#' @return The GLM prediction result.
#' @export Quantification.Glm
#' @importFrom stats glm
Quantification.Glm <- function(train, test = NULL) {
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
  glm_model <- glm(data_train ~ label_train)
  summary(glm_model)
  pre_result <- predict(glm_model, as.data.frame(data_val))
  return(pre_result)
}
