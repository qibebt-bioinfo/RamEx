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
#' @importFrom pls plsr
#' @importFrom stringr str_extract
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' quan_pls <- Quantification.Pls(data_cleaned)
Quantification.Pls <- function(train, test = NULL) {
  if (is.null(test)) {
    data_set <- train@datasets$normalized.data
    labels <- train@meta.data$group
    index <- stratified_partition(labels, p = 0.7)
    labels <- str_extract(labels, "\\d+")
    labels <- as.numeric(labels)
    data_train <- data_set[index,]
    data_train <- as.data.frame(data_train)
    data_train <- as.matrix(data_train)
    label_train <- labels[index]
    data_val <- data_set[-index,]
    data_val <- as.data.frame(data_val)
    data_val <- as.matrix(data_val)
    label_val <- labels[-index]
  } else {
    data_train <- train@datasets$normalized.data
    label_train <- train@meta.data$group
    label_train <- str_extract(label_train, "\\d+")
    label_train <- as.numeric(label_train)
    data_train <- as.data.frame(data_train)
    data_train <- as.matrix(data_train)
    data_val <- test@datasets$normalized.data
    data_val <- as.data.frame(data_val)
    data_val <- as.matrix(data_val)
    label_val <- test@meta.data$group
    label_val <- str_extract(label_val, "\\d+")
    label_val <- as.numeric(label_val)
  }
  pls_model <- plsr(label_train ~ data_train, ncomp = 8, scale = TRUE, validation = "none")
  pre_result <- predict(pls_model, data_val, ncomp = 8)
  if (is.null(test)) return(list(model=pls_model)) else return(list(model=pls_model, pred_test = pre_result))
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
#' @importFrom stringr str_extract
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' data_cleaned <- Feature.Reduction.Pca(data_cleaned,  draw=FALSE , save = FALSE)
#' quan_mlr <- Quantification.Mlr(data_cleaned)
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
  label_train <- str_extract(label_train, "\\d+")
  label_train <- as.numeric(label_train)
  label_val <- str_extract(label_val, "\\d+")
  label_val <- as.numeric(label_val)

  data.pca <- prcomp(data_train, scale = TRUE, retx = TRUE)
  data_20 <- scale(data_train, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame
  test_20 <- scale(data_val, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame

  mlr_model <- lm(label_train~ ., data = data_20)
  pre_result <- predict(mlr_model, test_20)
  if (is.null(test)) return(list(model=mlr_model)) else return(list(model=mlr_model, pred_test = pre_result))
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
#' @importFrom stringr str_extract
#' @examples
#' data(RamEx_data)
#' data_smoothed <- Preprocesssing.Smooth.Sg(RamEx_data)
#' data_baseline <- Preprocesssing.Baseline.Polyfit(data_smoothed)
#' data_baseline_bubble <- Preprocesssing.Baseline.Bubble(data_smoothed)
#' data_normalized <- Preprocesssing.Normalize(data_baseline, "ch")
#' data_cleaned <- Qualitycontrol.ICOD(data_normalized@datasets$normalized.data,var_tol = 0.4)
#' data_cleaned <- data_normalized[data_cleaned$index_good,]
#' data_cleaned <- Feature.Reduction.Intensity(data_cleaned, list(c(2000,2250),c(2750,3050), 1450, 1665))
#' data_cleaned <- Feature.Reduction.Pca(data_cleaned, draw=FALSE , save = FALSE)
#' quan_glm <- Quantification.Glm(data_cleaned)
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
  label_train <- str_extract(label_train, "\\d+")
  label_train <- as.numeric(label_train)
  label_val <- str_extract(label_val, "\\d+")
  label_val <- as.numeric(label_val)

  data.pca <- prcomp(data_train, scale = TRUE, retx = TRUE)
  data_20 <- scale(data_train, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame
  test_20 <- scale(data_val, center = data.pca$center, scale = data.pca$scale) %*% data.pca$rotation[, 1:20] %>% as.data.frame

  glm_model <- glm(label_train~ ., data_20, family = gaussian())
  pre_result <- predict(glm_model, test_20,  type = "response")
  if (is.null(test)) return(list(model=glm_model)) else return(list(model=glm_model, pred_test = pre_result))
}
