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
#' @examples
#' # Create sample Ramanome training object
#' wavenumbers <- seq(500, 3500, length.out=100)
#' train_spectra <- matrix(rnorm(1000), nrow=10)
#' train_metadata <- data.frame(
#'   group = factor(rep(c("Low", "High"), each=5))
#' )
#' train_obj <- new("Ramanome",
#'                  datasets=list(raw.data=train_spectra),
#'                  wavenumber=wavenumbers,
#'                  meta.data=train_metadata)
#'
#' # Create sample test object
#' test_spectra <- matrix(rnorm(500), nrow=5)
#' test_metadata <- data.frame(
#'   group = factor(rep(c("Low", "High"), c(3,2)))
#' )
#' test_obj <- new("Ramanome",
#'                 datasets=list(raw.data=test_spectra),
#'                 wavenumber=wavenumbers,
#'                 meta.data=test_metadata)
#'
#' # Perform PLS prediction
#' \dontrun{
#' pls_result <- Quantification.Pls(train_obj, test_obj)
#' print(pls_result)
#' }
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
#' @examples
#' # Create sample Ramanome training object with concentration data
#' wavenumbers <- seq(500, 3500, length.out=100)
#' train_spectra <- matrix(rnorm(1000), nrow=10)
#' train_metadata <- data.frame(
#'   group = rnorm(10, mean=50, sd=10) # Concentration values
#' )
#' train_obj <- new("Ramanome",
#'                  datasets=list(raw.data=train_spectra),
#'                  wavenumber=wavenumbers,
#'                  meta.data=train_metadata)
#'
#' # Create sample test object
#' test_spectra <- matrix(rnorm(500), nrow=5)
#' test_metadata <- data.frame(
#'   group = rnorm(5, mean=50, sd=10)
#' )
#' test_obj <- new("Ramanome",
#'                 datasets=list(raw.data=test_spectra),
#'                 wavenumber=wavenumbers,
#'                 meta.data=test_metadata)
#'
#' # Perform MLR prediction
#' \dontrun{
#' mlr_result <- Quantification.Mlr(train_obj, test_obj)
#' print(mlr_result)
#' }
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
#' @examples
#' # Create sample Ramanome training object with count data
#' wavenumbers <- seq(500, 3500, length.out=100)
#' train_spectra <- matrix(rnorm(1000), nrow=10)
#' train_metadata <- data.frame(
#'   group = rpois(10, lambda=5) # Count data following Poisson distribution
#' )
#' train_obj <- new("Ramanome",
#'                  datasets=list(raw.data=train_spectra),
#'                  wavenumber=wavenumbers,
#'                  meta.data=train_metadata)
#'
#' # Create sample test object
#' test_spectra <- matrix(rnorm(500), nrow=5)
#' test_metadata <- data.frame(
#'   group = rpois(5, lambda=5)
#' )
#' test_obj <- new("Ramanome",
#'                 datasets=list(raw.data=test_spectra),
#'                 wavenumber=wavenumbers,
#'                 meta.data=test_metadata)
#'
#' # Perform GLM prediction
#' \dontrun{
#' glm_result <- Quantification.Glm(train_obj, test_obj)
#' print(glm_result)
#' }
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
