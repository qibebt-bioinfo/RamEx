library(testthat)
library(RamEx)

test_that("data loading works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  expect_true(inherits(ramnome_obj, "Ramanome"))
  expect_equal(nrow(ramnome_obj@datasets$raw.data), 500)
})


test_that("data smoothing works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  smoothed_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  expect_true(inherits(smoothed_obj, "Ramanome"))
  expect_equal(nrow(smoothed_obj@datasets$smooth.data), 500)
})

test_that("baseline correction works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  expect_true(inherits(baseline_obj, "Ramanome"))
  expect_equal(nrow(baseline_obj@datasets$baseline.data), 500)
})

test_that("data normalization works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
  expect_true(inherits(normalized_obj, "Ramanome"))
  expect_equal(nrow(normalized_obj@datasets$normalized.data), 500)
})


test_that("quality control works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
  cleaned_obj <- Qualitycontrol.ICOD(normalized_obj@datasets$normalized.data, var_tol = 0.4)
  expect_true(is.list(cleaned_obj))
  #expect_true(nrow(cleaned_obj$interations) <= 500)
})

test_that("PCA reduction works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
  cleaned_obj <- Qualitycontrol.ICOD(normalized_obj@datasets$normalized.data, var_tol = 0.4)
  cleaned_obj <- normalized_obj[cleaned_obj$quality,] 
  pca_obj <- Feature.Reduction.Pca(cleaned_obj, draw = FALSE, save = FALSE)
  expect_true(inherits(pca_obj, "Ramanome"))
  expect_equal(ncol(pca_obj@reductions$PCA), 2)
})

test_that("PCA plot generation works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
  cleaned_obj <- Qualitycontrol.ICOD(normalized_obj@datasets$normalized.data, var_tol = 0.4)
  cleaned_obj <- normalized_obj[cleaned_obj$quality,] 
  pca_obj <- Feature.Reduction.Pca(cleaned_obj, draw = TRUE, save = FALSE)
  expect_true(file.exists("Rplots.pdf"))
})

# test_that("clustering analysis works correctly", {
#   data(RamEx_data)
#   ramnome_obj <- RamEx_data
#   ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
#   baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
#   normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
#   cleaned_obj <- Qualitycontrol.ICOD(normalized_obj@datasets$normalized.data, var_tol = 0.4)
#   cleaned_obj <- normalized_obj[cleaned_obj$quality,] 
#   pca_obj <- Feature.Reduction.Pca(cleaned_obj, draw = FALSE, save = FALSE)
#   cluster_info <- Phenotype.Analysis.Louvaincluster(object = pca_obj, resolutions = c(0.8))
#   expect_true(is.data.frame(cluster_info))
# })

test_that("t-SNE reduction works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
  cleaned_obj <- Qualitycontrol.ICOD(normalized_obj@datasets$normalized.data, var_tol = 0.4)
  cleaned_obj <- normalized_obj[cleaned_obj$quality,] 
  pca_obj <- Feature.Reduction.Pca(cleaned_obj, draw = FALSE, save = FALSE)
  tsne_obj <- Feature.Reduction.Tsne(pca_obj, draw = TRUE, save = FALSE)
  expect_true(inherits(tsne_obj, "Ramanome"))
  expect_equal(ncol(tsne_obj@reductions$tSNE), 2)
})

test_that("UMAP reduction works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
  cleaned_obj <- Qualitycontrol.ICOD(normalized_obj@datasets$normalized.data, var_tol = 0.4)
  cleaned_obj <- normalized_obj[cleaned_obj$quality,] 
  pca_obj <- Feature.Reduction.Pca(cleaned_obj, draw = FALSE, save = FALSE)
  umap_obj <- Feature.Reduction.Umap(pca_obj, draw = TRUE, save = FALSE)
  expect_true(inherits(umap_obj, "Ramanome"))
  expect_equal(ncol(umap_obj@reductions$UMAP), 2)
})


test_that("UMAP plot generation works correctly", {
  data(RamEx_data)
  ramnome_obj <- RamEx_data
  ramnome_obj <- Preprocessing.Smooth.Sg(ramnome_obj)
  baseline_obj <- Preprocessing.Baseline.Polyfit(ramnome_obj)
  normalized_obj <- Preprocessing.Normalize(baseline_obj, "ch")
  cleaned_obj <- Qualitycontrol.ICOD(normalized_obj@datasets$normalized.data, var_tol = 0.4)
  cleaned_obj <- normalized_obj[cleaned_obj$quality,] 
  pca_obj <- Feature.Reduction.Pca(cleaned_obj, draw = FALSE, save = FALSE)
  umap_obj <- Feature.Reduction.Umap(pca_obj, draw = TRUE, save = TRUE)
  expect_true(file.exists("Reduction.umap.png"))
})
