library(testthat)
library(RamEx)

test_that("read-spec function works correctly", {
  test_file_path <- system.file("extdata", "data", package = "RamEx")
  expect_true(file.exists(test_file_path))

  #expect_type(read.spec, "function")

  ramnome_obj <- read.spec.load(test_file_path)
  #expect_type(ramnome_obj, "Ramanome")

  expect_true(inherits(ramnome_obj, "Ramanome"))
  expect_equal(nrow(ramnome_obj@datasets$raw.data), 1215)


  #expect_error(read-spec("non_existent_file.txt"), class = "specific_error_class")

  #bad_format_file_path <- system.file("extdata", "bad_format_file.txt", package = "yourpackage")
 # expect_error(read-spec(bad_format_file_path), class = "specific_error_class")

})
