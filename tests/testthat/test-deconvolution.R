bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv",row.names = 1))
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- utils::read.csv("small_test_data/cell_annotations_small.csv",row.names = 1)$x

test_that("Bisque deconvolution works",{
  bisque_model <- as.matrix(read.csv("test_models/bisque_model_small.csv",row.names = 1,check.names = FALSE))
  deconvolution <- deconvolute(bulk_small,bisque_model, method = "bisque", sc_object_small, cell_annotations_small)
  expect_equal(info = "rows of deconvolution equal to columns of signature (same celltypes in same order)", object = sort(colnames(deconvolution)), expected = sort(colnames(bisque_model)))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution)) , expected = sort(colnames(bulk_small)))
})


test_that("MOMF deconvolution works",{
  momf_model <- as.matrix(read.csv("test_models/momf_model_small.csv",row.names = 1,check.names = FALSE))
  deconvolution <- deconvolute(bulk_small,momf_model, method = "momf", sc_object_small)
  expect_equal(info = "columns of deconvolution equal to columns of signature (same celltypes in same order)", object = sort(colnames(deconvolution)), expected = sort(colnames(momf_model)))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution)) , expected = sort(colnames(bulk_small)))
})

test_that("DWLS deconvolution works", {
  dwls_model <- as.matrix(read.csv("test_models/dwls_model_small.csv", row.names = 1,check.names = FALSE))
  deconvolution <- deconvolute(bulk_small,dwls_model, method = "dwls", sc_object_small)
  colnames(deconvolution)<- gsub("."," ", colnames(deconvolution),fixed = T)
  expect_equal(info = "rows of deconvolution equal to columns of signature (same celltypes, not same order)", object = sort(colnames(deconvolution)), expected = sort(colnames(dwls_model)))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution)) , expected = sort(colnames(bulk_small)))
})

test_that("Scaden deconvolution works", {
  model_dir <- paste0(tempdir(),"/model")
  skip_if_not(dir.exists(model_dir), message = "skipping deconvolution test")
  deconvolution <- deconvolute(bulk_small,model_dir, method = "scaden")
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution)) , expected = sort(colnames(bulk_small)))
})

