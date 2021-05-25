bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv",row.names = 1))
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- utils::read.csv("small_test_data/cell_annotations_small.csv",row.names = 1)$x

test_that("Bisque deconvolution works",{
  bisque_model <- as.matrix(read.csv("test_models/bisque_model_small.csv",row.names = 1,check.names = FALSE))
  deconvolution <- deconvolute(bulk_small,bisque_model, method = "bisque", sc_object_small, cell_annotations_small)
  expect_equal(info = "rows of deconvolution equal to columns of signature (same celltypes in same order)", object = sort(colnames(deconvolution)), expected = sort(colnames(bisque_model)))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution)) , expected = sort(colnames(bulk_small)))

  check_result <- as.matrix(read.csv("test_results/bisque_result_small.csv",row.names = 1,check.names = FALSE))
  expect_equal(info = "deconvolution result is correct", object = deconvolution, expected = check_result)
})


test_that("MOMF deconvolution works",{
  momf_model <- as.matrix(read.csv("test_models/momf_model_small.csv",row.names = 1,check.names = FALSE))
  deconvolution <- deconvolute(bulk_small,momf_model, method = "momf", sc_object_small)
  expect_equal(info = "columns of deconvolution equal to columns of signature (same celltypes in same order)", object = sort(colnames(deconvolution)), expected = sort(colnames(momf_model)))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution)) , expected = sort(colnames(bulk_small)))

  check_result <- as.matrix(read.csv("test_results/momf_result_small.csv",row.names = 1,check.names = FALSE))
  expect_equal(info = "deconvolution result is correct", object = deconvolution, expected = check_result)
})

test_that("DWLS deconvolution works", {
  dwls_model <- as.matrix(read.csv("test_models/dwls_model_small.csv", row.names = 1,check.names = FALSE))
  deconvolution_dwls <- deconvolute(bulk_small,dwls_model, method = "dwls", sc_object_small,dwls_submethod="DampenedWLS")
  deconvolution_ols <- deconvolute(bulk_small,dwls_model, method = "dwls", sc_object_small,dwls_submethod="OLS")
  deconvolution_svr <- deconvolute(bulk_small,dwls_model, method = "dwls", sc_object_small,dwls_submethod="SVR")
  expect_equal(info = "rows of deconvolution for dwls equal to columns of signature (same celltypes, not same order)", object = sort(colnames(deconvolution_dwls)), expected = sort(colnames(dwls_model)))
  expect_equal(info = "deconvolution of dwls contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution_dwls)) , expected = sort(colnames(bulk_small)))
  expect_equal(info = "rows of deconvolution for ols equal to columns of signature (same celltypes, not same order)", object = sort(colnames(deconvolution_ols)), expected = sort(colnames(dwls_model)))
  expect_equal(info = "deconvolution of ols contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution_ols)) , expected = sort(colnames(bulk_small)))
  expect_equal(info = "rows of deconvolution for svr equal to columns of signature (same celltypes, not same order)", object = sort(colnames(deconvolution_svr)), expected = sort(colnames(dwls_model)))
  expect_equal(info = "deconvolution of svr contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution_svr)) , expected = sort(colnames(bulk_small)))

  check_result_dwls <- as.matrix(read.csv("test_results/dwls_dwls_result_small.csv",row.names = 1,check.names = FALSE))
  check_result_ols <- as.matrix(read.csv("test_results/dwls_ols_result_small.csv",row.names = 1,check.names = FALSE))
  check_result_svr <- as.matrix(read.csv("test_results/dwls_svr_result_small.csv",row.names = 1,check.names = FALSE))
  expect_equal(info = "deconvolution result for dwls is correct", object = deconvolution_dwls, expected = check_result_dwls)
  expect_equal(info = "deconvolution result for dwls is correct", object = deconvolution_ols, expected = check_result_ols)
  expect_equal(info = "deconvolution result for dwls is correct", object = deconvolution_svr, expected = check_result_svr)
})

test_that("Scaden deconvolution works", {
  model_dir <- paste0(tempdir(),"/model")
  skip_if_not(dir.exists(model_dir), message = "skipping deconvolution test")
  deconvolution <- deconvolute(bulk_small,model_dir, method = "scaden")
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(rownames(deconvolution)) , expected = sort(colnames(bulk_small)))
})

