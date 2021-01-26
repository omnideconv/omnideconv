load("../../test_data/bulk.RData")
load("../../test_data/cell_type_annotations.RData")
load("../../test_data/single_cell_data.RData")

bisque_model <- as.matrix(read.csv("test_models/bisque_model_small.csv",row.names = 1))
momf_model <- as.matrix(read.csv("test_models/momf_model_small.csv",row.names = 1))
dwls_model <- as.matrix(read.csv("test_models/dwls_model_small.csv", row.names = 1))
scaden_model <- "test_models/scaden_model_small"
colnames(bisque_model)<- gsub("."," ", colnames(bisque_model),fixed = T)
colnames(momf_model)<- gsub("."," ", colnames(momf_model),fixed = T)
colnames(dwls_model)<- gsub("."," ", colnames(dwls_model),fixed = T)


test_that("Bisque deconvolution works",{
  deconvolution <- deconvolute(bulk_small,bisque_model, method = "bisque", sc_object_small, cell_annotations_small)
  expect_equal(info = "rows of deconvolution equal to columns of signature (same celltypes in same order)", object = rownames(deconvolution), expected = colnames(bisque_model))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(colnames(deconvolution)) , expected = sort(colnames(bulk_small)))
})


test_that("MOMF deconvolution works",{
  deconvolution <- deconvolute(bulk_small,momf_model, method = "momf", sc_object_small)
  expect_equal(info = "rows of deconvolution equal to columns of signature (same celltypes in same order)", object = rownames(deconvolution), expected = colnames(bisque_model))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(colnames(deconvolution)) , expected = sort(colnames(bulk_small)))
})

test_that("DWLS deconvolution works", {
  deconvolution <- deconvolute(bulk_small,dwls_model, method = "dwls", sc_object_small)
  expect_equal(info = "rows of deconvolution equal to columns of signature (same celltypes in same order)", object = rownames(deconvolution), expected = colnames(bisque_model))
  expect_equal(info = "deconvolution contains same samples as in bulk (not same order)", object =  sort(colnames(deconvolution)) , expected = sort(colnames(bulk_small)))
})

test_that("Scaden build model works", {
  expect_equal(info = "model folder is created and model assets are written", object= length(list.dirs(model, recursive = F)), expected = 3)
})

