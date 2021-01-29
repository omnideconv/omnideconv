bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv",row.names = 1))
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- utils::read.csv("small_test_data/cell_annotations_small.csv",row.names = 1)$x

test_that("Bisque GenerateSCReference works",{
  signature <- build_model(sc_object_small, cell_annotations_small, method = "bisque")
  expect_equal(info = "signature matrix has same amount of rows as single cell matrix", object = nrow(signature), expected = nrow(sc_object_small))
  expect_equal(info = "signature matrix has same amount of columns as unique cell types in single cell matrix", object = ncol(signature), expected = length(unique(cell_annotations_small)))
  check_signature <- as.matrix(read.csv("test_models/bisque_model_small.csv",row.names = 1))
  colnames(check_signature) <- gsub("."," ", colnames(check_signature),fixed = T)
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})


test_that("MOMF compute reference works",{
  signature <- build_model(sc_object_small, cell_annotations_small, method = "momf", bulk_gene_expression = bulk_small)
  expect_equal(info = "signature matrix has same amount of rows as single cell matrix", object = nrow(signature), expected = nrow(sc_object_small))
  expect_equal(info = "signature matrix has same amount of columns as unique cell types in single cell matrix", object = ncol(signature), expected = length(unique(cell_annotations_small)))
  check_signature <- as.matrix(read.csv("test_models/momf_model_small.csv",row.names = 1))
  colnames(check_signature)<- gsub("."," ", colnames(check_signature),fixed = T)
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})

test_that("DWLS build signature matrix works", {
  signature <- build_model(sc_object_small, cell_annotations_small, method = "dwls")
  expect_equal(info = "signature matrix has same amount of columns as unique cell types in single cell matrix",object = ncol(signature), expected = length(unique(cell_annotations_small)))
  check_signature <- as.matrix(read.csv("test_models/dwls_model_small.csv", row.names = 1))
  colnames(check_signature)<- gsub("."," ", colnames(check_signature),fixed = T)
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})

test_that("Scaden build model works", {
  model <- build_model(sc_object_small, cell_annotations_small, method="scaden", bulk_data = bulk_small, samples=10, cells=5, steps=150)
  expect_equal(info = "model folder is created and model assets are written", object= length(list.dirs(model, recursive = F)), expected = 3)
})
