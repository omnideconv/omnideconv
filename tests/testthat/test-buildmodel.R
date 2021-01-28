bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv",row.names = 1))
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- utils::read.csv("small_test_data/cell_annotations_small.csv",row.names = 1)$x

test_that("Bisque GenerateSCReference works",{
  expect_equal(info = "signature matrix has same amount of rows as single cell matrix", object = nrow(build_model(sc_object_small, cell_annotations_small, method = "bisque")), expected = nrow(sc_object_small))
  expect_equal(info = "signature matrix has same amount of columns as unique cell types in single cell matrix", object = ncol(build_model(sc_object_small, cell_annotations_small, method = "bisque")), expected = length(unique(cell_annotations_small)))
})


test_that("MOMF compute reference works",{
  expect_equal(info = "signature matrix has same amount of rows as single cell matrix", object = nrow(build_model(sc_object_small, cell_annotations_small, method = "momf")), expected = nrow(sc_object_small))
  expect_equal(info = "signature matrix has same amount of columns as unique cell types in single cell matrix", object = ncol(build_model(sc_object_small, cell_annotations_small, method = "momf")), expected = length(unique(cell_annotations_small)))
})

test_that("DWLS build signature matrix works", {
  expect_equal(info = "signature matrix has same amount of columns as unique cell types in single cell matrix",object = ncol(build_model(sc_object_small, cell_annotations_small, method = "dwls")), expected = length(unique(cell_annotations_small)))
})

test_that("Scaden build model works", {
  model <- build_model(sc_object_small, cell_annotations_small, method="scaden", bulk_data = bulk_small, samples=10, cells=5, steps=150, verbose=T)
  expect_equal(info = "model folder is created and model assets are written", object= length(list.dirs(model, recursive = F)), expected = 3)
})
