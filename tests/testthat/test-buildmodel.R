load("../../test_data/bulk.RData")
load("../../test_data/cell_type_annotations.RData")
load("../../test_data/single_cell_data.RData")

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
  model <- build_model(sc_object_small, cell_annotations_small, method="scaden", bulk_data = bulk_small, samples=10, cells=5, steps=150)
  expect_equal(info = "model folder is created and model assets are written", object= length(list.dirs(model, recursive = F)), expected = 3)
})
