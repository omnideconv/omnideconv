library(omnideconv)
py_config()


test_that("Python environment creation works", {
  # init_python()
  expect_identical(info = "Python available", object = reticulate::py_available(), expected = TRUE)
})

test_that("Python module Scaden is available", {
  # scaden_checkload()
  expect_identical(
    info = "Scaden available", object = reticulate::py_module_available("scaden"),
    expected = TRUE
  )
})

test_that("Python module anndata is available", {
  # anndata_checkload()
  expect_identical(
    info = "Anndata available", object = reticulate::py_module_available("anndata"),
    expected = TRUE
  )
})
test_that("Python module autogenes is available", {
  # autogenes_checkload()
  expect_identical(
    info = "AutoGeneS available",
    object = reticulate::py_module_available("autogenes"), expected = TRUE
  )
})



test_that("Scaden build model works", {
  model <- build_model(sc_object_small, cell_annotations_small,
    method = "scaden",
    bulk_gene_expression = bulk_small, samples = 10, cells = 5,
    steps = 150, verbose = F
  )
  expect_equal(
    info = "model folder is created and model assets are written",
    object = length(list.dirs(model, recursive = F)), expected = 3
  )
})

test_that("AutoGeneS build model works", {
  model <- build_model(sc_object_small, cell_annotations_small, "autogenes")
  expect_true(file.exists(model), "pickle file was created successfully")
})
