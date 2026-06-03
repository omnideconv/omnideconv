library(omnideconv)
# library(reticulate)
# reticulate::py_config()

omnideconv::install_all_python()

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

test_that("Rectangle conda environment is set up and rectanglepy is importable", {
  omnideconv::install_rectangle_python()
  expect_true(
    info = "r-omnideconv-rectangle conda env exists",
    "r-omnideconv-rectangle" %in% reticulate::conda_list()$name
  )
  python_bin <- omnideconv:::rectangle_python()
  exit_code <- system2(
    python_bin, c("-c", shQuote("import rectanglepy")),
    stdout = FALSE, stderr = FALSE
  )
  expect_equal(
    info = "rectanglepy is importable in the rectangle conda env",
    object = exit_code, expected = 0L
  )
})
