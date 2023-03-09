# test_that("Python environment creation works", {
#   init_python()
#   expect_identical(info = "Python available", object = reticulate::py_available(), expected = TRUE)
# })

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
