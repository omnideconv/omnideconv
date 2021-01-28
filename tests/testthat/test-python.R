test_that("Python environment creation works",{
  init_python()
  expect_identical(info = "Python available", object = reticulate::py_available(), expected = T)
})

test_that("Python module Scaden is available",{
  scaden_checkload()
  expect_identical(info = "Scaden available", object = reticulate::py_module_available("scaden"), expected = T)
})

test_that("Python module anndata is available",{
  scaden_checkload()
  expect_identical(info = "Anndata available", object = reticulate::py_module_available("anndata"), expected = T)
})
