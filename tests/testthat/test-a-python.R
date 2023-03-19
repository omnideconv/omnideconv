library(omnideconv)
# library(reticulate)
# reticulate::py_config()


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




bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv", row.names = 1))
# bulk_small <- bulk_small[, 1, drop = FALSE]
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- readr::read_lines("small_test_data/cell_annotations_small.txt")
batch_ids_small <- readr::read_lines("small_test_data/batch_ids_small.txt")
marker_genes <- readr::read_lines("small_test_data/marker_genes_small.txt")

markers_small <- list(marker_genes[1:9], marker_genes[10:14], marker_genes[15:20])
names(markers_small) <- sort(unique(cell_annotations_small))


test_that("Scaden build model works", {
  reticulate::use_python()
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
  reticulate::use_python(python = python)
  model <- build_model(sc_object_small, cell_annotations_small, "autogenes")
  expect_true(file.exists(model), "pickle file was created successfully")
})
