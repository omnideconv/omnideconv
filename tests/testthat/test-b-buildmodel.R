library(omnideconv)

bulk_small <- system.file("small_test_data", "bulk_small.csv",
                          package = "omnideconv", mustWork = TRUE) %>%
  as.matrix(utils::read.csv(., row.names = 1))

sc_object_small <- system.file("small_test_data", "sc_object_small.csv",
                               package = "omnideconv", mustWork = TRUE) %>%
  as.matrix(utils::read.csv(., row.names = 1))

cell_annotations_small <- system.file("small_test_data", "cell_annotations_small.txt",
                                      package = "omnideconv", mustWork = TRUE) %>%
  readr::read_lines(.)

batch_ids_small <- system.file("small_test_data", "batch_ids_small.txt",
                                      package = "omnideconv", mustWork = TRUE) %>%
  readr::read_lines(.)

marker_genes <- system.file("small_test_data", "marker_genes_small.txt",
                               package = "omnideconv", mustWork = TRUE) %>%

  readr::read_lines(.)

# sc_object_small <- as.matrix(utils::read.csv("small_test_data", "sc_object_small.csv", row.names = 1))
# cell_annotations_small <- readr::read_lines("small_test_data", "cell_annotations_small.txt")
# batch_ids_small <- readr::read_lines("small_test_data", "batch_ids_small.txt")
# marker_genes <- readr::read_lines("small_test_data", "marker_genes_small.txt")

markers_small <- list(marker_genes[1:9], marker_genes[10:14], marker_genes[15:20])
names(markers_small) <- sort(unique(cell_annotations_small))


test_that("Bisque GenerateSCReference works", {
  signature <- .bisque_patched_model(sc_object_small, cell_annotations_small, batch_ids_small)
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single
               cell matrix", object = ncol(signature),
    expected = length(unique(cell_annotations_small))
  )
  check_signature <- as.matrix(read.csv("test_models", "bisque_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)

  model <- build_model(sc_object_small, cell_annotations_small, "bisque",
    bulk_gene_expression = bulk_small
  )
  expect_null(info = "The Bisque Model is null (which it should be)", object = model)
})


test_that("MOMF compute reference works", {
  signature <- build_model(sc_object_small, cell_annotations_small, "momf",
    bulk_gene_expression = bulk_small
  )
  expect_equal(
    info = "signature matrix has same amount of rows as single cell matrix",
    object = nrow(signature), expected = nrow(sc_object_small)
  )
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
  check_signature <- system.file("test_models", "momf_model_small.csv",
    package = "omnideconv", mustWork = TRUE
  ) %>%
    read.csv(.,
      row.names = 1,
      check.names = FALSE
    ) %>%
  as.matrix(.)
  expect_equal(
    info = "signature matrix is correct", object = signature,
    expected = check_signature
  )
})

test_that("DWLS build signature matrix works", {
  signature <- build_model(sc_object_small, cell_annotations_small, "dwls")
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
  check_signature <- system.file("test_models", "dwls_model_small.csv",
    package = "omnideconv", mustWork = TRUE
  ) %>%
    read.csv(.,
             row.names = 1,
             check.names = FALSE
    ) %>%
    as.matrix(.)
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})

test_that("DWLS build signature matrix works with the optimized version", {
  signature <- build_model(sc_object_small, cell_annotations_small, "dwls", "mast_optimized")
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )

  check_signature <- system.file("test_models", "dwls_model_small.csv",
    package = "omnideconv", mustWork = TRUE
  ) %>%
    read.csv(.,
             row.names = 1,
             check.names = FALSE
    ) %>%
    as.matrix(.)
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})

test_that("CIBERSORTx build signature matrix works", {
  set_cibersortx_credentials(Sys.getenv("CIBERSORTX_EMAIL"), Sys.getenv("CIBERSORTX_TOKEN"))
  signature <- build_model(sc_object_small, cell_annotations_small, "cibersortx")
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single
               cell matrix", object = ncol(signature),
    expected = length(unique(cell_annotations_small))
  )
  check_signature <- system.file("test_models", "cibersortx_model_small.csv",
    package = "omnideconv", mustWork = TRUE
  ) %>%
    read.csv(.,
             row.names = 1,
             check.names = FALSE,
             sep='\t'
    ) %>%
    as.matrix(.)
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
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

test_that("AutoGeneS build model works, and the matrix can be extracted", {
  model <- build_model(sc_object_small, cell_annotations_small, "autogenes")
  expect_true(file.exists(model), "pickle file was created successfully")

  signature <- extract_signature_autogenes(model, sc_object_small, cell_annotations_small)
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single
               cell matrix", object = ncol(signature),
    expected = length(unique(cell_annotations_small))
  )
})

test_that("MuSiC build model works", {
  signature <- build_model(sc_object_small, cell_annotations_small, "music",
    batch_ids = batch_ids_small
  )
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
})
test_that("SCDC build model works", {
  signature <- build_model(sc_object_small, cell_annotations_small, "scdc",
    batch_ids = batch_ids_small
  )
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
})

test_that("CPM build model works", {
  model <- build_model(sc_object_small, cell_annotations_small, "cpm",
    bulk_gene_expression = bulk_small
  )
  expect_null(info = "The CPM Model is null (which it should be)", object = model)
})


test_that("BSeq-sc build model works", {
  signature <- build_model(sc_object_small, cell_annotations_small, "bseqsc",
    batch_ids = batch_ids_small, markers = markers_small
  )
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
  check_signature <- system.file("test_models", "bseqsc_model_small.csv",
    package = "omnideconv", mustWork = TRUE
  ) %>%
    read.csv(.,
             row.names = 1,
             check.names = FALSE
    ) %>%
    as.matrix(.)
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})

test_that("CDSeq build model works", {
  model <- build_model(sc_object_small, cell_annotations_small,
    method = "cdseq"
  )
  expect_null(info = "The CDSeq Model is null (which it should be)", object = model)
})

test_that("BayesPrism build model works", {
  model <- build_model(sc_object_small, cell_annotations_small, "bayesprism",
    bulk_gene_expression = bulk_small
  )
  expect_null(info = "The BayesPrism Model is null (which it should be)", object = model)
})
