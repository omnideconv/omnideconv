bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv", row.names = 1))
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- utils::read.csv("small_test_data/cell_annotations_small.csv",
  row.names = 1
)$x

test_that("Bisque GenerateSCReference works", {
  signature <- build_model(sc_object_small, cell_annotations_small, method = "bisque")
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single
               cell matrix", object = ncol(signature),
    expected = length(unique(cell_annotations_small))
  )
  check_signature <- as.matrix(read.csv("test_models/bisque_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})


test_that("MOMF compute reference works", {
  signature <- build_model(sc_object_small, cell_annotations_small,
    method = "momf",
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
  check_signature <- as.matrix(read.csv("test_models/momf_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(
    info = "signature matrix is correct", object = signature,
    expected = check_signature
  )
})

test_that("DWLS build signature matrix works", {
  signature <- build_model(sc_object_small, cell_annotations_small, method = "dwls")
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
  check_signature <- as.matrix(read.csv("test_models/dwls_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})

test_that("CibersortX build signature matrix works", {
  set_cibersortx_credentials("konstantin.pelz@tum.de", "27308ae0ef1458d381becac46ca7e480")
  signature <- build_model(sc_object_small, cell_annotations_small, method = "cibersortx")
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single
               cell matrix", object = ncol(signature),
    expected = length(unique(cell_annotations_small))
  )
  check_signature <- as.matrix(read.csv("test_models/cibersortx_model_small.tsv",
    row.names = 1,
    check.names = FALSE, sep = "\t"
  ))
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

test_that("AutoGeneS build model works", {
  model <- build_model(sc_object_small, cell_annotations_small, method = "autogenes")
  expect_true(file.exists(model), "pickle file was created successfully")
})

test_that("MuSiC build model works", {
  model <- build_model(sc_object_small, cell_annotations_small, method = "music",
                       bulk_gene_expression = bulk_small)
  expect_equal(info="The MuSiC Model contains five elements", object=length(model),expected=5)
})
