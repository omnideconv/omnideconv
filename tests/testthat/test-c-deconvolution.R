bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv", row.names = 1))
bulk_small_one_sample <- bulk_small[, 1, drop = FALSE]
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- readr::read_lines("small_test_data/cell_annotations_small.txt")
batch_ids_small <- readr::read_lines("small_test_data/batch_ids_small.txt")


test_that("Bisque deconvolution works", {
  bisque_model <- as.matrix(read.csv("test_models/bisque_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  deconvolution <- deconvolute(bulk_small, bisque_model,
    method = "bisque", sc_object_small, cell_annotations_small, batch_ids_small
  )
  expect_equal(
    info = "rows of deconvolution equal to columns of signature (same celltypes in same order)",
    object = sort(colnames(deconvolution)), expected = sort(colnames(bisque_model))
  )
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )

  check_result <- as.matrix(read.csv("test_results/bisque_result_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = deconvolution,
    expected = check_result
  )
  expect_error(
    info = "bisque is not appliable with just one bulk sample",
    object = deconvolute(bulk_small_one_sample, bisque_model,
      method = "bisque", sc_object_small,
      cell_annotations_small, batch_ids_small
    )
  )
})


test_that("MOMF deconvolution works", {
  momf_model <- as.matrix(read.csv("test_models/momf_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  deconvolution <- deconvolute(bulk_small, momf_model, method = "momf", sc_object_small)
  expect_equal(
    info = "columns of deconvolution equal to columns of signature (same celltypes in same order)",
    object = sort(colnames(deconvolution)), expected = sort(colnames(momf_model))
  )
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )
  check_result <- as.matrix(read.csv("test_results/momf_result_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = deconvolution,
    expected = check_result
  )
  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, momf_model,
      method = "momf", sc_object_small
    )),
    expected = 1
  )
})

test_that("DWLS deconvolution works", {
  dwls_model <- as.matrix(read.csv("test_models/dwls_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  deconvolution_dwls <- deconvolute(bulk_small, dwls_model,
    method = "dwls", sc_object_small,
    dwls_submethod = "DampenedWLS"
  )
  deconvolution_ols <- deconvolute(bulk_small, dwls_model,
    method = "dwls", sc_object_small,
    dwls_submethod = "OLS"
  )
  deconvolution_svr <- deconvolute(bulk_small, dwls_model,
    method = "dwls", sc_object_small,
    dwls_submethod = "SVR"
  )
  expect_equal(
    info = "rows of deconv for dwls equal to columns of signature (same celltypes, not same order)",
    object = sort(colnames(deconvolution_dwls)), expected = sort(colnames(dwls_model))
  )
  expect_equal(
    info = "deconvolution of dwls contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution_dwls)), expected = sort(colnames(bulk_small))
  )
  expect_equal(
    info = "rows of deconv for ols equal to columns of signature (same celltypes, not same order)",
    object = sort(colnames(deconvolution_ols)), expected = sort(colnames(dwls_model))
  )
  expect_equal(
    info = "deconvolution of ols contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution_ols)), expected = sort(colnames(bulk_small))
  )
  expect_equal(
    info = "rows of deconv for svr equal to columns of signature (same celltypes, not same order)",
    object = sort(colnames(deconvolution_svr)), expected = sort(colnames(dwls_model))
  )
  expect_equal(
    info = "deconvolution of svr contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution_svr)), expected = sort(colnames(bulk_small))
  )

  check_result_dwls <- as.matrix(read.csv("test_results/dwls_dwls_result_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  check_result_ols <- as.matrix(read.csv("test_results/dwls_ols_result_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  check_result_svr <- as.matrix(read.csv("test_results/dwls_svr_result_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result for dwls is correct", object = deconvolution_dwls,
    expected = check_result_dwls
  )
  expect_equal(
    info = "deconvolution result for dwls is correct", object = deconvolution_ols,
    expected = check_result_ols
  )
  expect_equal(
    info = "deconvolution result for dwls is correct", object = deconvolution_svr,
    expected = check_result_svr
  )

  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, dwls_model,
      method = "dwls", sc_object_small,
      dwls_submethod = "DampenedWLS"
    )),
    expected = 1
  )

  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, dwls_model,
      method = "dwls", sc_object_small,
      dwls_submethod = "OLS"
    )),
    expected = 1
  )

  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, dwls_model,
      method = "dwls", sc_object_small,
      dwls_submethod = "SVR"
    )),
    expected = 1
  )
})


test_that("CIBERSORTx deconvolution works", {
  set_cibersortx_credentials(Sys.getenv("CIBERSORTX_EMAIL"), Sys.getenv("CIBERSORTX_TOKEN"))
  cibersort_model <- as.matrix(read.csv("test_models/cibersortx_model_small.tsv",
    row.names = 1,
    check.names = FALSE, sep = "\t"
  ))
  deconvolution <- deconvolute(bulk_small, cibersort_model, method = "cibersortx")

  expect_equal(
    info = "columns of deconvolution equal to columns of signature (same celltypes in same order)",
    object = sort(colnames(deconvolution)), expected = sort(colnames(cibersort_model))
  )
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )

  check_result <- as.matrix(read.csv("test_results/cibersortx_result_small.tsv",
    row.names = 1,
    check.names = FALSE, sep = "\t"
  ))
  check_result <- check_result[, unique(cell_annotations_small)]
  check_result <- check_result[, sort(colnames(check_result))]
  expect_equal(
    info = "deconvolution result is correct", object = deconvolution,
    expected = check_result
  )

  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, cibersort_model, method = "cibersortx")),
    expected = 1
  )
})

test_that("Scaden deconvolution works", {
  model_dir <- paste0(tempdir(), "/model")
  skip_if_not(dir.exists(model_dir), message = "skipping scaden deconvolution test")
  deconvolution <- deconvolute(bulk_small, model_dir, method = "scaden")
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )
  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, model_dir, method = "scaden")),
    expected = 1
  )
})

test_that("Autogenes deconvolution works", {
  files <- file.info(list.files(tempdir(), full.names = T, pattern = "\\.pickle$"))
  skip_if(nrow(files) == 0, message = "skipping autogenes deconvolution")
  model <- rownames(files)[which.max(files$mtime)]
  deconvolution <- deconvolute(bulk_small, model, method = "autogenes")
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )

  check_result <- as.matrix(read.csv("test_results/autogenes_result_small.csv",
    row.names = 1,
    check.names = FALSE,
  ))
  expect_equal(
    info = "deconvolution result is correct", object = deconvolution,
    expected = check_result, tolerance = 1e-5
  )
  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, model, method = "autogenes")),
    expected = 1
  )
})

test_that("MuSiC deconvolution works", {
  deconvolution <- deconvolute(bulk_small, NULL,
    method = "music",
    single_cell_object = sc_object_small,
    cell_type_annotations = cell_annotations_small,
    batch_ids = batch_ids_small
  )
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )

  check_result <- as.matrix(read.csv("test_results/music_result_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = deconvolution,
    expected = check_result
  )
  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, NULL,
      method = "music",
      single_cell_object = sc_object_small,
      cell_type_annotations = cell_annotations_small,
      batch_ids = batch_ids_small
    )),
    expected = 1
  )
})


# This warning can be ignored (according to https://github.com/satijalab/seurat/issues/1249):
# Warning in irlba(A = t(x = object), nv = npcs, ...) :
# You're computing too large a percentage of total singular values, use a standard svd instead.
test_that("CPM deconvolution works", {
  deconvolution_pca <- deconvolute(bulk_small, NULL,
    method = "cpm",
    single_cell_object = sc_object_small,
    cell_type_annotations = cell_annotations_small,
    cell_space = "PCA"
  )

  deconvolution_umap <- deconvolute(bulk_small, NULL,
    method = "cpm",
    single_cell_object = sc_object_small,
    cell_type_annotations = cell_annotations_small,
    cell_space = "UMAP"
  )

  deconvolution_tsne <- deconvolute(bulk_small, NULL,
    method = "cpm",
    single_cell_object = sc_object_small,
    cell_type_annotations = cell_annotations_small,
    cell_space = "TSNE"
  )
  expect_equal(
    info = "deconvolution_pca contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution_pca)), expected = sort(colnames(bulk_small))
  )

  expect_equal(
    info = "deconvolution_umap contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution_umap)), expected = sort(colnames(bulk_small))
  )

  expect_equal(
    info = "deconvolution_tsne contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution_tsne)), expected = sort(colnames(bulk_small))
  )
  expect_error(
    info = "CPM is not appliable with just one bulk sample",
    object = deconvolute(bulk_small_one_sample, NULL,
      method = "cpm",
      single_cell_object = sc_object_small,
      cell_type_annotations = cell_annotations_small,
      cell_space = "PCA"
    )
  )
})


test_that("SCDC deconvolution works", {
  deconvolution <- deconvolute(bulk_small, NULL,
    method = "scdc", single_cell_object = sc_object_small,
    cell_type_annotations = cell_annotations_small,
    batch_ids = batch_ids_small
  )
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )

  check_result <- as.matrix(read.csv("test_results/scdc_result_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = deconvolution,
    expected = check_result
  )
  expect_error(
    info = "scdc is not appliable with just one bulk sample",
    object = deconvolute(bulk_small_one_sample, NULL,
      method = "scdc", single_cell_object = sc_object_small,
      cell_type_annotations = cell_annotations_small,
      batch_ids = batch_ids_small
    )
  )


  # deconvolution_ensemble <- deconvolute(bulk_small, NULL,
  #  method = "scdc",
  #  single_cell_object = list(
  #    single_cell_data,
  #    sc_object_small
  #  ),
  #  cell_type_annotations = list(
  #    cell_type_annotations,
  #    cell_annotations_small
  #  ),
  #  batch_ids = list(batch_ids, batch_ids_small)
  # )
  # expect_equal(
  #  info = "deconvolution contains same samples as in bulk (not same order)",
  #  object = sort(rownames(deconvolution_ensemble)), expected = sort(colnames(bulk_small))
  # )
  #
  # check_result <- as.matrix(read.csv("test_results/scdc_result_ensemble.csv",
  #  row.names = 1,
  #  check.names = FALSE
  # ))
  # expect_equal(
  #  info = "deconvolution result is correct", object = deconvolution_ensemble,
  #  expected = check_result
  # )
})

test_that("CDSeq deconvolution works", {
  deconvolution <- deconvolute(bulk_small, NULL,
    method = "cdseq",
    single_cell_object = sc_object_small,
    cell_type_annotations = cell_annotations_small,
    batch_ids = batch_ids_small
  )
  expect_equal(
    info = "deconvolution contains same samples as in bulk (not same order)",
    object = sort(rownames(deconvolution)), expected = sort(colnames(bulk_small))
  )
  expect_equal(
    info = "deconvolution result with one bulk sample throws no error",
    object = nrow(deconvolute(bulk_small_one_sample, NULL,
      method = "cdseq",
      single_cell_object = sc_object_small,
      cell_type_annotations = cell_annotations_small,
      batch_ids = batch_ids_small
    )),
    expected = 1
  )
})
