sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- utils::read.csv("small_test_data/cell_annotations_small.csv",row.names = 1)$x

test_that("Matrix/SingleCellExperiment conversion works",{


  sce <- matrix_to_singlecellexperiment(sc_object_small, cell_annotations_small)
  expect_identical(info = "Conversion from Matrix to SingleCellExperiment had no error", object = typeof(sce), expected = "S4")
  X <- as.list(SummarizedExperiment::assays(sce))[[1]]
  labels <- SingleCellExperiment::colData(sce)$label
  expect_equal(info = "SCE Matrix is correct", object = X, expected = sc_object_small)
  expect_equal(info = "SCE annotation vector is correct", object = labels, expected = cell_annotations_small)

  matrix_and_annotation <- singlecellexperiment_to_matrix(sce)
  expect_identical(info = "Conversion from SingleCellExperiment to Matrix had no error", object = typeof(matrix_and_annotation), expected = "list")
  matrix <- matrix_and_annotation$matrix
  annotation <- matrix_and_annotation$annotation_vector
  expect_equal(info = "Matrix is correct", object = matrix, expected = sc_object_small)
  expect_equal(info = "Annotation vector is correct", object = annotation, expected = cell_annotations_small)

})



test_that("SingleCellExperiment/Anndata conversion works",{
  ad <- anndata::AnnData(
    X = matrix(1:6, nrow = 2),
    obs = data.frame(group = c("a", "b"), row.names = c("s1", "s2")),
    var = data.frame(type = c(1L, 2L, 3L), row.names = c("var1", "var2", "var3")),
    layers = list(
      spliced = matrix(4:9, nrow = 2),
      unspliced = matrix(8:13, nrow = 2)
    ),
    obsm = list(
      ones = matrix(rep(1L, 10), nrow = 2),
      rand = matrix(rnorm(6), nrow = 2),
      zeros = matrix(rep(0L, 10), nrow = 2)
    ),
    varm = list(
      ones = matrix(rep(1L, 12), nrow = 3),
      rand = matrix(rnorm(6), nrow = 3),
      zeros = matrix(rep(0L, 12), nrow = 3)
    ),
    uns = list(
      a = 1,
      b = data.frame(i = 1:3, j = 4:6, value = runif(3)),
      c = list(c.a = 3, c.b = 4)
    )
  )

  sce_converted <- anndata_to_singlecellexperiment(ad)

  ad_converted <- singlecellexperiment_to_anndata(sce_converted)


  expect_equal(info = "Anndata conversion to SCE did not produce an error", object = typeof(sce_converted), expected = "S4")
  expect_equal(info = "SCE conversion to Anndata did not produce an error", object = typeof(ad), expected = "environment")

})

test_that("SingleCellExperiment/Anndata conversion does not lose information",{
  ad <- anndata::AnnData(
    X = matrix(1:6, nrow = 2),
    obs = data.frame(group = c("a", "b"), row.names = c("s1", "s2")),
    var = data.frame(type = c(1L, 2L, 3L), row.names = c("var1", "var2", "var3")),
    layers = list(
      spliced = matrix(4:9, nrow = 2),
      unspliced = matrix(8:13, nrow = 2)
    )
  )

  sce <- SingleCellExperiment::SingleCellExperiment(list(X=t(matrix(1:6, nrow = 2)), spliced=t(matrix(4:9, nrow = 2)), unspliced=t(matrix(8:13, nrow = 2))),
                              colData=data.frame(group = c("a", "b"), row.names = c("s1", "s2")),
                              rowData=data.frame(type = c(1L, 2L, 3L), row.names = c("var1", "var2", "var3"))
  )

  sce_converted <- anndata_to_singlecellexperiment(ad)

  ad_converted <- singlecellexperiment_to_anndata(sce)


  expect_equal(info = "Anndata conversion to SCE did not produce an error", object = typeof(sce_converted), expected = "S4")
  expect_equal(info = "SCE conversion to Anndata did not produce an error", object = typeof(ad), expected = "environment")
  expect_equal(info = "Conversion from Anndata to SCE is correct", object = sces_are_identical(sce, sce_converted), expected = TRUE)
  expect_equal(info = "Conversion from SCE to Anndata is correct", object = anndata_is_identical(ad,ad_converted), expected = TRUE)
})


