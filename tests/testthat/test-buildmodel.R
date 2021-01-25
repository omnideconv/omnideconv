test_that("Bisque GenerateSCReference works",{
  gene_expression = diag(5)
  gene_names = c("marker_A", "marker_B", "marker_C", "marker_D", "marker_E")
  cell_types = c("A", "B", "C", "D", "E")
  pdata_df = data.frame(cell_type=cell_types)
  rownames(pdata_df) = cell_types
  colnames(gene_expression) = cell_types
  rownames(gene_expression) = gene_names
  expect_equal(object = nrow(build_model(gene_expression, cell_types, method = "bisque")), expected = nrow(gene_expression))
  expect_equal(object = ncol(build_model(gene_expression, cell_types, method = "bisque")), expected = length(unique(cell_types)))
})


test_that("MOMF compute reference works",{
  gene_expression = diag(5)
  gene_names = c("marker_A", "marker_B", "marker_C", "marker_D", "marker_E")
  cell_types = c("A", "B", "C", "D", "E")
  pdata_df = data.frame(cell_type=cell_types)
  rownames(pdata_df) = cell_types
  colnames(gene_expression) = cell_types
  rownames(gene_expression) = gene_names
  expect_equal(object = nrow(build_model(gene_expression, cell_types, method = "momf")), expected = nrow(gene_expression))
  expect_equal(object = ncol(build_model(gene_expression, cell_types, method = "momf")), expected = length(unique(cell_types)))
})

# test_that("DWLS build signature matrix works", {
#   gene_expression = diag(5)
#   gene_names = c("marker_A", "marker_B", "marker_C", "marker_D", "marker_E")
#   cell_types = c("A", "B", "C", "D", "E")
#   pdata_df = data.frame(cell_type=cell_types)
#   rownames(pdata_df) = cell_types
#   colnames(gene_expression) = cell_types
#   rownames(gene_expression) = gene_names
#   expect_equal(object = ncol(build_model(gene_expression, cell_types, method = "dwls")), expected = length(unique(cell_types)))
# })

test_that("Scaden build model works", {
  gene_expression = diag(5)
  gene_names = c("marker_A", "marker_B", "marker_C", "marker_D", "marker_E")
  cell_types = c("A", "B", "C", "D", "E")
  pdata_df = data.frame(cell_type=cell_types)
  rownames(pdata_df) = cell_types
  colnames(gene_expression) = cell_types
  rownames(gene_expression) = gene_names
})
