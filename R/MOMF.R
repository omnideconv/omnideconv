#' Calculates the signature model with MOMF
#'
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples. Row and
#'   column names need to be set.
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
build_model_momf <- function(single_cell_object, cell_type_annotations, bulk_gene_expression) {
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  MOMF::momf.computeRef(
    single_cell_object[intersect(
      rownames(single_cell_object),
      rownames(bulk_gene_expression)
    ), ],
    cell_type_annotations
  )
}

#' Deconvolution Analysis using MOMF (via Nonnegative Factorization)
#'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param signature The signature matrix. Rows are genes, columns are cell types.
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param method Determines which divergence to use. Options: Kullback-Leibler "KL",
#'   Itakura-Saito "IS". Defaults to "KL"
#' @param verbose Whether to produce an output on the console.
#' @param ... additional parameters
#'
#' @return cell proportion matrix
#' @export
#'
deconvolute_momf <- function(bulk_gene_expression, signature, single_cell_object,
                             verbose = FALSE, method = "KL", ...) {
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  if (is.null(signature)) {
    stop("Parameter 'signature' is missing or null, but it is required.")
  }
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  # MOMF needs a list of the single_cell_object with cells x genes and the bulk RNA seq data with
  # individuals x genes
  # IMPORTANT: This line with the intersection is a difference from the original algorithm
  relevant_genes <-
    intersect(
      intersect(rownames(single_cell_object), rownames(bulk_gene_expression)),
      rownames(signature)
    )
  GList <- list(
    X1 = t(single_cell_object[relevant_genes, , drop = FALSE]),
    X2 = t(bulk_gene_expression[relevant_genes, , drop = FALSE])
  )
  signature <- signature[relevant_genes, , drop = FALSE]
  if (!verbose) {
    sink(tempfile())
    result <- tryCatch(MOMF::momf.fit(DataX = GList, DataPriorU = signature, method = method, ...),
      finally = sink()
    )
  } else {
    result <- MOMF::momf.fit(DataX = GList, DataPriorU = signature, method = method, ...)
  }
  if (is.null(result$cell.prop)) {
    stop("Something went wrong. Please switch on verbose mode")
  }
  # return slot in result with cell proportion matrix
  return(result)
}
