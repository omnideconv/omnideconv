#' No model is build as CDSeq does both steps in one.
#'
#' Please use the deconvolute method with your single cell and bulk rna seq data to use CDSeq.
#'
#'
#' @return NULL
#'
#' @export
build_model_cdseq <- function() {
  base::message(
    "The deconvolution with CDSeq is done in only one step. Please just use the",
    "deconvolute method."
  )

  return(NULL)
}

#' CDSeq Deconvolution
#'
#' This function is to calculate the CDSeq deconvolution proportions.
#' IMPORTANT: No model is needed. Everything is done inside this method.
#'
#' @param bulk_gene_expression  A matrix or dataframe with the bulk data. Rows
#'   are genes, columns are samples.
#' @param single_cell_object A Matrix with the single-cell data. Rows are genes
#'  and columns are samples.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to
#'  be in the same order as the samples in single_cell_object.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @return NULL
#' @export
deconvolute_cdseq <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                              batch_ids) {
  if (is.null(single_cell_object) || is.null(cell_type_annotations) || is.null(batch_ids)) {
    base::stop(
      "Single cell object or cell type annotations not provided. Call as: ",
      "deconvolute(bulk_gene_expression, NULL, \"music\", single_cell_object, ",
      "cell_type_annotations, batch_ids)"
    )
  }
  return(NULL)
}
