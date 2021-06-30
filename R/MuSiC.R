#' No model is build as MuSiC does both steps in one.
#'
#' Please use the deconvolute method with your single cell and bulk rna seq data to use MuSiC.
#'
#'
#' @return NULL
#'
#' @export
build_model_music <- function() {
  base::message(
    "The deconvolution with MuSiC is done in only one step. Please just use the",
    "deconvolute method."
  )

  return(NULL)
}

#' MuSiC Deconvolution
#'
#' This function is to calculate the MuSiC deconvolution proportions.
#' IMPORTANT: No model is needed. Everything is done inside this method.
#'
#' @param bulk_gene_expression  A matrix or dataframe with the bulk data. Rows
#'   are genes, columns are samples.
#' @param single_cell_object A Matrix with the single-cell data. Rows are genes
#'  and columns are samples.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to
#'  be in the same order as the samples in single_cell_object.
#' @param markers vector or list of gene names, default as NULL. If NULL, use all genes that
#'  provided by both bulk and single cell dataset.
#' @param clusters character, the phenoData of single cell dataset used as clusters.
#' @param samples character,the phenoData of single cell dataset used as samples.
#' @param select_ct vector of cell types. Default as NULL. If NULL, then use all cell types
#'  provided.
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd
#'  column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from
#'  data.
#' @param ct_cov logical. If TRUE, use the covariance across cell types.
#' @param verbose logical, default as TRUE.
#' @param iter_max numeric, maximum iteration number
#' @param nu regulation parameter, take care of weight when taking recipical.
#' @param eps Thredshold of convergence.
#' @param centered logic, substract avg of Y and D.
#' @param normalize logic, divide Y and D by their standard deviation.
#' @return a list with elements:
#' \itemize{
#'   \item Estimates of MuSiC
#'   \item Estimates of NNLS
#'   \item Weight of MuSiC
#'   \item r.squared of MuSiC
#'   \item Variance of MuSiC estimates
#' }
#' @importFrom Biobase exprs pData
#' @export
deconvolute_music <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                              markers = NULL, clusters = "cellType", samples = "SubjectName",
                              select_ct = NULL, cell_size = NULL, ct_cov = FALSE, verbose = FALSE,
                              iter_max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE,
                              normalize = FALSE) {
  if (is.null(single_cell_object) || is.null(cell_type_annotations)) {
    base::stop(
      "Single cell object or cell type annotations not provided. Call as: ",
      "deconvolute(bulk_gene_expression, NULL, \"music\", single_cell_object, ",
      "cell_type_annotations)"
    )
  }
  sc_eset <- get_single_cell_expression_set(
    single_cell_object, colnames(single_cell_object),
    rownames(single_cell_object), cell_type_annotations
  )
  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)

  return(MuSiC::music_prop(bulk_eset, sc_eset,
    markers = markers, clusters = clusters,
    samples = samples, select.ct = select_ct, cell_size = cell_size,
    ct.cov = ct_cov, verbose = verbose, iter.max = iter_max, nu = nu,
    eps = eps, centered = centered, normalize = normalize
  ))
}
