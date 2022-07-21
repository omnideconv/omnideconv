#' Signature creation with MuSiC.
#'
#' MuSiC does the signature creation in one step, not separated into build_model and deconvolute.
#' Please use the deconvolute method with your single cell and bulk RNA seq data to use MuSiC.
#'
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to
#'   be in the same order as the samples in single_cell_object.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param non_zero logical, default as TRUE. If true, remove all gene with zero expression.
#' @param markers vector or list of gene names, default as NULL. If NULL, use all genes.
#' @param clusters character, the name of the phenoData of the single cell dataset used as clusters.
#' @param samples character, the name of the phenoData of the single cell dataset used as samples.
#' @param select_ct vector of cell types. Default as NULL. If NULL, then use all cell types
#'   provided.
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd
#'   column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size
#'   from data.
#' @param ct_cov logical. If TRUE, use the covariance across cell types.
#' @param verbose Whether to produce an output on the console.
#' @return a list with elements:
#' \itemize{
#'   \item gene by cell type matrix of Design matrix
#'   \item subject by celltype matrix of Library size
#'   \item vector of average library size for each cell type
#'   \item gene by celltype matrix of average relative abundance
#'   \item gene by celltype matrix of cross-subject variation
#' }
#'
#' @importFrom Biobase exprs pData
#' @export
#'
build_model_music <- function(single_cell_object, cell_type_annotations, batch_ids, non_zero = TRUE,
                              markers = NULL, clusters = "cellType", samples = "batchId",
                              select_ct = NULL, cell_size = NULL, ct_cov = FALSE, verbose = FALSE) {
  message(
    "The deconvolution with MuSiC is done in only one step. Please just use the ",
    "deconvolute method. You can still calculate a signature matrix with MuSiC, ",
    "just not input one for the deconvolution step."
  )

  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (is.null(batch_ids)) {
    stop("Parameter 'batch_ids' is missing or null, but it is required.")
  }
  sc_eset <- get_single_cell_expression_set(
    single_cell_object, batch_ids,
    rownames(single_cell_object), cell_type_annotations
  )
  return(MuSiC::music_basis(sc_eset,
    non.zero = non_zero, markers = markers,
    clusters = clusters, samples = samples, select.ct = select_ct,
    cell_size = cell_size, ct.cov = ct_cov, verbose = verbose
  ))
}

#' MuSiC Deconvolution
#'
#' This function is to calculate the MuSiC deconvolution proportions.
#' IMPORTANT: No model is needed. Everything is done inside this method.
#'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to
#'   be in the same order as the samples in single_cell_object.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param markers vector or list of gene names, default as NULL. If NULL, use all genes that
#'   are provided by both bulk and single cell dataset.
#' @param clusters character, the phenoData of single cell dataset used as clusters.
#' @param samples character,the phenoData of single cell dataset used as samples.
#' @param select_ct vector of cell types. Default as NULL. If NULL, then use all cell types
#'   provided.
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd
#'   column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size
#'   from data.
#' @param ct_cov logical. If TRUE, use the covariance across cell types.
#' @param verbose Whether to produce an output on the console.
#' @param iter_max numeric, maximum iteration number.
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
                              batch_ids, markers = NULL, clusters = "cellType", samples = "batchId",
                              select_ct = NULL, cell_size = NULL, ct_cov = FALSE, verbose = FALSE,
                              iter_max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE,
                              normalize = FALSE) {
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  if (ncol(bulk_gene_expression) < 2) {
    stop("MuSiC requires at least two bulk samples.")
  }
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (is.null(batch_ids)) {
    stop("Parameter 'batch_ids' is missing or null, but it is required.")
  }
  sc_eset <- get_single_cell_expression_set(
    single_cell_object, batch_ids,
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
