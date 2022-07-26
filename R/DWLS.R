#' Signature matrix creation with DWLS using genes identified by a differential analysis
#'
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param dwls_method The method used to create the signature matrix. Options are "mast" and "seurat"
#' @param path The path where the generated files will be saved. If path=NULL, the generated files
#'   will be discarded.
#' @param verbose Whether to produce an output on the console.
#' @param diff_cutoff Cutoff to determine the FC-limit. How low can the lowest fold change be to
#'   still be considered differentially expressed?
#' @param pval_cutoff Cutoff to determine the pVal-limit. How high can the highest p-Value be to
#'   still be considered statistically significant?
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
build_model_dwls <- function(single_cell_object, cell_type_annotations,
                             dwls_method = c("mast", "seurat"), path = NULL, verbose = FALSE,
                             diff_cutoff = 0.5, pval_cutoff = 0.01) {
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (length(dwls_method) > 1) {
    dwls_method <- dwls_method[1]
    message(paste0(dwls_method, " was chosen because multiple values were supplied for \"dwls_method\""))
  }

  if (dwls_method == "mast") {
    return(DWLS::buildSignatureMatrixMAST(
      single_cell_object, cell_type_annotations, path,
      verbose, diff_cutoff, pval_cutoff
    ))
  } else if (dwls_method == "seurat") {
    return(DWLS::buildSignatureMatrixUsingSeurat(
      single_cell_object, cell_type_annotations, path,
      verbose, diff_cutoff, pval_cutoff
    ))
  } else {
    stop("Could not find dwls_method " + dwls_method + ". Please try \"mast\" or \"seurat\"")
  }
}
#' Calculates the decomposition using the dwls algorithm
#'
#' generates a reference profile based on single-cell data. Learns a transformation of bulk
#' expression based on observed single-cell proportions and performs NNLS regression on these
#' transformed values to estimate cell proportions.
#'
#' @param bulk_gene_expression An Expression Set containing bulk data.
#' @param signature The Signature matrix.
#' @param dwls_submethod Three alternative methods in DWLS: OLS, SVR, and DampenedWLS.
#' @param verbose Whether to produce an output on the console.

#'
#' @return A matrix of cell type proportion estimates with cell types as rows and individuals as
#'   columns.
#' @export
#'

deconvolute_dwls <- function(bulk_gene_expression, signature,
                             dwls_submethod = c("DampenedWLS", "OLS", "SVR"), verbose = FALSE) {
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  if (is.null(signature)) {
    stop("Parameter 'signature' is missing or null, but it is required.")
  }
  if (length(dwls_submethod) > 1) {
    dwls_submethod <- dwls_submethod[1]
    message(paste0(
      dwls_submethod, " was chosen because multiple values were supplied for ",
      "\"dwls_submethod\""
    ))
  }

  if (verbose) {
    message("Running DWLS deconvolution module")
  }

  # trim data
  genes <- intersect(rownames(signature), rownames(bulk_gene_expression))
  bulk <- bulk_gene_expression[genes, , drop = FALSE]
  sig <- signature[genes, , drop = FALSE]
  if (class(bulk)[[1]] == "numeric" || class(sig)[[1]] == "numeric") {
    stop("Either bulk data or signature matrix just contains one row!")
  }

  # perform reconvolution in different sub_methods
  res <- NULL

  if (dwls_submethod == "OLS") {
    solutions_ols <- NULL
    for (i in 1:ncol(bulk)) {
      bulk_i <- bulk[, i]
      sol <- DWLS::solveOLS(sig, bulk_i, verbose)
      # sol<-round(sol,5)
      solutions_ols <- cbind(solutions_ols, sol)
    }
    colnames(solutions_ols) <- colnames(bulk)
    res <- solutions_ols
  } else if (dwls_submethod == "SVR") {
    solutions_svr <- NULL
    for (i in 1:ncol(bulk)) {
      bulk_i <- bulk[, i]
      sol <- DWLS::solveSVR(sig, bulk_i, verbose)
      # sol<-round(sol,5)
      solutions_svr <- cbind(solutions_svr, sol)
    }
    colnames(solutions_svr) <- colnames(bulk)
    res <- solutions_svr
  } else if (dwls_submethod == "DampenedWLS") {
    solutions_dampened_wls <- NULL
    for (i in 1:ncol(bulk)) {
      bulk_i <- bulk[, i]
      sol <- DWLS::solveDampenedWLS(sig, bulk_i)
      # sol<-round(sol,5)
      solutions_dampened_wls <- cbind(solutions_dampened_wls, sol)
    }
    colnames(solutions_dampened_wls) <- colnames(bulk)
    res <- solutions_dampened_wls
  } else {
    stop("Submethod " + dwls_submethod + " not found. Please provide a valid one.")
  }

  if (verbose) {
    message("Deconvolution sucessful!")
  }
  result <- t(res)
  return(result)
}
