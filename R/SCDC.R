# Requires raw read counts.


#' No model is build as SCDC does both steps in one.
#'
#' Please use the deconvolute method with your single cell and bulk rna seq data to use SCDC.
#'
#'
#' @return NULL
#'
#' @export
build_model_scdc <- function() {
  base::message(
    "The deconvolution with SCDC is done in only one step. Please just use the",
    "deconvolute method."
  )

  return(NULL)
}

#' SCDC Deconvolution
#'
#' This function is to calculate the SCDC deconvolution proportions.
#' IMPORTANT: No model is needed. Everything is done inside this method.
#'
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset.list list of ExpressionSet objects for single cell reference datasets. Note that
#'   the variable names of multiple ExpressionSet objects should be consistent. Not required if
#'   prop.input is specified.
#' @param ct_sub vector. a subset of cell types that are selected to construct basis matrix.
#' @param ct_varname character string specifying the variable name for 'cell types'.
#' @param sample character string specifying the variable name for subject/samples.
#' @param iter_max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @param ct_cell_size default is NULL, which means the "library size" is calculated based on the
#'   data. Users can specify a vector of cell size factors corresponding to the ct.sub according to
#'   prior knowledge. The vector should be named: names(ct_cell_size input) should not be NULL.
#' @param Transform_bisque The bulk sample transformation from bisqueRNA. Aiming to reduce the
#'   systematic difference between single cells and bulk samples.
#' @param grid_search logical. whether to allow grid search method to derive the ENSEMBLE weights.
#' @param search_length a number between 0 to 0.5. if using "Grid search", the step length used.
#'   Smaller search.length derives more accurate optimization results.
#' @param verbose Whether to create any output.
#'
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
deconvolute_scdc <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                             ct_varname = "cellType", sample = "SubjectName", ct_sub = NULL,
                             iter_max = 1000, nu = 1e-04, epsilon = 0.01, truep = NULL,
                             weight_basis = T, ct_cell_size = NULL, Transform_bisque = F,
                             grid_search = F, search_length = 0.05, names_sc_objects = NULL,
                             verbose = FALSE) {
  if (is.null(single_cell_object) || is.null(cell_type_annotations)) {
    base::stop(
      "Single cell object or cell type annotations not provided. Call as: ",
      "deconvolute(bulk_gene_expression, NULL, \"scdc\", single_cell_object, ",
      "cell_type_annotations)"
    )
  }
  sc_eset <- mapply(function(sc_obj, cell_anno) {
    return(get_single_cell_expression_set(sc_obj, colnames(sc_obj), rownames(sc_obj), cell_anno))
  }, single_cell_object, cell_type_annotations)
  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)

  if (lenght(sc_eset) == 1) {
    sc_eset <- unlist(sc_eset)

    return(SCDC::SCDC_prop(bulk_eset, sc_eset,
      ct.varname = ct_varname, sample = sample,
      ct.sub = ct_sub, iter.max = iter_max, nu = nu, epsilon = epsilon,
      truep = truep, weight.basis = weight_basis, ct.cell.size = ct_cell_size,
      Transform_bisque = Transform_bisque
    ))
  } else {
    if (!is.null(names_sc_objects)) {
      if (length(sc_eset) != length(names_sc_objects)) {
        base::stop(
          "The", length(names_sc_objects), "names supplied are not the same number of",
          "names as the number of single cell objects (", length(sc_eset), ")"
        )
      }
      names(sc_eset) <- names_sc_objects
    } else {
      names(sc_eset) <- seq_len(length(sc_eset))
    }
    return(SCDC::SCDC_ENSEMBLE(bulk_eset, sc_eset,
      ct.varname = ct_varname, sample = sample,
      ct.sub = ct_sub, iter.max = iter_max, nu = nu, epsilon = epsilon,
      truep = truep, weight.basis = weight_basis,
      ct.cell.size = ct_cell_size, Transform_bisque = Transform_bisque,
      grid.search = grid_search, search.length = search_length
    ))
  }

  # Check number of indiviuals -> ONE method?
  # if qc do that
  # check if input single cell is list or single expression set


  return(MuSiC::music_prop(bulk_eset, sc_eset,
    markers = markers, clusters = clusters,
    samples = samples, select.ct = select_ct, cell_size = cell_size,
    ct.cov = ct_cov, verbose = verbose, iter.max = iter_max, nu = nu,
    eps = eps, centered = centered, normalize = normalize
  ))
}
