#' Signature creation with SCDC.
#'
#' SCDC does the signature creation in one step, not separated into build_model and deconvolute.
#' Please use the deconvolute method with your single cell and bulk RNA seq data to use SCDC.
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes,
#'   columns are samples. Row and column names need to be set. This can also be a list of objects,
#'   if SCDC_ENSEMBLE should be used.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object. This can also be a list of vectors, if SCDC_ENSEMBLE
#'   should be used.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param ct_sub vector. a subset of cell types that are selected to construct basis matrix. NULL
#'   means that all are used.
#' @param ct_varname character string specifying the variable name for 'cell types'.
#' @param sample character string specifying the variable name for subject/samples.
#' @param ct_cell_size default is NULL, which means the "library size" is calculated based on the
#'   data. Users can specify a vector of cell size factors corresponding to the ct.sub according to
#'   prior knowledge. The vector should be named: names(ct_cell_size input) should not be NULL.
#' @param verbose Whether to produce an output on the console.
#' @return a list with elements:
#' \itemize{
#'   \item basis matrix
#'   \item sum of cell-type-specific library size
#'   \item sample variance matrix
#'   \item basis matrix by mvw
#'   \item mvw matrix
#' }
#'
#' @export
build_model_scdc <- function(single_cell_object, cell_type_annotations, batch_ids, ct_sub = NULL,
                             ct_varname = "cellType", sample = "batchId", ct_cell_size = NULL,
                             verbose = FALSE) {
  message(
    "The deconvolution with SCDC is done in only one step. Please just use the ",
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
  if (length(unique(sc_eset@phenoData@data[, sample])) > 1) {
    return(SCDC::SCDC_basis(
      x = sc_eset, ct.sub = ct_sub, ct.varname = ct_varname,
      sample = sample, ct.cell.size = ct_cell_size
    ))
  } else {
    return(SCDC::SCDC_basis_ONE(
      x = sc_eset, ct.sub = ct_sub,
      ct.varname = ct_varname, sample = sample, ct.cell.size = ct_cell_size
    ))
  }
}

#' SCDC Deconvolution
#'
#' This function is to calculate the SCDC deconvolution proportions.
#' IMPORTANT: No model is needed. Everything is done inside this method.
#'
#' SCDC_ENSEMBLE can be used by supplying lists to the parameters single_cell_object and
#' cell_type_annotations. To name the single cell data sets, supply a vector with their
#' corresponding names to names_sc_objects.
#'
#' Requires raw read counts.
#' Works best with multiple cells per single cell patient/subject
#'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes,
#'   columns are samples. Row and column names need to be set. This can also be a list of objects,
#'   if SCDC_ENSEMBLE should be used.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object. This can also be a list of vectors, if SCDC_ENSEMBLE
#'   should be used.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param ct_sub vector. a subset of cell types that are selected to construct basis matrix. NULL
#'   means that all are used.
#' @param ct_varname character string specifying the variable name for 'cell types'.
#' @param sample character string specifying the variable name for subject/samples.
#' @param iter_max the maximum number of iteration in WNNLS. If the parameter is NULL, the default
#'   value of 1000 for a single single cell object and 2000 for a list is chosen.
#' @param nu a small constant to facilitate the calculation of variance.
#' @param epsilon a small constant number used for convergence criteria. If the parameter is NULL,
#'   the default value of 0.01 for a single single cell object and 0.001 for a list is chosen.
#' @param truep true cell-type proportions for bulk samples if known.
#' @param weight_basis Whether to use the Basis Matrix adjusted for maximal variance weight,
#'   created by the SCDC_basis function.
#' @param ct_cell_size default is NULL, which means the "library size" is calculated based on the
#'   data. Users can specify a vector of cell size factors corresponding to the ct.sub according to
#'   prior knowledge. The vector should be named: names(ct_cell_size input) should not be NULL.
#' @param transform_bisque The bulk sample transformation from bisqueRNA. Aiming to reduce the
#'   systematic difference between single cells and bulk samples.
#' @param grid_search logical. whether to allow grid search method to derive the ENSEMBLE weights.
#' @param search_length a number between 0 to 0.5. if using "Grid search", the step length used.
#'   Smaller search.length derives more accurate optimization results.
#' @param names_sc_objects A vector with the names of the single cell objects. Only used if
#'   a list of single cell objects is supplied. If it remains NULL, the objects are named by their
#'   index.
#' @param qcthreshold The probability threshold used to filter out questionable cells, only used if
#'   quality_control = TRUE.
#' @param verbose Whether to produce an output on the console.
#' @param quality_control Whether to perform the SCDC_qc quality control method.
#'
#' @return Depends on whether one or multiple single cell sets are used.\cr\cr
#' One:
#' \item{prop.est.mvw}{A matrix of cell type proportion estimates with cell types as rows and
#'   individuals as columns.}
#' \item{basis.mvw}{The signature matrix. Rows are genes, columns are cell types.}
#' \item{yhat}{The predicted gene expression levels for the bulk samples.}
#' \item{yeval}{The evaluation of the predicted gene expression levels.}
#' \item{peval}{The evaluation of the deconvoluted proportions. Since we dont have a ground truth,
#'   this is always NULL.}
#'
#' ENSEMBLE:
#' \item{w_table}{A matrix with the suggested weights for each single cell dataset and some
#'   statistical evaluation.}
#' \item{prop.list}{A list of the "One:" outputs as seen above for each single cell dataset.}
#' \item{prop.only}{A list of the prop.est.mvw values for each single cell dataset.}
#' \item{gridres}{A matrix with the results of the gridsearch. NULL if grid_search = FALSE.}
#'
#' @export
deconvolute_scdc <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                             batch_ids, ct_varname = "cellType", sample = "batchId",
                             ct_sub = NULL, iter_max = NULL, nu = 1e-04, epsilon = NULL,
                             truep = NULL, weight_basis = TRUE, ct_cell_size = NULL,
                             transform_bisque = FALSE, grid_search = FALSE, search_length = 0.05,
                             names_sc_objects = NULL, qcthreshold = 0.7, verbose = FALSE,
                             quality_control = FALSE) {
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
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

  if (ncol(bulk_gene_expression) < 2) {
    stop("SCDC requires at least two bulk samples.")
  }

  if (!"list" %in% class(single_cell_object)) {
    single_cell_object <- list(single_cell_object)
    cell_type_annotations <- list(cell_type_annotations)
    batch_ids <- list(batch_ids)
    if (is.null(ct_sub)) {
      ct_sub <- unique(unlist(cell_type_annotations))
    }
  }
  if (is.null(ct_sub)) {
    ct_sub <- Reduce(intersect, cell_type_annotations)
  }

  sc_eset <- mapply(function(sc_obj, cell_anno, batch_ids_curr) {
    if (!"ExpressionSet" %in% class(sc_obj)) {
      if (length(batch_ids_curr) == length(unique(batch_ids_curr))) {
        message(
          "Each batch_id of at least one list item only contained one cell. If an error regarding ",
          "the number of valid cell types occurs, try with \"weight_basis=FALSE\" or ",
          "with \"quality_control=FALSE\""
        )
      }
      sc_obj <- get_single_cell_expression_set(
        sc_obj, batch_ids_curr, rownames(sc_obj),
        cell_anno
      )
    }
    if (quality_control) {
      # SCDC always supplies two types of each method, one normal one and one if all the single
      # cells are from the same individual (_ONE)
      if (length(unique(sc_obj@phenoData@data[, sample])) > 1) {
        sc_obj <- SCDC::SCDC_qc(sc_obj, ct_varname, sample, ct_sub,
          iter_max = iter_max, nu = nu,
          epsilon = epsilon, qcthreshold = qcthreshold, generate.figure = FALSE,
          ct.cell.size = ct_cell_size
        )$sc.eset.qc
      } else {
        sc_obj <- SCDC::SCDC_qc_ONE(sc_obj, ct_varname, sample, ct_sub,
          iter_max = iter_max,
          nu = nu, epsilon = epsilon, weight.basis = weight_basis,
          qcthreshold = qcthreshold, generate.figure = FALSE,
          ct.cell.size = ct_cell_size
        )$sc.eset.qc
      }
    }
    return(sc_obj)
  }, single_cell_object, cell_type_annotations, batch_ids)
  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)

  # If multiple sc sets are supplied, the _ENSEMBLE method is used, otherwise the _prop one
  if (length(sc_eset) == 1) {
    # The SCDC_prop method and the SCDC_ENSEMBLE method have different default values for iter_max
    # and epsilon. That is my solution to indicate the different default values in the method
    # header while still supplying them to the final method
    if (is.null(iter_max)) {
      iter_max <- 1000
    }
    if (is.null(epsilon)) {
      epsilon <- 0.01
    }
    sc_eset <- sc_eset[[1]]

    # SCDC always supplies two types of each method, one normal one and one if all the single
    # cells are from the same individual (_ONE)
    if (length(unique(sc_eset@phenoData@data[, sample])) > 1) {
      return(SCDC::SCDC_prop(bulk_eset, sc_eset,
        ct.varname = ct_varname, sample = sample,
        ct.sub = ct_sub, iter.max = iter_max, nu = nu, epsilon = epsilon,
        truep = truep, weight.basis = weight_basis, ct.cell.size = ct_cell_size,
        Transform_bisque = transform_bisque
      ))
    } else {
      return(SCDC::SCDC_prop_ONE(bulk_eset, sc_eset,
        ct.varname = ct_varname, sample = sample,
        ct.sub = ct_sub, iter.max = iter_max, nu = nu, epsilon = epsilon,
        truep = truep, weight.basis = weight_basis, ct.cell.size = ct_cell_size,
        Transform_bisque = transform_bisque
      ))
    }
  } else {
    # Using the normal defaults of the method
    # The SCDC_prop method and the SCDC_ENSEMBLE method have different default values for iter_max
    # and epsilon. That is my solution to indicate the different default values in the method
    # header while still supplying them to the final method
    if (is.null(iter_max)) {
      iter_max <- 2000
    }
    if (is.null(epsilon)) {
      epsilon <- 0.001
    }

    if (!is.null(names_sc_objects)) {
      if (length(sc_eset) != length(names_sc_objects)) {
        stop(
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
      ct.cell.size = ct_cell_size, Transform_bisque = transform_bisque,
      grid.search = grid_search, search.length = search_length
    ))
  }
}
