#' List of supported deconvolution methods
#'
#' The methods currently supported are
#' `AutoGeneS`, `BayesPrism`, `Bisque`, `BSeq-sc`, `CIBERSORTx`, `CPM`, `DWLS`, `MOMF`, `MuSiC`,
#' `Scaden`, `SCDC`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods <- c(
  "AutoGeneS" = "autogenes", "BayesPrism" = "bayesprism", "Bisque" = "bisque", "BSeq-sc" = "bseqsc",
  "CIBERSORTx" = "cibersortx", "CDSeq" = "cdseq", "CPM" = "cpm", "DWLS" = "dwls", "MOMF" = "momf",
  "MuSiC" = "music", "Scaden" = "scaden", "SCDC" = "scdc"
)


#' Building the signature matrix
#'
#' The single_cell_object is expected to have rownames() and colnames()
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes,
#'   columns are samples. Row and column names need to be set. Alternatively a SingleCellExperiment
#'   or an AnnData object can be provided. In that case, note that cell-type labels need to be
#'   indicated either directly providing a vector (cell_type_annotations) or by indicating the
#'   column name that indicates the cell-type labels (cell_type_column_name). (Anndata: obs object,
#'   SingleCellExperiment: colData object).
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param method A string specifying the method.
#'   Supported methods are the ones in the variable 'deconvolution_methods'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples. Necessary
#'   for MOMF and Scaden, defaults to NULL. Row and column names need to be set
#' @param verbose Whether to produce an output on the console.
#' @param cell_type_column_name Name of the column in (Anndata: obs, SingleCellExperiment: colData),
#'   that contains the cell-type labels. Is only used if no cell_type_annotations vector is
#'   provided.
#' @param markers Named list of cell type marker genes. This parameter is only used by BSeq-sc.
#'   The type of gene identifiers (names(markers)) must be the same as the ones used as feature/row
#'   names in the single_cell_object.
#' @param ... Additional parameters, passed to the algorithm used
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
#' @examples
#' # More examples can be found in the unit tests at tests/testthat/test-b-buildmodel.R
#' data("single_cell_data_1")
#' data("cell_type_annotations_1")
#' data("batch_ids_1")
#' data("bulk")
#'
#' single_cell_data <- single_cell_data_1[1:2000, 1:500]
#' cell_type_annotations <- cell_type_annotations_1[1:500]
#' batch_ids <- batch_ids_1[1:500]
#' bulk <- bulk[1:2000, ]
#'
#' signature_matrix_momf <- build_model(
#'   single_cell_data, cell_type_annotations, "momf",
#'   bulk_gene_expression = bulk
#' )
#' pickle_path_autogenes <- build_model(single_cell_data, cell_type_annotations, "autogenes",
#'   population_size = 500, offspring_size = 30,
#'   crossover_pb = 0.3
#' )
build_model <- function(single_cell_object, cell_type_annotations = NULL,
                        method = deconvolution_methods, batch_ids = NULL,
                        bulk_gene_expression = NULL, verbose = FALSE,
                        cell_type_column_name = NULL, markers = NULL, ...) {
  if (length(method) > 1) {
    stop(
      "Please only specify one method and not ", length(method), ": ",
      paste(method, collapse = ", ")
    )
  }
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)
  check_and_install(method)


  # Converting all other data types into a matrix
  matrix_and_annotation <- convert_to_matrix(
    single_cell_object, cell_type_annotations,
    cell_type_column_name
  )
  single_cell_object <- matrix_and_annotation$matrix
  cell_type_annotations <- matrix_and_annotation$cell_type_annotations

  if (is.null(rownames(single_cell_object))) {
    stop("The single cell object does not have any rownames!")
  }
  if (is.null(colnames(single_cell_object))) {
    stop("The single cell object does not have any colnames!")
  }

  # Check the input data for problems like different numbers of cells in the object and the
  # annotation or strings in the data
  check_data(single_cell_object, cell_type_annotations, bulk_gene_expression)

  cell_type_annotations <- escape_special_chars(cell_type_annotations)
  rownames(single_cell_object) <- escape_special_chars(rownames(single_cell_object))
  colnames(single_cell_object) <- escape_special_chars(colnames(single_cell_object))
  if (!is.null(bulk_gene_expression)) {
    rownames(bulk_gene_expression) <- escape_special_chars(rownames(bulk_gene_expression))
    colnames(bulk_gene_expression) <- escape_special_chars(colnames(bulk_gene_expression))
  }
  if (!is.null(batch_ids)) {
    batch_ids <- escape_special_chars(batch_ids)
  }
  if (!is.null(markers)) {
    names(markers) <- escape_special_chars(names(markers))
  }

  signature <- switch(method,
    bisque = build_model_bisque(),
    # momf needs bulk set and signature matrix containing the same genes
    momf = build_model_momf(single_cell_object, cell_type_annotations, bulk_gene_expression, ...),
    scaden = build_model_scaden(single_cell_object, cell_type_annotations, bulk_gene_expression,
      verbose = verbose, ...
    ),
    dwls = build_model_dwls(as.data.frame(single_cell_object), cell_type_annotations,
      path = NULL,
      verbose = verbose, ...
    ),
    cibersortx = build_model_cibersortx(single_cell_object, cell_type_annotations,
      verbose = verbose, ...
    ),
    autogenes = build_model_autogenes(single_cell_object, cell_type_annotations,
      verbose = verbose, ...
    ),
    music = build_model_music(single_cell_object, cell_type_annotations, batch_ids, ...)$Disgn.mtx,
    scdc = build_model_scdc(single_cell_object, cell_type_annotations, batch_ids, ...)$basis,
    cpm = build_model_cpm(),
    bseqsc = build_model_bseqsc(single_cell_object, cell_type_annotations, markers, batch_ids, ...),
    cdseq = build_model_cdseq(),
    bayesprism = build_model_bayesprism()
  )



  # Only do if it is a matrix or dataframe
  if ("matrix" %in% class(signature) || "data.frame" %in% class(signature)) {
    rownames(signature) <- deescape_special_chars(rownames(signature))
    colnames(signature) <- deescape_special_chars(colnames(signature))
  }

  return(signature)
}


#' Deconvolution
#'
#' @param bulk_gene_expression A matrix or dataframe with the bulk data. Rows are genes, columns
#'   are samples.
#' @param signature The signature matrix.
#' @param method A string specifying the method.
#'   Supported methods are 'bisque', 'momf', 'dwls', 'scaden', 'cibersortx', 'autogenes' and 'bayesprism'
#' @param single_cell_object Needed for deconvolution with MOMF, Bisque and BayesPrism. Defaults to NULL.
#'   Alternatively a SingleCellExperiment or an AnnData object can be provided. In that case, note
#'   that cell-type labels need to be indicated either directly providing a vector
#'   (cell_type_annotations) or by indicating the column name that indicates the cell-type labels
#'   (cell_type_column_name). (Anndata: obs object, SingleCellExperiment: colData object)
#' @param cell_type_annotations Needed for deconvolution with Bisque, MuSiC, SCDC and BayesPrism
#'   Defaults to NULL.
#' @param batch_ids A vector of the ids of the samples or individuals. Defaults to NULL.
#' @param verbose Whether to produce an output on the console.
#' @param ... Additional parameters, passed to the algorithm used.
#' @param cell_type_column_name Name of the column in (Anndata: obs, SingleCellExperiment: colData),
#'   that contains the cell-type labels. Is only used if no cell_type_annotations vector
#'   is provided.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are
#' individuals, columns are cell types.
#' @export
#'
#' @examples
#' # More examples can be found in the unit tests at tests/testthat/test-c-deconvolute.R
#' data("single_cell_data_1")
#' data("cell_type_annotations_1")
#' data("batch_ids_1")
#' data("bulk")
#'
#' single_cell_data <- single_cell_data_1[1:2000, 1:500]
#' cell_type_annotations <- cell_type_annotations_1[1:500]
#' batch_ids <- batch_ids_1[1:500]
#' bulk <- bulk[1:2000, ]
#'
#' signature_matrix_autogenes <- build_model(
#'   single_cell_data, cell_type_annotations, "autogenes",
#'   batch_ids
#' )
#' deconv_autogenes <- deconvolute(
#'   bulk, signature_matrix_autogenes, "autogenes"
#' )
#'
#' deconv_bisque <- deconvolute(
#'   bulk, NULL, "bisque", single_cell_data,
#'   cell_type_annotations, batch_ids
#' )
deconvolute <- function(bulk_gene_expression, signature, method = deconvolution_methods,
                        single_cell_object = NULL, cell_type_annotations = NULL, batch_ids = NULL,
                        cell_type_column_name = NULL, verbose = FALSE, ...) {
  if (length(method) > 1) {
    stop(
      "Please only specify one method and not ", length(method), ": ",
      paste(method, collapse = ", ")
    )
  }
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)
  check_and_install(method)


  # Converting all other data types into a matrix
  matrix_and_annotation <- convert_to_matrix(
    single_cell_object, cell_type_annotations,
    cell_type_column_name
  )
  single_cell_object <- matrix_and_annotation$matrix
  cell_type_annotations <- matrix_and_annotation$cell_type_annotations


  # Converting all other data types into a matrix
  bulk_gene_expression <- convert_to_matrix(bulk_gene_expression, "bulk")$matrix

  # Check the input data for problems like different numbers of cells in the object and the
  # annotation or strings in the data
  check_data(single_cell_object, cell_type_annotations, bulk_gene_expression)

  rownames(bulk_gene_expression) <- escape_special_chars(rownames(bulk_gene_expression))
  colnames(bulk_gene_expression) <- escape_special_chars(colnames(bulk_gene_expression))
  # Only do if it is a matrix or dataframe
  if ("matrix" %in% class(signature) || "data.frame" %in% class(signature)) {
    rownames(signature) <- escape_special_chars(rownames(signature))
    colnames(signature) <- escape_special_chars(colnames(signature))
  }
  if (!is.null(single_cell_object)) {
    if ("matrix" %in% class(single_cell_object) || "data.frame" %in% class(single_cell_object)) {
      rownames(single_cell_object) <- escape_special_chars(rownames(single_cell_object))
      colnames(single_cell_object) <- escape_special_chars(colnames(single_cell_object))
    } else if ("list" %in% class(single_cell_object)) {
      single_cell_object <- lapply(single_cell_object, function(sc) {
        rownames(sc) <- escape_special_chars(rownames(sc))
        colnames(sc) <- escape_special_chars(colnames(sc))
        sc
      })
    }
  }
  if (!is.null(cell_type_annotations)) {
    if ("character" %in% class(cell_type_annotations)) {
      cell_type_annotations <- escape_special_chars(cell_type_annotations)
    } else if ("list" %in% class(cell_type_annotations)) {
      cell_type_annotations <- lapply(cell_type_annotations, escape_special_chars)
    }
  }

  if (!is.null(batch_ids)) {
    if ("character" %in% class(batch_ids)) {
      batch_ids <- escape_special_chars(batch_ids)
    } else if ("list" %in% class(batch_ids)) {
      batch_ids <- lapply(batch_ids, escape_special_chars)
    }
  }

  if (verbose && method %in% c("bisque", "music", "scdc", "cpm", "cdseq", "bayesprism") && !is.null(signature)) {
    message(
      "A signature was provided, even though you chose a method that does not use ",
      "an external one."
    )
  }

  deconv <- switch(method,
    bisque = t(deconvolute_bisque(bulk_gene_expression, single_cell_object, cell_type_annotations,
      batch_ids,
      verbose = verbose, ...
    )$bulk.props),
    momf = deconvolute_momf(bulk_gene_expression, signature, single_cell_object,
      verbose = verbose, ...
    )$cell.prop,
    scaden = deconvolute_scaden(signature, bulk_gene_expression, verbose = verbose, ...),
    dwls = deconvolute_dwls(bulk_gene_expression, signature, verbose = verbose, ...),
    cibersortx = deconvolute_cibersortx(bulk_gene_expression, signature, verbose = verbose, ...),
    autogenes = deconvolute_autogenes(bulk_gene_expression, signature,
      verbose = verbose, ...
    )$proportions,
    music = deconvolute_music(bulk_gene_expression, single_cell_object, cell_type_annotations,
      batch_ids,
      verbose = verbose, ...
    )$Est.prop.weighted,
    scdc = {
      res <- deconvolute_scdc(bulk_gene_expression, single_cell_object, cell_type_annotations,
        batch_ids,
        verbose = verbose, ...
      )
      if ("prop.est.mvw" %in% names(res)) {
        res$prop.est.mvw
      } else if ("w_table" %in% names(res)) {
        SCDC::wt_prop(res$w_table, res$prop.only)
      } else {
        message(
          "There seems to be an error, as the result of deconvolute_scdc did not ",
          "contain prop.est.mvw or w_table"
        )
        res
      }
    },
    cpm = deconvolute_cpm(bulk_gene_expression, single_cell_object, cell_type_annotations,
      verbose = verbose, ...
    )$cellTypePredictions,
    bseqsc = t(deconvolute_bseqsc(bulk_gene_expression, signature,
      verbose = verbose, ...
    )$coefficients),
    cdseq = t(deconvolute_cdseq(bulk_gene_expression, single_cell_object, cell_type_annotations,
      batch_ids,
      verbose = verbose, ...
    )$cdseq_prop_merged),
    bayesprism = deconvolute_bayesprism(
      bulk_gene_expression, single_cell_object, cell_type_annotations,
      ...
    )$theta
  )

  if (!is.null(deconv)) {
    # Normalize the results
    deconv <- normalize_deconv_results(deconv)
    # Alphabetical order of celltypes
    rownames(deconv) <- deescape_special_chars(rownames(deconv))
    colnames(deconv) <- deescape_special_chars(colnames(deconv))
    deconv <- deconv[, order(colnames(deconv)), drop = FALSE]
  }
  return(deconv)
}


#' The dependencies for each method
#'
required_packages <- list(
  "autogenes" = c("reticulate"),
  "bayesprism" = c("Danko-Lab/BayesPrism/BayesPrism"),
  "bisque" = c("BisqueRNA"),
  "bseqsc" = c("shenorrlab/bseqsc"),
  "cdseq" = c("PelzKo/CDSeq_R_Package"),
  "cibersortx" = c("uuid"),
  "cpm" = c("amitfrish/scBio"),
  "dwls" = c("PelzKo/dwls"),
  "momf" = c("grst/MOMF"),
  "music" = c("xuranw/MuSiC"),
  "scaden" = c("reticulate"),
  "scdc" = c("grst/SCDC")
)

#' Checking and installing all dependencies for the specific methods
#'
#' @param method The name of the method that is used
#'
#' @importFrom utils askYesNo
check_and_install <- function(method) {
  if (!(method %in% deconvolution_methods)[[1]]) {
    stop(
      paste(
        "Method", method,
        "not recognized. Please refer to 'deconvolution_methods' for the integrated methods."
      )
    )
  }
  method <- method[[1]]
  packages <- required_packages[[method]]
  github_pkgs <- grep("^.*?/.*?$", packages, value = TRUE)
  cran_pkgs <- packages[!(packages %in% github_pkgs)]
  repositories_set <- FALSE
  package_download_allowed <- FALSE
  sapply(cran_pkgs, function(pkgname) {
    if (!requireNamespace(pkgname, quietly = TRUE)) {
      if (!repositories_set) {
        utils::setRepositories(graphics = FALSE, ind = c(1, 2, 3, 4, 5))
        repositories_set <<- TRUE
        package_download_allowed <<- askYesNo(
          paste0(
            "You requested to run ", method,
            " which is currently not installed. Do you want ",
            "to install the packages required for it: ", packages
          )
        )
        message(
          "To install the dependencies for all methods at once, run ",
          "devtools::install_github(\"omnideconv/omnideconv\", ",
          "dependencies = c(\"Imports\", \"Suggests\"))"
        )
      }
      if (package_download_allowed) {
        utils::install.packages(pkgname)
      }
    }
  })
  sapply(github_pkgs, function(pkgname) {
    bare_pkgname <- sub(".*?/", "", pkgname)
    if (bare_pkgname == "CDSeq_R_Package") {
      bare_pkgname <- "CDSeq"
    } else if (bare_pkgname == "dwls") {
      bare_pkgname <- "DWLS"
    }
    if (!requireNamespace(bare_pkgname, quietly = TRUE)) {
      if (!repositories_set) {
        utils::setRepositories(graphics = FALSE, ind = c(1, 2, 3, 4, 5))
        repositories_set <<- TRUE
        package_download_allowed <<- askYesNo(
          paste0(
            "You requested to run ", method,
            " which is currently not installed. Do you want ",
            "to install the packages required for it: ", packages
          )
        )
        message(
          "To install the dependencies for all methods at once, run ",
          "devtools::install_github(\"omnideconv/omnideconv\", ",
          "dependencies = c(\"Imports\", \"Suggests\"))"
        )
      }
      if (package_download_allowed) {
        remotes::install_github(pkgname)
      }
    }
  })
  if (repositories_set && !package_download_allowed) {
    message(
      "To install the dependencies for all methods at once, run ",
      "devtools::install_github(\"omnideconv/omnideconv\", ",
      "dependencies = c(\"Imports\", \"Suggests\"))"
    )
    stop(paste0(method, " can not be run without installing the required packages: ", packages))
  }
}


#' Calculation of the condition number
#'
#' A problem with a low condition number is said to be well-conditioned, while a problem with a
#' high condition number is said to be ill-conditioned. An ill-conditioned problem is one where, for
#' a small change in the inputs (the independent variables) there is a large change in the answer or
#' dependent variable.
#'
#' @param signature_matrix A signature matrix created with the build_model method
#'
#' @return The condition number
#' @export
#'
#'
#' @examples
#' data("single_cell_data_1")
#' data("cell_type_annotations_1")
#' data("batch_ids_1")
#' data("bulk")
#'
#' single_cell_data <- single_cell_data_1[1:2000, 1:500]
#' cell_type_annotations <- cell_type_annotations_1[1:500]
#' batch_ids <- batch_ids_1[1:500]
#' bulk <- bulk[1:2000, ]
#'
#' signature_matrix_momf <- build_model(
#'   single_cell_data, cell_type_annotations, "momf",
#'   bulk_gene_expression = bulk
#' )
#' cond_num <- calc_condition_number(signature_matrix_momf)
#' cond_num
calc_condition_number <- function(signature_matrix) {
  return(kappa(signature_matrix, exact = TRUE))
}

#' Install all python packages
#'
#' This makes sure a valid python installation exists and all needed packages are pulled and
#' installed.
#'
#' @export
#'
install_all_python <- function() {
  init_python()
  anndata_checkload()
  autogenes_checkload()
  scaden_checkload()
}
