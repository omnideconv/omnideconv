#' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `Bisque`, `MOMF`, `DWLS`, `Scaden`, `CibersortX`, `AutoGeneS`, `MuSiC`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods <- c(
  "Bisque" = "bisque", "MOMF" = "momf", "DWLS" = "dwls",
  "Scaden" = "scaden", "CibersortX" = "cibersortx",
  "AutoGeneS" = "autogenes", "MuSiC" = "music", "SCDC" = "scdc", "CPM" = "cpm",
  "BSEQ-sc" = "bseqsc"
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
#'   SingleCellExperiment: colData object)
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param method A string specifying the method.
#'   Supported methods are 'bisque', 'momf', 'dwls', 'scaden', 'cibersortx' and 'autogenes'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples. Necessary
#'   for MOMF, defaults to NULL. Row and column names need to be set
#' @param verbose Whether the algorithms should print out what they are doing
#' @param cell_type_column_name Name of the column in (Anndata: obs, SingleCellExperiment: colData),
#'   that contains the cell-type labels. Is only used if no cell_type_annotations vector is
#'   provided.
#' @param markers Named list of cell type marker genes. This parameter is only used by BSEQ-sc.
#'   The type of gene identifiers (names(markers)) must be the same as the ones used as feature/row
#'   names in the single_cell_object.
#' @param ... Additional parameters, passed to the algorithm used
#'
#' @return The signature matrix. Rows are genes, columns are cell types
#' @export
#'
#' @examples
build_model <- function(single_cell_object, cell_type_annotations = NULL, batch_ids = NULL,
                        method = deconvolution_methods, bulk_gene_expression = NULL, verbose = TRUE,
                        cell_type_column_name = NULL, markers = NULL, ...) {
  method <- tolower(method)
  check_and_install(method)

  if (class(single_cell_object)[[1]] == "AnnDataR6") {
    single_cell_object <- anndata_to_singlecellexperiment(single_cell_object)
  }

  if (class(single_cell_object)[[1]] == "SingleCellExperiment") {
    matrix_and_annotation <-
      singlecellexperiment_to_matrix(single_cell_object,
        cell_type_column_name = cell_type_column_name
      )
    single_cell_object <- matrix_and_annotation$matrix
    if (is.null(cell_type_annotations)) {
      if (is.null(cell_type_column_name)) {
        base::stop(
          "Either provide cell type annotations as vector (cell_type_annotations) or the ",
          "name of the column that stores label information!"
        )
      } else {
        cell_type_annotations <- matrix_and_annotation$annotation_vector
      }
    }
  }

  if (class(single_cell_object)[[1]] != "matrix") {
    single_cell_object <- as.matrix(single_cell_object)
  }

  cell_type_annotations <- escape_blanks(cell_type_annotations)
  rownames(single_cell_object) <- escape_blanks(rownames(single_cell_object))
  colnames(single_cell_object) <- escape_blanks(colnames(single_cell_object))
  if (!is.null(bulk_gene_expression)) {
    rownames(bulk_gene_expression) <- escape_blanks(rownames(bulk_gene_expression))
    colnames(bulk_gene_expression) <- escape_blanks(colnames(bulk_gene_expression))
  }
  if (!is.null(batch_ids)) {
    batch_ids <- escape_blanks(batch_ids)
  }
  if (!is.null(markers)) {
    names(markers) <- escape_blanks(names(markers))
  }

  signature <- switch(method,
    bisque = build_model_bisque(single_cell_object, cell_type_annotations, batch_ids,
      verbose = verbose, ...
    ),
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
    music = build_model_music(),
    scdc = build_model_scdc(),
    cpm = build_model_cpm(),
    bseqsc = build_model_bseqsc(single_cell_object, cell_type_annotations, markers, batch_ids, ...)
  )



  # Only do if it is a matrix or dataframe
  if ("matrix" %in% class(signature) || "data.frame" %in% class(signature)) {
    rownames(signature) <- deescape_blanks(rownames(signature))
    colnames(signature) <- deescape_blanks(colnames(signature))
  }

  return(signature)
}


#' Deconvolution
#'
#' @param bulk_gene_expression A matrix or dataframe with the bulk data. Rows are genes, columns
#'   are samples.
#' @param signature The signature matrix.
#' @param method A string specifying the method.
#'   Supported methods are 'bisque', 'momf', 'dwls', 'scaden', 'cibersortx' and 'autogenes'
#' @param single_cell_object Needed for deconvolution with MOMF and Bisque. Defaults to NULL.
#'   Alternatively a SingleCellExperiment or an AnnData object can be provided. In that case, note
#'   that cell-type labels need to be indicated either directly providing a vector
#'   (cell_type_annotations) or by indicating the column name that indicates the cell-type labels
#'   (cell_type_column_name). (Anndata: obs object, SingleCellExperiment: colData object)
#' @param cell_type_annotations Needed for deconvolution with Bisque, MuSiC and SCDC.
#'   Defaults to NULL.
#' @param batch_ids A vector of the ids of the samples or individuals. Defaults to NULL.
#' @param verbose Whether the algorithms should print out what they are doing.
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
deconvolute <- function(bulk_gene_expression, signature, method = deconvolution_methods,
                        single_cell_object = NULL, cell_type_annotations = NULL, batch_ids = NULL,
                        cell_type_column_name = NULL, verbose = FALSE, ...) {
  method <- tolower(method)
  check_and_install(method)

  if (class(single_cell_object)[[1]] == "AnnDataR6") {
    single_cell_object <- anndata_to_singlecellexperiment(single_cell_object)
  }

  if (class(single_cell_object)[[1]] == "SingleCellExperiment") {
    matrix_and_annotation <-
      singlecellexperiment_to_matrix(single_cell_object,
        cell_type_column_name = cell_type_column_name
      )
    single_cell_object <- matrix_and_annotation$matrix
    if (is.null(cell_type_annotations)) {
      if (is.null(cell_type_column_name)) {
        base::stop(
          "Either provide cell type annotations as vector (cell_type_annotations) or the ",
          "name of the column that stores label information!"
        )
      } else {
        cell_type_annotations <- matrix_and_annotation$annotation_vector
      }
    }
  }


  if (class(bulk_gene_expression)[[1]] != "matrix") {
    bulk_gene_expression <- base::as.matrix(bulk_gene_expression)
  }


  rownames(bulk_gene_expression) <- escape_blanks(rownames(bulk_gene_expression))
  colnames(bulk_gene_expression) <- escape_blanks(colnames(bulk_gene_expression))
  # Only do if it is a matrix or dataframe
  if ("matrix" %in% class(signature) || "data.frame" %in% class(signature)) {
    rownames(signature) <- escape_blanks(rownames(signature))
    colnames(signature) <- escape_blanks(colnames(signature))
  }
  if (!is.null(single_cell_object)) {
    if ("matrix" %in% class(single_cell_object) || "data.frame" %in% class(single_cell_object)) {
      rownames(single_cell_object) <- escape_blanks(rownames(single_cell_object))
      colnames(single_cell_object) <- escape_blanks(colnames(single_cell_object))
    } else if ("list" %in% class(single_cell_object)) {
      single_cell_object <- lapply(single_cell_object, function(sc) {
        rownames(sc) <- escape_blanks(rownames(sc))
        colnames(sc) <- escape_blanks(colnames(sc))
        sc
      })
    }
  }
  if (!is.null(cell_type_annotations)) {
    if ("character" %in% class(cell_type_annotations)) {
      cell_type_annotations <- escape_blanks(cell_type_annotations)
    } else if ("list" %in% class(cell_type_annotations)) {
      cell_type_annotations <- lapply(cell_type_annotations, escape_blanks)
    }
  }

  if (!is.null(batch_ids)) {
    if ("character" %in% class(batch_ids)) {
      batch_ids <- escape_blanks(batch_ids)
    } else if ("list" %in% class(batch_ids)) {
      batch_ids <- lapply(batch_ids, escape_blanks)
    }
  }

  deconv <- switch(method,
    bisque = deconvolute_bisque(bulk_gene_expression, signature, single_cell_object,
      cell_type_annotations, batch_ids,
      verbose = verbose, ...
    )$bulk_props,
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
        base::message(
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
    )$coefficients)
  )

  if (!is.null(deconv)) {
    # Alphabetical order of celltypes
    deconv <- deconv[, order(colnames(deconv))]
    rownames(deconv) <- deescape_blanks(rownames(deconv))
    colnames(deconv) <- deescape_blanks(colnames(deconv))
  }
  return(deconv)
}


#' The dependencies for each method
#'
required_packages <- list(
  "bisque" = c("BisqueRNA", "limSolve"),
  "momf" = c("grst/MOMF"),
  "dwls" = c("quadprog", "reshape", "e1071", "ROCR", "varhandle", "MAST", "magrittr"),
  "scaden" = c("reticulate"),
  "cibersortx" = c(),
  "autogenes" = c("reticulate"),
  "music" = c("xuranw/MuSiC"),
  "scdc" = c("grst/SCDC"),
  "cpm" = c("scBio", "dotCall64"),
  "bseqsc" = c("shenorrlab/bseqsc")
)

#' Checking and installing all dependencies for the specific methods
#'
#' @param method The name of the method that is used
check_and_install <- function(method) {
  if (!(method %in% deconvolution_methods)[[1]]) {
    base::stop(
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
  sapply(cran_pkgs, function(pkgname) {
    if (!requireNamespace(pkgname, quietly = TRUE)) {
      if (!repositories_set) {
        utils::setRepositories(graphics = FALSE, ind = c(1, 2, 3, 4, 5))
        repositories_set <- TRUE
      }
      utils::install.packages(pkgname)
    }
  })
  sapply(github_pkgs, function(pkgname) {
    bare_pkgname <- sub(".*?/", "", pkgname)
    if (!requireNamespace(bare_pkgname, quietly = TRUE)) {
      if (!repositories_set) {
        utils::setRepositories(graphics = FALSE, ind = c(1, 2, 3, 4, 5))
        repositories_set <- TRUE
      }
      remotes::install_github(pkgname)
    }
  })
}
