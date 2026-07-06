RECTANGLE_ENV <- "r-omnideconv-rectangle"

#' @noRd
rectangle_python <- function() {
  reticulate::conda_python(RECTANGLE_ENV)
}

#' @noRd
rectangle_checkload <- function() {
  if (!(RECTANGLE_ENV %in% reticulate::conda_list()$name)) {
    message("Creating ", RECTANGLE_ENV, " conda environment (Python 3.10)...")
    reticulate::conda_create(envname = RECTANGLE_ENV, python_version = "3.10")
  }

  python_bin <- rectangle_python()
  is_installed <- system2(
    python_bin, c("-c", shQuote("import rectanglepy")),
    stdout = FALSE, stderr = FALSE
  ) == 0

  if (!is_installed) {
    message("Installing rectanglepy into ", RECTANGLE_ENV, "...")
    # Use --ignore-installed so all deps land inside the conda env and are
    # not skipped because they happen to exist in the user's site-packages.
    ret <- system2(
      python_bin,
      c("-m", "pip", "install", "--ignore-installed", "rectanglepy"),
      stdout = FALSE, stderr = FALSE
    )
    if (ret != 0) {
      stop(
        "Failed to install rectanglepy into ", RECTANGLE_ENV,
        ". Run install_rectangle_python() with verbose output for details."
      )
    }
  }
}

#' Install Python dependencies for Rectangle
#'
#' Creates the r-omnideconv-rectangle conda environment with Python 3.10 and
#' installs rectanglepy. This is called automatically on first use but can be
#' run in advance to avoid the setup delay.
#'
#' @return Invisibly NULL. Called for its side effect of setting up the conda environment.
#' @export
install_rectangle_python <- function() {
  rectangle_checkload()
}

#' Build a Rectangle signature model
#'
#' Builds a hierarchical signature model from single-cell RNA-seq data using
#' DESeq2-based marker gene selection and optional hierarchical clustering.
#' The resulting model can be passed to [deconvolute()] to avoid rebuilding
#' the signature for every bulk dataset.
#'
#' @param single_cell_object A matrix with the single-cell data. Rows are genes,
#'   columns are cells. Row and column names must be set. Raw counts expected.
#' @param cell_type_annotations A vector of cell type labels, one per column of
#'   single_cell_object.
#' @param bulk_gene_expression Optional. A matrix of bulk data (rows = genes,
#'   columns = samples). When provided, the signature is restricted to genes
#'   present in the bulk data.
#' @param optimize_cutoffs Whether to optimize p-value and log-fold-change
#'   cutoffs via grid search. Recommended but slow. Default: TRUE.
#' @param p P-value cutoff for DE analysis. Only used if optimize_cutoffs = FALSE.
#' @param lfc Log fold-change cutoff for DE analysis. Only used if
#'   optimize_cutoffs = FALSE.
#' @param n_cpus Number of CPUs for DE analysis. NULL uses all available.
#' @param gene_expression_threshold Minimum fraction of cells that must express
#'   a gene for it to be included in DE analysis. Default: 0.5.
#' @param model_path Optional path where the signature pickle file should be
#'   saved. If NULL (default), a temporary file is used. Providing a permanent
#'   path allows the signature to be reused across R sessions by passing the
#'   same path to [deconvolute()] or [extract_signature_rectanglepy()].
#' @param verbose Whether to print progress to the console. Default: FALSE.
#' @param ... Additional arguments, silently ignored. Allows extra parameters to pass through
#'   the [deconvolute()] dispatcher without error.
#'
#' @return Path to a pickle file containing the RectangleSignatureResult.
#'   Pass this path as the `signature` argument to [deconvolute()] or
#'   [extract_signature_rectanglepy()].
#' @export
build_model_rectanglepy <- function(single_cell_object, cell_type_annotations,
                                    bulk_gene_expression = NULL,
                                    model_path = NULL,
                                    optimize_cutoffs = TRUE,
                                    p = 0.015,
                                    lfc = 1.5,
                                    n_cpus = NULL,
                                    gene_expression_threshold = 0.5,
                                    verbose = FALSE, ...) {
  rectangle_checkload()

  python_bin <- rectangle_python()
  script <- system.file("python", "rectanglepy_wrapper.py", package = "omnideconv")

  sc_h5ad <- save_as_h5ad(single_cell_object, cell_type_annotations)
  output_pickle <- if (!is.null(model_path)) model_path else tempfile(fileext = ".pkl")

  cmd_args <- c(
    script, "build_model",
    "--sc_h5ad", sc_h5ad,
    "--output_pickle", output_pickle,
    "--optimize_cutoffs", tolower(as.character(optimize_cutoffs)),
    "--p", p,
    "--lfc", lfc,
    "--gene_expression_threshold", gene_expression_threshold
  )

  if (!is.null(bulk_gene_expression)) {
    bulk_csv <- tempfile(fileext = ".csv")
    write.csv(bulk_gene_expression, file = bulk_csv)
    cmd_args <- c(cmd_args, "--bulk_csv", bulk_csv)
  }

  if (!is.null(n_cpus)) {
    cmd_args <- c(cmd_args, "--n_cpus", as.integer(n_cpus))
  }

  ret <- system2(
    python_bin, cmd_args,
    stdout = if (verbose) "" else FALSE,
    stderr = if (verbose) "" else FALSE
  )
  if (ret != 0) {
    stop("Rectangle model building failed. Set verbose = TRUE for details.")
  }

  return(output_pickle)
}

#' Extract the signature matrix from a Rectangle model
#'
#' Calls `get_signature_matrix()` on the `RectangleSignatureResult` object
#' stored in the signature pickle, returning a standard R matrix. This allows
#' inspection of the selected marker genes and comparison with signatures from
#' other methods. The matrix can be saved for later use with [saveRDS()].
#'
#' @param signature Path to a .pkl file created by [build_model_rectanglepy()].
#' @param include_mrna_bias Whether to apply mRNA bias correction when
#'   computing the signature matrix. Default: TRUE.
#' @param verbose Whether to print progress to the console. Default: FALSE.
#'
#' @return A matrix with rows = genes and columns = cell types.
#' @export
extract_signature_rectanglepy <- function(signature, include_mrna_bias = TRUE, verbose = FALSE) {
  rectangle_checkload()

  if (!file.exists(signature)) {
    stop(
      "Signature file not found: ", signature, "\n",
      "The signature argument must be a path returned by build_model_rectanglepy()."
    )
  }

  python_bin <- rectangle_python()
  script <- system.file("python", "rectanglepy_wrapper.py", package = "omnideconv")
  output_csv <- tempfile(fileext = ".csv")

  ret <- system2(
    python_bin,
    c(
      script, "get_signature_matrix",
      "--signature_pickle", signature,
      "--output_csv", output_csv,
      "--include_mrna_bias", tolower(as.character(include_mrna_bias))
    ),
    stdout = if (verbose) "" else FALSE,
    stderr = if (verbose) "" else FALSE
  )
  if (ret != 0) {
    stop("Failed to extract Rectangle signature matrix. Set verbose = TRUE for details.")
  }

  as.matrix(read.csv(output_csv, row.names = 1, check.names = FALSE))
}

#' Deconvolution with Rectangle
#'
#' Estimates cell-type fractions from bulk RNA-seq data using the Rectangle
#' method. Accepts either a pre-built signature (from [build_model_rectanglepy()])
#' or single-cell data for a one-step run.
#'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns
#'   are samples. Row and column names must be set. TPM normalized.
#' @param signature Path to a .pkl file created by [build_model_rectanglepy()].
#'   When NULL, single_cell_object and cell_type_annotations are required and
#'   the signature is built internally.
#' @param single_cell_object A matrix with the single-cell data. Rows are genes,
#'   columns are cells. Row and column names must be set. Raw counts expected.
#'   Only used when signature is NULL.
#' @param cell_type_annotations A vector of cell type labels. Only used when
#'   signature is NULL.
#' @param correct_mrna_bias Whether to correct for mRNA content bias.
#'   Default: TRUE.
#' @param optimize_cutoffs Whether to optimize DE cutoffs during signature
#'   building. Only used when signature is NULL. Default: TRUE.
#' @param p P-value cutoff for DE analysis. Only used when signature is NULL
#'   and optimize_cutoffs = FALSE.
#' @param lfc Log fold-change cutoff. Only used when signature is NULL and
#'   optimize_cutoffs = FALSE.
#' @param n_cpus Number of CPUs. NULL uses all available.
#' @param gene_expression_threshold Gene expression threshold for DE analysis.
#'   Only used when signature is NULL. Default: 0.5.
#' @param verbose Whether to print progress to the console. Default: FALSE.
#'
#' @return A matrix with rows = samples and columns = cell types (including an
#'   "Unknown" column representing uncharacterized cellular content).
#' @export
deconvolute_rectanglepy <- function(bulk_gene_expression,
                                    signature = NULL,
                                    single_cell_object = NULL,
                                    cell_type_annotations = NULL,
                                    correct_mrna_bias = TRUE,
                                    optimize_cutoffs = TRUE,
                                    p = 0.015,
                                    lfc = 1.5,
                                    n_cpus = NULL,
                                    gene_expression_threshold = 0.5,
                                    verbose = FALSE) {
  rectangle_checkload()

  python_bin <- rectangle_python()
  script <- system.file("python", "rectanglepy_wrapper.py", package = "omnideconv")

  bulk_csv <- tempfile(fileext = ".csv")
  write.csv(bulk_gene_expression, file = bulk_csv)
  output_csv <- tempfile(fileext = ".csv")

  n_cpus_args <- if (!is.null(n_cpus)) c("--n_cpus", as.integer(n_cpus)) else character(0)

  if (!is.null(signature)) {
    if (!file.exists(signature)) {
      stop(
        "Signature file not found: ", signature, "\n",
        "The signature argument must be a path returned by build_model_rectanglepy()."
      )
    }
    cmd_args <- c(
      script, "deconvolute",
      "--bulk_csv", bulk_csv,
      "--signature_pickle", signature,
      "--output_csv", output_csv,
      "--correct_mrna_bias", tolower(as.character(correct_mrna_bias)),
      n_cpus_args
    )
  } else {
    if (is.null(single_cell_object) || is.null(cell_type_annotations)) {
      stop(
        "single_cell_object and cell_type_annotations are required ",
        "when no signature is provided."
      )
    }
    sc_h5ad <- save_as_h5ad(single_cell_object, cell_type_annotations)
    cmd_args <- c(
      script, "rectangle",
      "--sc_h5ad", sc_h5ad,
      "--bulk_csv", bulk_csv,
      "--output_csv", output_csv,
      "--correct_mrna_bias", tolower(as.character(correct_mrna_bias)),
      "--optimize_cutoffs", tolower(as.character(optimize_cutoffs)),
      "--p", p,
      "--lfc", lfc,
      "--gene_expression_threshold", gene_expression_threshold,
      n_cpus_args
    )
  }

  ret <- system2(
    python_bin, cmd_args,
    stdout = if (verbose) "" else FALSE,
    stderr = if (verbose) "" else FALSE
  )
  if (ret != 0) {
    stop("Rectangle deconvolution failed. Set verbose = TRUE for details.")
  }

  result <- as.matrix(read.csv(output_csv, row.names = 1, check.names = FALSE))
  return(result)
}
