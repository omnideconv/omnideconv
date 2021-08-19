#' Calculates the signature model with BSEQ-sc
#'
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param markers Named list of cell type marker genes.
#'   The type of gene identifiers (names(markers)) must be the same as the ones used as feature/row
#'   names in the single_cell_object.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param ct_scale logical that indicates if the single cell expresson profiles
#'   should be rescaled according to their cell type average total count.
#' @param limit logical that indicates if the returned basis matrix should only contain
#'   cell types listed in `markers`.
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
build_model_bseqsc <- function(single_cell_object, cell_type_annotations, markers, batch_ids,
                               ct_scale = TRUE, limit = TRUE) {
  if (is.null(markers) || is.null(batch_ids)) {
    stop("'markers' and 'batch_ids' argument is required for BSEQ-sc")
  }
  return(bseqsc::bseqsc_basis(single_cell_object, markers, cell_type_annotations,
    batch_ids,
    ct.scale = ct_scale, limit = limit
  ))
}

#' Deconvolution Analysis using BSEQ-sc
#'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param signature The signature matrix. Rows are genes, columns are cell types.
#'   If `NULL` then a built-in reference for pancreatic islet cell sub-population
#'   is used (data PancreasIslet in bseqsc package).
#' @param log Logical that indicates if the data `x` is in in log-scale.
#'   If `NULL`, then log scale is inferred by [xbioc::is_logscale].
#' @param ... Other arguments passed to `CIBERSORT`.
#' @param verbose Whether to produce an output on the console.
#'
#' @return A list including:
#' \item{coef}{The matrix of estimated proportions (cell type x samples).}
#' \item{stats}{The statistics computed by `CIBERSORT`.}
#'
#' @export
#'
deconvolute_bseqsc <- function(bulk_gene_expression, signature = NULL, log = NULL, ...,
                               verbose = FALSE) {
  return(bseqsc::bseqsc_proportions(bulk_gene_expression,
    reference = signature, log = log, verbose = verbose, ...
  ))
}


#' Setup BSeq-SC External Dependencies
#'
#' Configure `bseqsc` by setting up `CIBERSORT` source code.
#'
#' `BSeq-sc` uses `CIBERSORT` to estimate cell type proportions, based on reference
#'  expression profiles.
#'  Due to licensing requirements, source code for this algorithm needs to be
#'  downloaded separately from its website http://cibersort.stanford.edu.
#'  It is released under the Stanford Non-Commercial License.
#' In order to use it with `bseqsc` you will need to:
#'
#'   1. Got to http://cibersort.stanford.edu
#'   2. Register and log in
#'   3. Download the latest R source code from the [Download
#'      section](http://cibersort.stanford.edu/download.php).
#'   4. Configure `bseqsc` by pointing it to the downloaded file. This is
#'      done using the function `bseqsc_config`, which will copy the given R
#'      source file into the `R-data/bseqsc` sub-directory in the user's home
#'      directory for subsequent usage:
#'
#' ```
#' bseqsc_config('path/to/downloaded/source/CIBERSORT.R')
#' ```
#'
#' @param file Path to the CIBERSORT source R code.
#' @param error Logical that indicates if an error should be thrown
#' if configuration failed.
#'
#' @return The path where the file was copied, or `NULL` if `bseqsc` is not correctly
#' configured.
#'
#' @export
bseqsc_config <- function(file = NULL, error = FALSE) {
  return(bseqsc::bseqsc_config(file, error))
}
