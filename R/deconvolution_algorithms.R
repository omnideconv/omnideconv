#' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `bisque`, `dwls`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods = c("Bisque"="bisque", "MOMF"="momf", "DWLS" = "dwls", "Scaden"="scaden")


#' Building the signature matrix
#'
#' The single_cell_object is expected to have rownames() and colnames()
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes, columns are samples. Row and column names need to be set.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object
#' @param method A string specifying the method.
#'   Supported methods are \"bisque\", \"momf\", \"dwls\", \"...\"
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples. Necessary for MOMF, defaults to NULL.Row and column names need to be set.
#' @param verbose Whether the algorithms should print out what they are doing.
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
#' @examples
build_model <- function(single_cell_object, cell_type_annotations, method = deconvolution_methods, bulk_gene_expression = NULL, verbose = TRUE, ...){

  if (class(single_cell_object)[[1]]!="matrix")
    single_cell_object <- as.matrix(single_cell_object)


  sc_eset <- get_single_cell_expression_set(single_cell_object, colnames(single_cell_object), rownames(single_cell_object), cell_type_annotations)


  signature <- switch(tolower(method),
                      bisque = BisqueRNA::GenerateSCReference(sc_eset,"cellType"),
                      #momf needs bulk set and signature matrix containing the same genes
                      momf = {
                        if (is.null(bulk_gene_expression)){
                          base::stop("'bulk_gene_expression' argument is required for MOMF")
                        }
                        MOMF::momf.computeRef(single_cell_object[intersect(rownames(single_cell_object), rownames(bulk_gene_expression)),], cell_type_annotations)
                      },
                      scaden = {
                        if (is.null(bulk_gene_expression)){
                          base::stop("'bulk_gene_expression' argument is required for Scaden")
                        }
                        scaden_build_model(single_cell_object,cell_type_annotations, bulk_data = bulk_gene_expression, verbose = verbose, ...)
                      },
                      dwls = buildSignatureMatrixMAST(as.data.frame(single_cell_object), cell_type_annotations, path = NULL, verbose = verbose, ...)
  )

  return(signature)
}


#' Deconvolution
#'
#' @param bulk_gene_expression A matrix or dataframe with the bulk data. Rows are genes, columns are samples.
#' @param signature The signature matrix.
#' @param method A string specifying the method.
#'   Supported methods are \"bisque\", \"momf\", \"dwls\", \"...\"
#' @param single_cell_object Needed for deconvolution with MOMF and Bisque. Defaults to NULL.
#' @param cell_type_annotations Needed for deconvolution with Bisque. Defaults to NULL.
#' @param verbose Whether the algorithms should print out what they are doing.
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are individuals, columns are cell types.
#' @export
#'
#' @examples
deconvolute <- function(bulk_gene_expression, signature, method = deconvolution_methods, single_cell_object = NULL, cell_type_annotations = NULL, verbose = TRUE, ...){

  if (class(bulk_gene_expression)[[1]]!="matrix")
    bulk_gene_expression <- base::as.matrix(bulk_gene_expression)

  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)

  deconv <- switch(tolower(method),
                   bisque = {
                     #Necessary for bisque, because bisqueReferenceDecomp needs to access internal bisque-package methods
                     base::environment(bisque_reference_decomp) <- base::environment(BisqueRNA::SimulateData)
                     bisque_reference_decomp(bulk_eset, signature, single_cell_object,
                                             cell_type_annotations, verbose = verbose)$bulk.props
                   },
                   momf=deconvolute_MOMF(bulk_gene_expression, signature, single_cell_object, verbose = verbose, ...),
                   scaden = scaden_deconvolute(signature, bulk_gene_expression, verbose = verbose, ...),
                   dwls = deconvolute_dwls(bulk_gene_expression, signature, verbose = verbose, ...)
  )

  #Alphabetical order of celltypes
  deconv <- deconv[,order(colnames(deconv))]

  return(deconv)
}

