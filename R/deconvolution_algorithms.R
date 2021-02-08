#' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `Bisque`, `MOMF`, `DWLS`, `Scaden`, `CibersortX`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods = c("Bisque" = "bisque", "MOMF" = "momf", "DWLS" = "dwls",
                          "Scaden" = "scaden", "CibersortX" = "cibersortx")


#' Building the signature matrix
#'
#' The single_cell_object is expected to have rownames() and colnames()
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes, columns are samples. Row and column names need to be set.
#' Alternatively a SingleCellExperiment or an AnnData object can be provided. In that case, note that cell-type labels need to be indicated either directly providing a vector (cell_type_annotations)
#' or by indicating the column name that indicates the cell-type labels (cell_type_column_name). (Anndata: obs object, SingleCellExperiment: colData object)
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object
#' @param method A string specifying the method.
#'   Supported methods are \"bisque\", \"momf\", \"dwls\", \"scaden\", \"cibersortx\", \"...\"
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples. Necessary for MOMF, defaults to NULL.Row and column names need to be set.
#' @param verbose Whether the algorithms should print out what they are doing.
#' @param cell_type_column_name Name of the column in (Anndata: obs , SingleCellExperiment: colData), that contains the cell-type labels. Is only used if no cell_type_annotations vector is provided.
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
#' @examples

build_model <- function(single_cell_object, cell_type_annotations = NULL, method = deconvolution_methods, bulk_gene_expression = NULL, verbose = TRUE, cell_type_column_name = NULL,...){


  if (class(single_cell_object)[[1]]=="AnnDataR6"){
    single_cell_object <- anndata_to_singlecellexperiment(ad)
  }

  if (class(single_cell_object)[[1]]=="SingleCellExperiment"){
    matrix_and_annotation <- singlecellexperiment_to_matrix(single_cell_object,cell_type_column_name = cell_type_column_name)
    single_cell_object <- matrix_and_annotation$matrix
    if (is.null(cell_type_annotations)){
      if (is.null(cell_type_column_name)){
        base::stop("Either provide cell type annotations as vector (cell_type_annotations) or the name of the column that stores label information!")
      }
      else{
        cell_type_annotations <- matrix_and_annotation$annotation_vector
      }
    }
  }

  if (class(single_cell_object)[[1]]!="matrix"){
    single_cell_object <- as.matrix(single_cell_object)
  }


  signature <- switch(tolower(method),
                      bisque = build_model_bisque(single_cell_object,cell_type_annotations, ...),
                      #momf needs bulk set and signature matrix containing the same genes
                      momf = build_model_momf(single_cell_object,cell_type_annotations,bulk_gene_expression, ...),
                      scaden = build_model_scaden(single_cell_object,cell_type_annotations, bulk_gene_expression, verbose = verbose, ...),
                      dwls = build_model_dwls(as.data.frame(single_cell_object), cell_type_annotations, path = NULL, verbose = verbose, ...),
                      cibersortx = build_model_cibersortx(single_cell_object,cell_type_annotations,verbose = verbose, ...)
  )

  return(signature)
}


#' Deconvolution
#'
#' @param bulk_gene_expression A matrix or dataframe with the bulk data. Rows are genes, columns are samples.
#' @param signature The signature matrix.
#' @param method A string specifying the method.
#'   Supported methods are \"bisque\", \"momf\", \"dwls\", \"scaden\", \"cibersortx\", \"...\"
#' @param single_cell_object Needed for deconvolution with MOMF and Bisque. Defaults to NULL.
#' @param cell_type_annotations Needed for deconvolution with Bisque. Defaults to NULL.
#' @param verbose Whether the algorithms should print out what they are doing.
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are individuals, columns are cell types.
#' @export
#'
#' @examples
deconvolute <- function(bulk_gene_expression, signature, method = deconvolution_methods, single_cell_object = NULL, cell_type_annotations = NULL, verbose = FALSE, ...){

  if (class(bulk_gene_expression)[[1]]!="matrix")
    bulk_gene_expression <- base::as.matrix(bulk_gene_expression)

  if (! "character" %in% class(signature)){
    colnames(signature) <- make.names(colnames(signature))
  }

  deconv <- switch(tolower(method),
                   bisque = {
                     cell_type_annotations <- make.names(cell_type_annotations)
                     bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)
                     #Necessary for bisque, because bisqueReferenceDecomp needs to access internal bisque-package methods
                     base::environment(deconvolute_bisque) <- base::environment(BisqueRNA::SimulateData)
                     deconvolute_bisque(bulk_eset, signature, single_cell_object,
                                        cell_type_annotations, verbose = verbose)$bulk.props
                   },
                   momf=deconvolute_momf(bulk_gene_expression, signature, single_cell_object, verbose = verbose, ...),
                   scaden = deconvolute_scaden(signature, bulk_gene_expression, verbose = verbose, ...),
                   dwls = deconvolute_dwls(bulk_gene_expression, signature, verbose = verbose, ...),
                   cibersortx = deconvolute_cibersortx(bulk_gene_expression, signature,verbose = verbose, ...)
  )

  #Alphabetical order of celltypes
  deconv <- deconv[,order(colnames(deconv))]
  colnames(deconv) <- gsub("\\.", " ", colnames(deconv))
  return(deconv)
}

