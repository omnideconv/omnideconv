#' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `bisque`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export

deconvolution_methods = c("Bisque"="bisque","MOMF"="momf","Scaden"="scaden")





#' Building the signature matrix
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes, columns are samples.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object
#' @param method A string specifying the method.
#'   Supported methods are \"bisque\", \"momf\", \"...\"
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
#' @examples
build_model <- function(single_cell_object, cell_type_annotations, method = deconvolution_methods, ...){

  if (class(single_cell_object)[[1]]!="matrix")
    single_cell_object <- as.matrix(single_cell_object)

  sc_eset <- get_single_cell_expression_set(single_cell_object, colnames(single_cell_object), rownames(single_cell_object), cell_type_annotations)

  signature <- switch(method,
                      bisque = BisqueRNA::GenerateSCReference(sc_eset,"cellType"),
                      momf = MOMF::momf.computeRef(single_cell_object, cell_type_annotations)
  )

  return(signature)
}

#Input: Rows are genes, cols are samples
#Input datatype: Matrix or dataframe
#Output: Rows are genes, cols are cell types
#Output datatype: Matrix

#' Deconvolution
#'
#' @param bulk_gene_expression A matrix or dataframe with the bulk data. Rows are genes, columns are samples.
#' @param signature The signature matrix.
#' @param method A string specifying the method.
#'   Supported methods are \"bisque\", \"momf\", \"...\"
#' @param single_cell_object Needed for deconvolution with MOMF. Defaults to NULL.
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are individuals, columns are cell types.
#' @export
#'
#' @examples
deconvolute <- function(bulk_gene_expression, signature, method = deconvolution_methods, single_cell_object = NULL, ...){

  if (class(bulk_gene_expression)[[1]]!="matrix")
    bulk_gene_expression <- base::as.matrix(bulk_gene_expression)

  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)


  deconv <- switch(method,
                   bisque = {
                     #Necessary for bisque, because bisqueReferenceDecomp needs to access internal bisque-package methods
                     base::environment(bisque_reference_decomp) <- base::environment(BisqueRNA::SimulateData)
                     bisque_reference_decomp(bulk_eset, signature, ...)$bulk.props
                   },
                   momf=deconvolute_MOMF(bulk_gene_expression, signature, single_cell_object, ...)
  )
  return(deconv)
}


##Remark: Scaden now callable by one single function. Can later be split up and added to the build_model() and deconvolute() functions.
##        Then build_model() needs to inform the user to provide bulk_data. Additionally deconvolute() should then inform user that scaden can't take up a signature matrix, but only the scaden model.

#' Scaden model creation and deconvolution
#'
#' @param single_cell_object A matrix or dataframe with the single-cell expression data. Rows are cells, columns are genes.
#' @param cell_type_annotations Vector of celltype labels. Same order as rows in single_cell_object.
#' @param gene_labels Vector of gene symbols. Same order as columns in single_cell_object.
#' @param bulk_data Matrix or dataframe of bulk RNA expression. Rows are genes, columns are samples.
#' @param anndata_object Single cell expression in AnnData format. If provided, single_cell_object, cell_type_annotations and gene_labels are ignored.
#' @param ... Parameters that can be forwarded to Scaden
#'
#' @return
#' @export
#'
#' @examples
#' Calling scaden with single cell expression matrix + annotations:
#' cell_type_proportions <- scaden(single_cell_matrix ,cell_type_annotation ,gene_labels ,bulk_data)
#'
#' Calling scaden with .h5ad file:
#' cell_type_proportions <- scaden(anndata_object = h5ad_file, bulk_data = bulk_data)
#'
#' Caution: The process step of scaden fails if inputting a single_cell experiment that is log-transformed.
#'
scaden <- function(single_cell_object, cell_type_annotations, gene_labels , bulk_data, anndata_object=NULL,verbose=F, ...){
  init_python(verbose = verbose)

  if (is.null(anndata_object)){
    model <- scaden_build_model(single_cell_object, cell_type_annotations,gene_labels, bulk_data = bulk_data, ...)
  }
  else{
    model <- scaden_build_model_from_h5ad(anndata_object,bulk_data = bulk_data,...)
  }


  deconvolution <- scaden_deconvolute(model,bulk_data, ...)
  return(deconvolution)
}
