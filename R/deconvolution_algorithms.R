source("DWLS.R")
source("decon_dwls.R")

#' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `bisque`, `dwls`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods = c("Bisque"="bisque", "DWLS" = "dwls")


#' Building the signature matrix for dwls
#'
#' @param single_cell_object A matrix or dataframe of single cell data, where genes as rows and cells as columns.
#' @param gene_annotations A Vector of genes' id.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order as cells in sc_data.
#' @param method The deconvoluation_method
#' @param result_path A string of the path, where the result should be saved
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
#' @examples
build_signature <- function(single_cell_object, gene_annotations, cell_type_annotations, method = deconvolution_methods, path = result_path){

  if (class(single_cell_object)[[1]]!="matrix")
    single_cell_object <- as.matrix(single_cell_object)
  if(method == "bisque"){
    sc_eset <- get_single_cell_expression_set(single_cell_object, colnames(single_cell_object), rownames(single_cell_object), cell_type_annotations)
  }

  if(method=="dwls"){
    cell_type_annotations<-str_replace(cell_type_annotations, " ", "_")
    row.names(single_cell_object)<-gene_annotations
    colnames(single_cell_object)<-cell_type_annotations
  }

  signature <- switch(method,
                      bisque = BisqueRNA::GenerateSCReference(sc_eset,"cellType"),
                      dwls = buildSignatureMatrixMAST(single_cell_object, cell_type_annotations, path),
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
#'   Supported methods are \"bisque\", \"...\"
#' @param ... Additional parameters, passed to the algorithm used.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are individuals, columns are cell types.
#' @export
#'
#' @examples
deconvolute <- function(bulk_gene_expression, signature, method = deconvolution_methods, dwls_submethod = dwls_submethod){

  if (class(bulk_gene_expression)[[1]]!="matrix")
    bulk_gene_expression <- base::as.matrix(bulk_gene_expression)

  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)


  deconv <- switch(method,
                   bisque = {
                     #Necessary for bisque, because bisqueReferenceDecomp needs to access internal bisque-package methods
                     base::environment(bisque_reference_decomp) <- base::environment(BisqueRNA::SimulateData)
                     bisque_reference_decomp(bulk_eset, signature, ...)$bulk.props
                   },
                   dwls = deconvolute_dwls(bulk_gene_expression, signature, dwls_submethod),


  )
  return(deconv)
}
