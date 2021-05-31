#' Calculates the signature model with AutoGeneS
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes, columns are samples. Row and column names need to be set.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples. Necessary for MOMF, defaults to NULL.Row and column names need to be set.
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
build_model_autogenes <- function(single_cell_object, cell_type_annotations, bulk_gene_expression = NULL){
  if (!is.null(bulk_gene_expression)){
    single_cell_object <- single_cell_object[intersect(rownames(single_cell_object),rownames(bulk_gene_expression)),]
  }
  saved_h5ad <- save_as_h5ad(single_cell_object,cell_type_annotations)

  path_to_build_model_python_script <- paste0(here(),"/R/train_model.py")
  parameters <- get_buildmodel_parameters()

  system(paste("python3",path_to_build_model_python_script,parameters), ignore.stdout = !verbose, ignore.stderr = !verbose)
  return()
}

#' Deconvolution Analysis using AutoGeneS
#'
#' @param bulk_gene_expression Dataframe or matrix of bulk RNA-seq data (genes x individuals)
#' @param signature Signature Matrix (genes x individuals from scRNA-seq)
#' @param single_cell_object scRNA-seq Object (genes x cells)
#' @param method Determines which divergence to use. Options: Kullback-Leibler "KL", Itakura-Saito "IS". Defaults to "KL"
#' @param verbose Whether the algorithm should print out what it is doing.
#' @param ... additional parameters
#'
#' @return cell proportion matrix
#' @export
#'
deconvolute_autogenes <- function(bulk_gene_expression, signature, single_cell_object,
                             verbose = FALSE, method = "KL", ...){
}




#' Install AutoGeneS
#'
#' Install AutoGeneS into a previously defined python environment.
#' A python environment can be defined by set_virtualenv() or set_python() command.
#' Alternatively a new environment can be created via create_virtualenv() method.
#'
install_autogenes <- function() {
  reticulate::py_install("autogenes",pip=TRUE)

}

#' Checks if scaden is installed.
#'
#' If it is available, the python module is imported.
#'
autogenes_checkload <- function(){
  if (!python_available()){
    base::message("Setting up python environment..")
    init_python()
    if (!python_available()){
      base::stop("Could not initiate miniconda python environment. Please set up manually with init_python(python=your/python/version)")
    }
  }
  if (!reticulate::py_module_available("autogenes")){
    install_autogenes()
  }
  reticulate::import("autogenes")
}

get_buildmodel_parameters <- function(){

}
