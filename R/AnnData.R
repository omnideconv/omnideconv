#' Create anndata object from single cell expression matrix, cell type labels and gene symbols
#'
#' @param x single cell expression matrix, rows = genes, columns = samples
#' @param obs vector of cell type labels
#' @param var vector of gene symbols
#' @param obsm Key-indexed multi-dimensional observations annotation of length #observations.
#' @param varm Key-indexed multi-dimensional variables annotation of length #variables.
#'
#' @return AnnData object
#'
build_anndata <- function(x, obs, var, obsm = NULL, varm = NULL){

  anndata_check_load()

  x <- as.matrix(x)

  ad <- anndata::AnnData(X=x, obs = obs, var = var, obsm=obsm, varm=varm)

  rownames(transPosed)<-var
  colnames(transPosed)<-obs

  return(ad)
}

#' Read anndata object (.h5ad format)
#'
#' @param path path to .h5ad formatted file
#'
#' @return AnnData object
#'
read_anndata <- function(path){
  anndata_check_load()
  data <- anndata::read_h5ad(path)
  return(data)
}

#' Write anndata object (.h5ad format)
#'
#' @param data AnnData object
#' @param path path where AnnData object should be written to (.h5ad format)
#'
write_anndata <- function(data, path){
  anndata_check_load()
  data$write_h5ad(path)
}

#' Checks if anndata package is loaded
#'
#' If called and python environment is not set up, this is realized. Else, it checks if the anndata package is loaded, and if not, it does this.
#'
anndata_check_load <- function(){
  if (!python_available()){
    init_python()
    anndata_check_load()
  }

  if (!reticulate::py_module_available("anndata")){
    anndata::install_anndata()
  }
}

anndata_is_identical <- function(a,b){
  same <- TRUE

  if (!identical(as.matrix(a$X),as.matrix(b$X))){
    base::message("X object is different")
    same <- FALSE
  }
  if (!identical(as.matrix(a$obs),as.matrix(b$obs))){
    base::message("obs object is different")
    same <- FALSE
  }
  if (!identical(as.matrix(a$var),as.matrix(b$var))){
    base::message("var object is different")
    same <- FALSE
  }
  if (!all.equal(a$uns,b$uns)){
    base::message("uns object is different")
    same <- FALSE
  }
  if (!all.equal(a$obsm,b$obsm)){
    base::message("obsm object is different")
    same <- FALSE
  }
  if (!all.equal(a$varm,b$varm)){
    base::message("varm object is different")
    same <- FALSE
  }
  if (!all.equal(a$layers,b$layers)){
    base::message("layers object is different")
    same <- FALSE
  }
  return(same)

}


