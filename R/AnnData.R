anndata_load <<- F

#' Create anndata object from single cell expression matrix, cell type labels and gene symbols
#'
#' @param x single cell expression matrix, rows = cells, columns = genes
#' @param obs vector of cell type labels
#' @param var vector of gene symbols
#' @param obsm
#' @param varm
#'
#' @return
#' @export
#'
#' @examples
build_anndata <- function(x, obs, var, obsm=NULL, varm=NULL){

  anndata_check_load()

  x <- as.matrix(x)

  ad <- anndata::AnnData(X=x, obs = obs, var = var, obsm=obsm, varm=varm)
  return(ad)
}

#' Read anndata object (.h5ad format)
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_anndata <- function(path){
  anndata_check_load()
  data <- anndata::read_h5ad(path)
  return(data)
}

#' Write anndata object (.h5ad format)
#'
#' @param data
#' @param path
#'
#' @return
#' @export
#'
#' @examples
write_anndata <- function(data,path){
  anndata_check_load()
  o <- data$write_h5ad(path)
}

#' Checks if anndata package is loaded
#'
#' If called and python environment is not set up, this is realized. Else, it checks if the anndata package is loaded, and if not, it does this.
#'
#' @return
#' @export
#'
#' @examples
anndata_check_load <- function(){
  if (!python_available()){
    init_python()
  }

  if (!anndata_load){
    anndata::install_anndata()
    anndata_load <<- T
  }
}


