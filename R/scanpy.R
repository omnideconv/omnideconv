#' Install scanpy into python environment.
#'
#' @param method
#' @param conda
#'
#' @return
#' @export
#'
#' @examples
install_scanpy <- function(method = "auto", conda = "auto") {
  reticulate::py_install("scanpy", method = method, conda = conda)
}

#' Read .h5ad files into R with scanpy
#'
#' @param path Path to .h5ad file
#'
#' @return AnnData Object
#' @export
#'
#' @examples
read_h5ad <- function(path){
  scanpy_checkload()
  file <- scanpy$read(path)
  return(file)
}

#' Title
#'
#' @param AnnData AnnData Object
#' @param path Path to where .h5ad file should be saved
#'
#' @return
#' @export
#'
#' @examples
write_h5ad <- function(AnnData,path){
  scanpy_checkload()
  AnnData$write(path)
}

#' Checks if scanpy is installed
#'
#' @return
#' @export
#'
#' @examples
scanpy_checkload <- function(){
  if (reticulate::py_module_available("scanpy")){
    reticulate::import("scanpy")
  }
  else{
    base::stop("python module scanpy not available in environment! Run install_scanpy() to install it.")
  }
}
