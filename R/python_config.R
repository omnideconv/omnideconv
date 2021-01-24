#' Defines python path.
#'
#' Python3 is needed for scanpy and scaden to work.
#'
#' @param path Path to the Python binaries.
#'
#' @return
#' @export
#'
#' @examples
set_python <- function(path_to_python_binaries){
  reticulate::use_python(python = path_to_python_binaries)
}


init_python <- function(python=NULL){
  if (!reticulate::py_available()){
    if (is.null(python)){
      if(!dir.exists(reticulate::miniconda_path())){
        base::message("Setting up miniconda environment..")
        suppressMessages(reticulate::install_miniconda())
      }
      reticulate::use_miniconda(condaenv = "r-reticulate",required = T)
      config <- reticulate::py_config()
      if (!python_available()){
        base::message("Python not available")
        print(config)
        base::message("Please indicate your version of python calling init_python(python=your/python)")
      }
    }
    else{
      reticulate::use_python(python= python)
      reticulate::py_config()
    }
  }
}

python_available<- function(){
  return(reticulate::py_available())
}




