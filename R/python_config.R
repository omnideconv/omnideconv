#' Defines python path.
#'
#' Python3 is needed for scanpy and scaden to work.
#'
#' @param path_to_python_binaries Path to python binaries
#'
#' @export
#'
#' @examples
set_python <- function(path_to_python_binaries){
  reticulate::use_python(python = path_to_python_binaries)
}


#' Initiates python environment
#'
#' @param python (optional) If own python should be used please indicate it's binaries
#'
#' @export
#'
#' @examples
init_python <- function(python=NULL){
  if (!reticulate::py_available()){
    if (is.null(python)){
      if(!dir.exists(reticulate::miniconda_path())){
        base::message("Setting up miniconda environment..")
        suppressMessages(reticulate::install_miniconda())
      }
      reticulate::use_miniconda(condaenv = "r-reticulate",required = TRUE)
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

#' Checks if python is available in environment
#'
#' @return boolean
#' @export
#'
#' @examples
python_available<- function(){
  return(reticulate::py_available())
}




