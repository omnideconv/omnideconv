#' Defines a virtual environment that should be used.
#'
#'
#' This environment needs to contain python3 for scanpy and scaden to work.
#'
#' @param env_path Path to virtual environment.
#'
#' @return
#' @export
#'
#' @examples
set_virtualenv <- function(env_path){
  reticulate::use_virtualenv(env_path,required = T)
}

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

#' Creates a new virtual environment.
#'
#' Creates a virtual environment and activates it.
#' The pip version is automatically upgraded to the newest version.
#' Python3 required for scaden and scanpy to work.
#'
#'
#' @param path_to_python_binaries Path to the Python binaries.
#'
#' @return
#' @export
#'
#' @examples
create_virtualenv <- function(path_to_python_binaries){
  reticulate::virtualenv_create(python = path_to_python_binaries,envname = "r-reticulate")
  reticulate::use_virtualenv("r-reticulate",required = T)
}

init_python <- function(python=NULL,verbose=F){
  if (is.null(python)){
    if (!python_available()){
      python3 <- NULL
      disc <- reticulate::py_discover_config()
      for (python_version in disc$python_versions) {
        if (grepl('python3', tolower(python_version))){
          python3 <- python_version
        }
      }

      if (is.null(python3)){
        base::warning("No python3 version could detected! Deconvolution with Scaden not realizable.")
        base::warning("Please run init_python(python= /path/to/python3) to set up an environment manually.")
      }
      else{
        suppressMessages(create_virtualenv(python3))
      }
    }
  }
  else{
    suppressMessages(create_virtualenv(python))
  }
}

python_available<- function(){
  return(reticulate::py_available())
}


