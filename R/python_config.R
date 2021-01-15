#' Defines a virtual environment that should be used.
#'
#'
#' This environment needs to contain python>=3.7 for scanpy and scaden to work.
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
#' Python>=3.7 is needed for scanpy and scaden to work.
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
#' Python>=3.7 required for scaden and scanpy to work.
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
  reticulate::use_virtualenv("~/.virtualenvs/r-reticulate",required = T)
  base::message("Uprgading pip in environment:")
  system("pip install -U pip")
}


