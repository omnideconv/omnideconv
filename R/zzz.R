#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#'
NULL

.onLoad <- function(libname, pkgname) {

  # We ensure to have reticulate
  if (!dir.exists(reticulate::miniconda_path())) {
    message("Setting python version in miniconda to be 3.9")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.9)
    message("Setting up miniconda environment..")
    suppressMessages(reticulate::install_miniconda())
  }


  # We ensure to have the r-reticulate env
  if (!file.exists(reticulate::conda_python("r-reticulate"))) {
    reticulate::conda_create(envname = "r-reticulate")
    reticulate::py_config()
  }

  paths <- reticulate::conda_list()
  path <- paths[paths$name == 'r-reticulate', 2]
  if(.Platform$OS.type == "windows"){
    path <- gsub("\\\\", "/", path)
  }
  path.bin <- gsub("/envs/r-reticulate/python.exe", "/library/bin", path)
  Sys.setenv(PATH= paste(path.bin,Sys.getenv()["PATH"],sep=";"))
  Sys.setenv(RETICULATE_PYTHON = path)
  #library(reticulate)


  reticulate::use_miniconda(condaenv = "r-reticulate", required = TRUE)
  reticulate::py_config()
  reticulate::configure_environment(pkgname, force = TRUE)

  if (!reticulate::py_module_available("anndata")) {
    anndata::install_anndata()
  }
}
