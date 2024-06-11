#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#'
NULL

.onLoad <- function(libname, pkgname) {
  cli::cli_alert("checking omnideconv environment and dependencies")

  # We ensure to have reticulate
  if (!dir.exists(reticulate::miniconda_path())) {
    message("Setting python version in miniconda to be 3.8")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.8)
    message("Setting up miniconda environment..")
    suppressMessages(reticulate::install_miniconda())
  }


  # We ensure to have the r-reticulate env
  # if (!file.exists(reticulate::conda_python("r-reticulate"))) {
  if (!("r-omnideconv" %in% reticulate::conda_list()$name)) {
    reticulate::conda_create(envname = "r-omnideconv")
  }

  paths <- reticulate::conda_list()
  path <- paths[paths$name == "r-omnideconv", 2]
  if (.Platform$OS.type == "windows") {
    path <- gsub("\\\\", "/", path)
  }
  path.bin <- gsub("/envs/omnideconv/python.exe", "/library/bin", path)
  Sys.setenv(PATH = paste(path.bin, Sys.getenv()["PATH"], sep = ";"))
  Sys.setenv(RETICULATE_PYTHON = path)


  reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
  reticulate::py_config()
  reticulate::configure_environment(pkgname, force = TRUE)
  reticulate::py_install("pip==23.1")

  if (!reticulate::py_module_available("anndata")) {
    anndata::install_anndata()
    # reticulate::py_install("anndata", pip = TRUE)
  }
}
