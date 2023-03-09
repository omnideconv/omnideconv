#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#'
NULL

.onLoad <- function(libname, pkgname) {
  reticulate::py_config()
  reticulate::configure_environment(pkgname, force = TRUE)
}
