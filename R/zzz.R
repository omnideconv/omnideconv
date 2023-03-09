#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#'
#' @name fix_dependencies
NULL

.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname)
  reticulate::py_config()
}
