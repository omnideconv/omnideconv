#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name omnideconvstartup
NULL

.onLoad <- function(libname, pkgname) {
  #no installatins have to be checked when running omnideconv inside the Docker container
  message('Skipping python installation checks ...')

}
