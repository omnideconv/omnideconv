# global reference to Scaden,scanpy (will be initialized in .onLoad)
scaden <- NULL
scanpy <- NULL



.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to Scaden.
  scaden <<- reticulate::import("scaden", delay_load = TRUE)
  scanpy <<- reticulate::import("scanpy", delay_load = TRUE)
}
