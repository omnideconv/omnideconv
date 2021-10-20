#' A wrapper function whether to suppress messages
#'
#' @param verbose Whether to produce an output on the console.
#'
#' @return A function which will suppress messages or not, depending on the verbose parameter
#'
verbose_wrapper <- function(verbose) {
  return(function(method) {
    if (!verbose) {
      suppressMessages(method)
    } else {
      method
    }
  })
}


#' Docker availability check
#'
#' @return A boolean value whether docker is available on the system
#'
docker_available <- function() {
  return(system("docker", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0)
}

#' Docker connectability check
#'
#' @return A boolean value whether it is possible to connect to docker
#'
docker_connectable <- function() {
  return(system("docker ps", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0)
}

#' Removes blanks by substituting them with "._._._" which should not be used naturally
#'
#' @param string The string to be escaped
#'
#' @return The String without blanks
#'
escape_blanks <- function(string) {
  if (is.null(string)) {
    return(NULL)
  }
  return(gsub(" ", "._._._", string))
}

#' Removes the substitutions "._._._" and turns them back into blanks
#'
#' @param string The string to be de-escaped
#'
#' @return The String with blanks
#'
deescape_blanks <- function(string) {
  if (is.null(string)) {
    return(NULL)
  }
  return(gsub("._._._", " ", string))
}


### Python helper methods

#' Defines python path.
#'
#' Python3 is needed for scanpy and scaden to work.
#'
#' @param path_to_python_binaries Path to python binaries
#'
set_python <- function(path_to_python_binaries) {
  reticulate::use_python(python = path_to_python_binaries)
}


#' Initiates python environment
#'
#' @param python (optional) If own python should be used please indicate it's binaries
#'
init_python <- function(python = NULL) {
  if (!reticulate::py_available()) {
    if (is.null(python)) {
      if (!dir.exists(reticulate::miniconda_path())) {
        message("Setting up miniconda environment..")
        suppressMessages(reticulate::install_miniconda())
      }
      # ensure using a current Python version
      retuculate::conda_create(envname='omnideconv', python_version='3.8', conda=reticulate::miniconda_path())
      reticulate::use_condaenv(condaenv = "omnideconv", required = TRUE, conda=reticulate::miniconda_path())
      config <- reticulate::py_config()
      if (!python_available()) {
        message("Python not available")
        print(config)
        message(
          "Please indicate your version of python calling init_python(python=your/python)"
        )
      }
    } else {
      reticulate::use_python(python = python)
      reticulate::py_config()
    }
  }
}

#' Checks if python is available in environment
#'
#' @return boolean
#'
python_available <- function() {
  return(reticulate::py_available())
}
