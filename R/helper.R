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
        message("Setting python version in miniconda to be 3.8")
        Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.8)
        message("Setting up miniconda environment..")
        suppressMessages(reticulate::install_miniconda())
      }
      reticulate::use_miniconda(condaenv = "r-reticulate", required = TRUE)
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

#' Normalize deconvolution result
#'
#' @param deconv_result The original deconvolution result
#'
#' @return A matrix with the rowsums of one and no negative values
#' @export
normalize_deconv_results <- function(deconv_result) {
  celltypes <- colnames(deconv_result)
  deconv_result[deconv_result < 0] <- 0
  deconv_result <- t(apply(deconv_result, 1, function(row) row / sum(row)))
  # Apply returns a vector when only supplied with one celltype. To counter it and return a matrix
  # and not a vector, this operation is needed
  if (length(celltypes) == 1) {
    deconv_result <- t(deconv_result)
    colnames(deconv_result) <- celltypes
  }
  return(deconv_result)
}

check_data <- function(single_cell_object, cell_type_annotations, bulk_gene_expression) {
  if (!is.null(single_cell_object)) {
    if ("character" %in% unique(apply(single_cell_object, 1, class))) {
      stop(
        "The single cell object matrix contains entries with the class 'character'. Please make ",
        "sure that it only contains numerics."
      )
    }
    if (!is.null(cell_type_annotations)) {
      if (ncol(single_cell_object) != length(cell_type_annotations)) {
        stop(
          "The single cell object contains ", ncol(single_cell_object), " cells while your cell ",
          "type annotations contain ", length(cell_type_annotations), " cells."
        )
      }
    }
  }
  if (!is.null(bulk_gene_expression)) {
    if ("character" %in% unique(apply(bulk_gene_expression, 1, class))) {
      stop(
        "The bulk gene expression matrix contains entries with the class 'character'. Please ",
        "make sure that it only contains numerics."
      )
    }
  }
}
