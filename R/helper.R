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


#' Checks wether docker/singularity are available and can be used
#' @param container The container for which the commands are tested
#' @return A boolean value
#'
check_container <- function(container = c("docker", "singularity")) {
  container.available <- (system(container, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0)

  if (!container.available) {
    message(paste0(
      "Installation of ", container, " can not be found. Please check whether you can ",
      "call 'docker' in the command line and get a help menu"
    ))
    return(FALSE)
  }

  if (container == "docker") {
    command <- "docker ps"
  } else {
    command <- "singularity instance list"
  }



  container.connectable <- (system(command, ignore.stdout = TRUE, ignore.stderr = TRUE) == 0)

  if (!container.connectable) {
    message(paste0(
      "Error during connection to ", container, ". Please check whether you can ",
      "call \'", command, "\' in the command line and get a (possibly empty) list and not an error ",
      "message"
    ))
    return(FALSE)
  }

  return(TRUE)
}



#' Setup of the singularity container
#' @param container_path the path where the singularity .sif file should be stored (optional)
#'   If the file 'fractions_latest.sif' is already present, it will be used
#'
#' @return the path to the singularity container
#'
setup_singularity_container <- function(container_path = NULL) {
  if (is.null(container_path)) {
    container_path <- file.path(path.pexpand("~"), ".local/share/omnideconv")
    dir.create(container_path, showWarnings = FALSE)
    message(paste0("singularity container written to `", container_path, "/cibersortx_fractions.sif`.
            Set the `container_path` directory to choose a different location"))
  }

  # We assume that, even in case of user provided file, the file name will
  # be 'fractions_latest.sif'
  container_file <- file.path(continer_path, "fractions_latest.sif")

  if (!file.exists(container_file)) {
    system(paste0("singularity pull --dir ", container_path, " docker://cibersortx/fractions"))
  }

  return(container_file)
}



#' Removes special characters by substituting them with unique string which should not be used naturally
#'
#' @param string The string to be escaped
#'
#' @return The String without special characters
#'
escape_special_chars <- function(string) {
  if (is.null(string)) {
    return(NULL)
  }
  string <- gsub("\u0020", "21b29fb07f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # " "
  string <- gsub("\u0021", "21b29fb17f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "!"
  string <- gsub("\u0022", "21b29fb27f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # """
  string <- gsub("\u0023", "21b29fb37f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "#"
  string <- gsub("\u0024", "21b29fb47f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "$"
  string <- gsub("\u0025", "21b29fb57f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "%"
  string <- gsub("\u0026", "21b29fb67f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "&"
  string <- gsub("\u0027", "21b29fb77f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "'"
  string <- gsub("\u0028", "21b29fb87f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "("
  string <- gsub("\u0029", "21b29fb97f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ")"
  string <- gsub("\u002A", "21b29fba7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "*"
  string <- gsub("\u002B", "21b2c6e87f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "+"
  string <- gsub("\u002C", "21b2c6e97f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ","
  string <- gsub("\u002D", "21b2c7567f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "-"
  string <- gsub("\u002E", "21b2c7577f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "."
  string <- gsub("\u002F", "21b2c7587f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "/"
  string <- gsub("\u003A", "21b2c7597f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ":"
  string <- gsub("\u003B", "21b2c75a7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ";"
  string <- gsub("\u003C", "21b2c75b7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "<"
  string <- gsub("\u003D", "21b2c75c7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "="
  string <- gsub("\u003E", "21b2c75d7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ">"
  string <- gsub("\u003F", "21b2c75e7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "?"
  string <- gsub("\u0040", "21b2c75f7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "@"
  string <- gsub("\u005B", "21b2c7607f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "["
  string <- gsub("\u005C", "21b2ee347f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "\"
  string <- gsub("\u005D", "21b2ee357f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "]"
  string <- gsub("\u005E", "21b2ee367f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "^"
  string <- gsub("\u005F", "21b2ee377f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "_"
  string <- gsub("\u0060", "21b2ee387f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "`"
  string <- gsub("\u00B4", "21b2ee397f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "´"
  string <- gsub("\u007B", "21b2ee3a7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "{"
  string <- gsub("\u007C", "21b2ee3b7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "|"
  string <- gsub("\u007D", "21b2ee3c7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "}"
  string <- gsub("\u007E", "21b2ee3d7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "~"
  string <- gsub("\u00A7", "21b2ee3e7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "§"
  string <- gsub("\u00DF", "21b315447f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "ß"

  return(string)
}


#' Removes the substitutions and turns them back into the special characters
#'
#' @param string The string to be de-escaped
#'
#' @return The String with special characters
#'
deescape_special_chars <- function(string) {
  if (is.null(string)) {
    return(NULL)
  }
  string <- gsub("21b29fb07f8711ec9bf265fb9bf6ab9c", "\u0020", string, fixed = TRUE) # " "
  string <- gsub("21b29fb17f8711ec9bf265fb9bf6ab9c", "\u0021", string, fixed = TRUE) # "!"
  string <- gsub("21b29fb27f8711ec9bf265fb9bf6ab9c", "\u0022", string, fixed = TRUE) # """
  string <- gsub("21b29fb37f8711ec9bf265fb9bf6ab9c", "\u0023", string, fixed = TRUE) # "#"
  string <- gsub("21b29fb47f8711ec9bf265fb9bf6ab9c", "\u0024", string, fixed = TRUE) # "$"
  string <- gsub("21b29fb57f8711ec9bf265fb9bf6ab9c", "\u0025", string, fixed = TRUE) # "%"
  string <- gsub("21b29fb67f8711ec9bf265fb9bf6ab9c", "\u0026", string, fixed = TRUE) # "&"
  string <- gsub("21b29fb77f8711ec9bf265fb9bf6ab9c", "\u0027", string, fixed = TRUE) # "'"
  string <- gsub("21b29fb87f8711ec9bf265fb9bf6ab9c", "\u0028", string, fixed = TRUE) # "("
  string <- gsub("21b29fb97f8711ec9bf265fb9bf6ab9c", "\u0029", string, fixed = TRUE) # ")"
  string <- gsub("21b29fba7f8711ec9bf265fb9bf6ab9c", "\u002A", string, fixed = TRUE) # "*"
  string <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "\u002B", string, fixed = TRUE) # "+"
  string <- gsub("21b2c6e97f8711ec9bf265fb9bf6ab9c", "\u002C", string, fixed = TRUE) # ","
  string <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9c", "\u002D", string, fixed = TRUE) # "-"
  string <- gsub("21b2c7577f8711ec9bf265fb9bf6ab9c", "\u002E", string, fixed = TRUE) # "."
  string <- gsub("21b2c7587f8711ec9bf265fb9bf6ab9c", "\u002F", string, fixed = TRUE) # "/"
  string <- gsub("21b2c7597f8711ec9bf265fb9bf6ab9c", "\u003A", string, fixed = TRUE) # ":"
  string <- gsub("21b2c75a7f8711ec9bf265fb9bf6ab9c", "\u003B", string, fixed = TRUE) # ";"
  string <- gsub("21b2c75b7f8711ec9bf265fb9bf6ab9c", "\u003C", string, fixed = TRUE) # "<"
  string <- gsub("21b2c75c7f8711ec9bf265fb9bf6ab9c", "\u003D", string, fixed = TRUE) # "="
  string <- gsub("21b2c75d7f8711ec9bf265fb9bf6ab9c", "\u003E", string, fixed = TRUE) # ">"
  string <- gsub("21b2c75e7f8711ec9bf265fb9bf6ab9c", "\u003F", string, fixed = TRUE) # "?"
  string <- gsub("21b2c75f7f8711ec9bf265fb9bf6ab9c", "\u0040", string, fixed = TRUE) # "@"
  string <- gsub("21b2c7607f8711ec9bf265fb9bf6ab9c", "\u005B", string, fixed = TRUE) # "["
  string <- gsub("21b2ee347f8711ec9bf265fb9bf6ab9c", "\u005C", string, fixed = TRUE) # "\"
  string <- gsub("21b2ee357f8711ec9bf265fb9bf6ab9c", "\u005D", string, fixed = TRUE) # "]"
  string <- gsub("21b2ee367f8711ec9bf265fb9bf6ab9c", "\u005E", string, fixed = TRUE) # "^"
  string <- gsub("21b2ee377f8711ec9bf265fb9bf6ab9c", "\u005F", string, fixed = TRUE) # "_"
  string <- gsub("21b2ee387f8711ec9bf265fb9bf6ab9c", "\u0060", string, fixed = TRUE) # "`"
  string <- gsub("21b2ee397f8711ec9bf265fb9bf6ab9c", "\u00B4", string, fixed = TRUE) # "´"
  string <- gsub("21b2ee3a7f8711ec9bf265fb9bf6ab9c", "\u007B", string, fixed = TRUE) # "{"
  string <- gsub("21b2ee3b7f8711ec9bf265fb9bf6ab9c", "\u007C", string, fixed = TRUE) # "|"
  string <- gsub("21b2ee3c7f8711ec9bf265fb9bf6ab9c", "\u007D", string, fixed = TRUE) # "}"
  string <- gsub("21b2ee3d7f8711ec9bf265fb9bf6ab9c", "\u007E", string, fixed = TRUE) # "~"
  string <- gsub("21b2ee3e7f8711ec9bf265fb9bf6ab9c", "\u00A7", string, fixed = TRUE) # "§"
  string <- gsub("21b315447f8711ec9bf265fb9bf6ab9c", "\u00DF", string, fixed = TRUE) # "ß"
  return(string)
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
        message("Setting python version in miniconda to be 3.9")
        Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.9)
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
