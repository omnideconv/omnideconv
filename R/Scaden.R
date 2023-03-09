#' Builds Scaden Model from scRNA data
#'
#' The model is saved in a defined directory.
#'
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param model_path Path where model directory should be created (optional).
#' @param temp_dir The temporary directory to use for the computations (optional)
#' @param samples Bulk simulation: Number of samples to simulate (default: 1000).
#' @param cells Bulk simulation: Number of cells per sample (default: 100).
#' @param dataset_name Bulk simulation: Name of simulated dataset (default scaden).
#' @param var_cutoff Training data processing: Filter out genes with a variance less than the
#'   specified cutoff. A low cutoff is recommended,this should only remove genes that are obviously
#'   uninformative. (default NULL).
#' @param batch_size Training  of model: Batch size to use for training (default: 128).
#' @param learning_rate Training of model: Learning rate used for training (default: 1E-4).
#' @param steps Training of model: Number of training steps (default: 5000).
#' @param verbose Whether to produce an output on the console (default: false).
#'
#' @return The path to the scaden model.
#' @export
#'
build_model_scaden <- function(single_cell_object, cell_type_annotations, bulk_gene_expression,
                               model_path = NULL, temp_dir = NULL, batch_size = 128, learning_rate = 1E-4,
                               steps = 5000, var_cutoff = NULL, cells = 100, samples = 1000,
                               dataset_name = "scaden", verbose = FALSE) {
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }


  if (ncol(bulk_gene_expression) < 2) {
    stop("Scaden requires at least two bulk samples.")
  }

  if (!verbose) {
    if (Sys.info()["sysname"] == "Windows") {
      message("The windows implementation requires verbose mode. It is now switched on.")
      verbose <- TRUE
    }
  }

  # check if Scaden is installed and loaded.
  #scaden_checkload()
  reticulate::import("scaden")
  single_cell_object <- t(single_cell_object)
  if (nrow(single_cell_object) != length(cell_type_annotations)) {
    stop(
      "Celltype labels must be same length as samples (number of columns) ",
      "in single_cell_object!"
    )
  }
  gene_labels <- colnames(single_cell_object)

  training_h5ad <- scaden_simulate(
    cell_type_annotations = cell_type_annotations, gene_labels = gene_labels,
    single_cell_object = single_cell_object, temp_dir = temp_dir, cells = cells,
    samples = samples, dataset_name = dataset_name,
    verbose = verbose
  )
  processed <- scaden_process(training_h5ad,
    temp_dir = temp_dir, bulk_gene_expression,
    verbose = verbose,
    var_cutoff = var_cutoff
  )
  output_model_path <- scaden_train(processed,
    temp_dir = temp_dir,
    model_path = model_path, verbose = verbose,
    batch_size = batch_size, learning_rate = learning_rate,
    steps = steps
  )
  return(output_model_path)
}


#' Performs deconvolution with Scaden
#'
#' @param signature Path to the model directory
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param verbose Whether to produce an output on the console (default: false).
#' @param temp_dir The temporary directory to use for the computations (optional)
#' @return A matrix of cell type proportion estimates with cell types as rows and
#'   individuals as columns.
#' @export
#'
deconvolute_scaden <- function(signature, bulk_gene_expression, temp_dir = NULL, verbose = FALSE) {
  if (is.null(signature)) {
    stop(
      "Parameter 'signature' is missing or null, but it is required. The path to the model ",
      "directory needs to be specified as the parameter signature for the deconvolute function."
    )
  }
  if ("matrix" %in% class(signature)) {
    stop(
      "Parameter 'signature' requires the path to the model directory created in the Scaden ",
      "build model method, not a matrix of values."
    )
  }
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  if (!verbose) {
    if (Sys.info()["sysname"] == "Windows") {
      message("The windows implementation requires verbose mode. It is now switched on.")
      verbose <- TRUE
    }
  }

  #scaden_checkload()
  reticulate::import("scaden")
  prediction <- scaden_predict(signature, bulk_gene_expression, temp_dir = temp_dir, verbose = verbose)

  return(t(prediction))
}

#' Install Scaden
#'
#' Install Scaden into a previously defined python environment.
#' A python environment can be defined by set_virtualenv() or set_python() command.
#' Alternatively a new environment can be created via create_virtualenv() method.
#'
install_scaden <- function() {
  reticulate::py_install("scaden", pip = TRUE)
}


#' Trains a Scaden model
#'
#' Trains a Scaden model on pre-processed training data.
#'
#' @param h5ad_processed pre-processed training data
#' @param temp_dir The temporary directory to use for the computations (optional)
#' @param batch_size Batch size to use for training. Default: 128
#' @param learning_rate Learning rate used for training. Default: 0.0001
#' @param model_path Path to store the model in
#' @param steps Number of training steps
#' @param verbose Whether to produce an output on the console. (default: false)
#'
#'
#' @return Scaden model
#'
scaden_train <- function(h5ad_processed, temp_dir = NULL, batch_size = 128, learning_rate = 0.0001,
                         model_path = NULL, steps = 5000, verbose = FALSE) {
  if (verbose) message("Training model")

  # create temporary directory where Scaden input files should be saved at.
  tmp_dir <- temp_dir
  if (is.null(temp_dir)) {
    tmp_dir <- tempdir()
    dir.create(tmp_dir, showWarnings = FALSE)
  }


  # the file is only created temporarily, later deleted at unlink(tmp)
  h5ad_processed_tmp <- tempfile(tmpdir = tmp_dir)
  write_anndata(h5ad_processed, h5ad_processed_tmp)

  if (is.null(model_path)) {
    currentwd <- getwd()
    model_path <- tryCatch(
      {
        setwd(tmp_dir)
        model_path <- paste0(tmp_dir, "/model")
        if (dir.exists(model_path)) {
          unlink(model_path, recursive = TRUE)
        }
        dir.create(model_path, showWarnings = FALSE)
        model_path
      },
      error = function(cond) {
        setwd(currentwd)
        stop(cond)
      },
      warning = function(cond) {
        warning(cond)
      },
      finally = {
        setwd(currentwd)
      }
    )
  }

  # Calling Scaden command
  system(paste(
    "scaden train", h5ad_processed_tmp, "--batch_size", batch_size, "--learning_rate",
    learning_rate, "--steps", steps, "--model_dir", model_path
  ),
  ignore.stdout = !verbose, ignore.stderr = !verbose
  )

  model_dirs <- list.dirs(path = model_path, full.names = FALSE, recursive = FALSE)
  model_names <- c("m1024", "m256", "m512")
  for (model in model_names) {
    if (!(model %in% model_dirs)) {
      stop("Model generation failed!")
    }
  }

  return(model_path)
}



#' Pre-processes training data
#'
#' Training data needs to be processed before a model can be trained.
#'
#' @param h5ad File that should be processed. Must be in AnnData format (.h5ad)
#' @param temp_dir The temporary directory to use for the computations (optional)
#' @param bulk_gene_expression Bulk RNA-seq data. (genes x individuals)
#' @param var_cutoff Filter out genes with a variance less than the specified cutoff. A low cutoff
#'   is recommended,this should only remove genes that are obviously uninformative.
#' @param verbose Whether to produce an output on the console. (default: false)
#'
#' @return processed training file. (.h5ad format)
#'
scaden_process <- function(h5ad, temp_dir = NULL, bulk_gene_expression, var_cutoff = NULL, verbose = FALSE) {
  if (verbose) message("Processing training data for model creation ...")

  # create temporary directory where Scaden input files should be saved at.
  tmp_dir <- temp_dir
  if (is.null(temp_dir)) {
    tmp_dir <- tempdir()
    dir.create(tmp_dir, showWarnings = FALSE)
  }

  out <- tryCatch(
    {
      dir.create(tmp_dir, showWarnings = FALSE)

      h5ad_tmp <- tempfile(tmpdir = tmp_dir, fileext = ".h5ad")
      write_anndata(h5ad, h5ad_tmp)

      bulk_gene_expression_tmp <- tempfile(tmpdir = tmp_dir)
      utils::write.table(bulk_gene_expression,
        file = bulk_gene_expression_tmp, sep = "\t", row.names = TRUE,
        col.names = NA, quote = FALSE
      )

      processed_h5ad <- tempfile(fileext = ".h5ad", tmpdir = tmp_dir)


      if (is.null(var_cutoff)) {
        system(paste("scaden process", h5ad_tmp, bulk_gene_expression_tmp, "--processed_path", processed_h5ad),
          ignore.stdout = !verbose, ignore.stderr = !verbose
        )
      } else {
        system(paste(
          "scaden process", h5ad_tmp, bulk_gene_expression_tmp, "--processed_path", processed_h5ad,
          "--var_cutoff", var_cutoff
        ),
        ignore.stdout = !verbose, ignore.stderr = !verbose
        )
      }
      read_anndata(processed_h5ad)
    },
    error = function(cond) {
      message(
        "Error preprocessing training data! Make sure training data is not in ",
        "logarithmic space!"
      )
      message(cond)
    },
    warning = function(cond) {
      if (verbose) {
        message(cond)
      }
    }
  )
  return(out)
}

#' Predict cell proportions
#'
#' Predicts cell proportions in bulk RNA sample.
#'
#' @param model_dir Directory where model is saved
#' @param bulk_gene_expression Bulk RNA-seq data. (genes x bulk_RNA samples)
#' @param temp_dir The temporary directory to use for the computations (optional)
#' @param verbose Whether to produce an output on the console. (default: false)
#'
#' @return Cell type fractions per sample
#'
scaden_predict <- function(model_dir, bulk_gene_expression, temp_dir = NULL, verbose = FALSE) {
  if (verbose) message("Predicting cell type proportions")

  current_wd <- getwd()

  predictions <- tryCatch(
    {
      # create temporary directory where Scaden input files should be saved at.
      tmp_dir <- temp_dir
      if (is.null(temp_dir)) {
        tmp_dir <- tempdir()
        dir.create(tmp_dir, showWarnings = FALSE)
      }
      setwd(tmp_dir)

      bulk_gene_expression_tmp <- tempfile(tmpdir = tmp_dir)
      utils::write.table(bulk_gene_expression,
        file = bulk_gene_expression_tmp, sep = "\t", row.names = TRUE,
        col.names = NA, quote = FALSE
      )

      system(paste("scaden predict --model_dir", model_dir, bulk_gene_expression_tmp),
        ignore.stdout = !verbose, ignore.stderr = !verbose
      )

      t(utils::read.table(paste0(tmp_dir, "/scaden_predictions.txt"),
        sep = "\t", header = TRUE,
        row.names = 1
      ))
    },
    error = function(cond) {
      stop(cond)
    },
    warning = function(cond) {
      warning(cond)
    },
    finally = {
      setwd(current_wd)
    }
  )

  return(predictions)
}

#' Simulates example Data provided by Scaden
#'
#' Used for testing.
#'
#' @param example_data_path Path to where example data should be saved. (directory)
#' @param verbose Whether to produce an output on the console. (default: false)
#' @return List with list$simulated_h5ad =  example training data
#' and list$bulk = example bulk data.
#'
scaden_simulate_example <- function(example_data_path = NULL, verbose = FALSE) {
  current_wd <- getwd()

  if (is.null(example_data_path)) {
    tmp_dir <- tempdir()
    dir.create(tmp_dir, showWarnings = FALSE)
    setwd(tmp_dir)
  } else {
    setwd(example_data_path)
  }


  if (!verbose) {
    logfile <- tempfile(tmpdir = tmp_dir, fileext = "train.log")
    sink(file = logfile)
  }


  system("mkdir example_data")
  system("scaden example --out example_data/")
  system(paste0("scaden simulate --data ", tmp_dir, "/example_data/ -n 100 --pattern *_counts.txt"))

  simulated_h5ad <- read_anndata(paste0(tmp_dir, "/data.h5ad"))
  bulk <- utils::read.table(paste0(tmp_dir, "/example_data/example_bulk_gene_expression.txt"))

  setwd(current_wd)
  unlink(tmp_dir)

  output <- list("simulated_h5ad" = simulated_h5ad, "bulk" = bulk)
  return(output)
}

#' Simulates training data from scRNA data
#'
#' @param cell_type_annotations Vector of celltype labels. Order corresponds to rows in
#'   single_cell_object matrix.
#' @param gene_labels Vector of gene labels. Order corresponds to columns in
#'   single_cell_object matrix.
#' @param single_cell_object Matrix or dataframe of scRNA data. Rows=cells and columns=genes
#' @param temp_dir The temporary directory to use for the computations (optional)
#' @param samples Bulk simulation: Number of samples to simulate (default: 1000)
#' @param cells Bulk simulation: Number of cells per sample (default: 100)
#' @param temp_dir The temporary directory to use for the computations (optional)
#' @param dataset_name Name of dataset
#' @param verbose Whether to produce an output on the console. (default: false)
#'
#' @return Simulated bulk data of known cell type fractions
#'
scaden_simulate <- function(cell_type_annotations, gene_labels, single_cell_object, temp_dir = NULL, cells = 100,
                            samples = 1000, dataset_name = "scaden", verbose = FALSE) {
  if (verbose) {
    message(
      "Simulating training data from single cell experiment: ", samples,
      " samples of ", cells, " cells"
    )
  }
  current_wd <- getwd()


  output <- tryCatch(
    {
      # create temporary directory where Scaden input files should be saved at.
      tmp_dir <- temp_dir
      if (is.null(temp_dir)) {
        tmp_dir <- tempdir()
        dir.create(tmp_dir, showWarnings = FALSE)
      }
      setwd(tmp_dir)
      if (dir.exists(dataset_name)) {
        unlink(dataset_name, recursive = TRUE)
      }
      dir.create(dataset_name, showWarnings = FALSE)
      setwd(paste0(tmp_dir, "/", dataset_name))

      colnames(single_cell_object) <- gene_labels
      rownames(single_cell_object) <- 0:(length(rownames(single_cell_object)) - 1)
      cell_types <- data.frame("Celltype" = cell_type_annotations)

      utils::write.table(format(single_cell_object, digits = 1),
        paste0(tmp_dir, "/", dataset_name, "/", dataset_name, "_counts.txt"),
        sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE
      )
      utils::write.table(cell_types, paste0(
        tmp_dir, "/", dataset_name, "/", dataset_name,
        "_celltypes.txt"
      ),
      quote = FALSE, row.names = FALSE, col.names = TRUE
      )
      setwd(tmp_dir)

      system(paste(
        "scaden simulate --data", paste0(tmp_dir, "/", dataset_name), "-n", samples,
        "-c", cells, "--pattern *_counts.txt"
      ),
      ignore.stdout = !verbose, ignore.stderr = !verbose
      )

      # Workaround to not have any Inf values in the simulated data
      temp_output <- read_anndata(paste0(tmp_dir, "/data.h5ad"))
      value_to_set_infinities_to <- max(temp_output$X[is.finite(temp_output$X)])
      number_of_infs <- sum(temp_output$X == Inf)
      temp_output$X[temp_output$X == Inf] <- value_to_set_infinities_to * 2
      if (verbose & number_of_infs > 0) {
        message(paste0(
          number_of_infs, " Inf values were replaced by twice the maximum value (",
          value_to_set_infinities_to, "*2)"
        ))
      }
      temp_output
    },
    error = function(cond) {
      stop(cond)
    },
    warning = function(cond) {
      warning(cond)
    },
    finally = {
      setwd(current_wd)
    }
  )
  return(output)
}

#' Checks if scaden is installed.
#'
#' If it is available, the python module is imported.
#'
scaden_checkload <- function(python = NULL) {
  if (!python_available()) {
    message("Setting up python environment..")
    init_python(python)
    if (!python_available()) {
      stop(
        "Could not initiate miniconda python environment. Please set up manually with ",
        "init_python(python=your/python/version)"
      )
    }
  }
  if (!reticulate::py_module_available("scaden")) {
    install_scaden()
  }
  reticulate::import("scaden")
}
