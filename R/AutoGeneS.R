#' Calculates the signature model with AutoGeneS.
#' This function can be used in case that the signature should be examined. The final signature is stored as a pickle file.
#'
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param bulk_gene_expression OPTIONAL: A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set. If the bulk data is supplied, the single cell object
#'   loses all gene rows not contained in the bulk data as not to create solution sets with them.
#' @param ngen Number of generations. The higher, the longer it takes.
#' @param mode In standard mode, the number of genes of a selection is allowed to vary arbitrarily.
#'   In fixed mode, the number of selected genes is fixed (using nfeatures).
#' @param nfeatures Number of genes to be selected in fixed mode.
#' @param weights Weights applied to the objectives. For the optimization, only the sign is
#'   relevant: 1 means to maximize the respective objective, -1 to minimize it and 0 means to ignore
#'   it. The weight supplied here will be the default weight for selection. There must be as many
#'   weights as there are objectives.
#' @param objectives The objectives to maximize or minimize. Must have the same length as weights.
#'   The default objectives (correlation, distance) can be referred to using strings. For custom
#'   objectives, a function has to be passed.
#' @param seed Seed for random number generators.
#' @param population_size Size of every generation (mu parameter).
#' @param offspring_size Number of individuals created in every generation (lambda parameter).
#' @param crossover_pb Crossover probability.
#' @param mutation_pb Mutation probability.
#' @param mutate_flip_pb Mutation flipping probability (fixed mode).
#' @param crossover_thres Crossover threshold (standard mode).
#' @param ind_standard_pb Probability used to generate initial population in standard mode.
#' @param plot_weights Plotting: Weights with which to weight the objective values. For example,
#'   (-1,2) will minimize the first objective and maximize the the second (with higher weight).
#' @param plot_objectives Plotting: The objectives to be plotted. Contains indices of objectives.
#'   The first index refers to the objective that is plotted on the x-axis. For example, (2,1) will
#'   plot the third objective on the x-axis and the second on the y-axis.
#' @param index Plotting: If one int is passed, return pareto\[index\] If two ints are passed, the
#'   first is an objective (0 for the first). The second is the nth element if the solutions have
#'   been sorted by the objective in ascending order. For example, (0,1) will return the solution
#'   that has the second-lowest value in the first objective. (1,-1) will return the solution with
#'   the highest value in the second objective.
#' @param close_to Plotting: Select the solution whose objective value is closest to a certain
#'   value. Assumes (objective,value). For example, (0,100) will select the solution whose value
#'   for the first objective is closest to 100.
#' @param plot Whether to produce a plot at all. This just hands over the reticulate plot and
#'   has some visualization problems. To get a normal plot, use the pickle file, open it in python
#'   and use the plot method there.
#' @param output_dir path to directory where the picke output file will be saved. Default is tempdir().
#' @param verbose Whether to produce an output on the console.
#'
#' @return The path to the pickle file needed for the deconvolution with AutoGeneS.
#' @export
#'
#'
build_model_autogenes <- function(single_cell_object, cell_type_annotations,
                                  bulk_gene_expression = NULL, ngen = 5000,
                                  mode = c("fixed", "standard"), nfeatures = 500,
                                  weights = list(-1, 1), objectives = list("correlation", "distance"),
                                  seed = 0, population_size = 100, offspring_size = 100,
                                  crossover_pb = 0.7, mutation_pb = 0.3, mutate_flip_pb = 1E-3,
                                  crossover_thres = 1000, ind_standard_pb = 0.1,
                                  plot_weights = NULL, plot_objectives = c(0, 1), index = NULL,
                                  close_to = NULL, plot = FALSE, output_dir = tempdir(),
                                  verbose = FALSE) {

  message(
    "The deconvolution with AutoGeneS is done in only two steps, however due to runtime reasons we chose to implement
    it as a one-step method. Please just use the deconvolute method.
    If you only want to calculate the signature, without deconvolution, please use this function."
  )

  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (!is.null(bulk_gene_expression)) {
    single_cell_object <- single_cell_object[intersect(
      rownames(single_cell_object),
      rownames(bulk_gene_expression)
    ), ]
  }
  if (length(mode) > 1) {
    mode <- mode[1]
    message(paste0(mode, " was chosen because multiple values were supplied for \"mode\""))
  }

  # anndata_checkload()

  sce <- matrix_to_singlecellexperiment(single_cell_object, cell_type_annotations)
  ad <- singlecellexperiment_to_anndata(sce)

  # autogenes_checkload()
  ag <- reticulate::import("autogenes")
  ag$init(ad, celltype_key = "label")

  params <- list(
    ngen = as.integer(ngen), weights = weights, objectives = objectives,
    verbose = verbose, seed = as.integer(seed), mode = mode,
    population_size = as.integer(population_size),
    offspring_size = as.integer(offspring_size), crossover_pb = crossover_pb,
    mutation_pb = mutation_pb
  )

  if (mode == "standard") {
    params <- c(params,
      crossover_thres = as.integer(crossover_thres),
      ind_standard_pb = ind_standard_pb
    )
  } else if (mode == "fixed") {
    params <- c(params, nfeatures = as.integer(nfeatures), mutate_flip_pb = mutate_flip_pb)
  } else {
    stop(paste0("Mode ", mode, " not recognized. Please try 'standard' or 'fixed'"))
  }

  do.call(ag$optimize, params)

  if (plot) {
    if (sum(!sapply(list(weights, index, close_to), is.null)) > 1) {
      stop(
        "When selecting a solution in autogenes only one of the parameters 'plot_weights', ",
        "'index' and 'close_to' can be used. It is also possible to use none of them"
      )
    }
    plot_objectives <- as.integer(plot_objectives)
    if (!is.null(plot_weights)) {
      ag$plot(objectives = plot_objectives, weights = as.integer(plot_weights))
    } else if (!is.null(index)) {
      ag$plot(objectives = plot_objectives, index = as.integer(index))
    } else if (!is.null(close_to)) {
      ag$plot(objectives = plot_objectives, close_to = as.integer(close_to))
    } else {
      ag$plot(objectives = plot_objectives)
    }
  }
  filename <- tempfile(fileext = ".pickle", tmpdir = output_dir)
  ag$save(filename)
  if (verbose) {
    message("Successfully saved")
  }
  return(filename)
}

#' Extraction of the signature matrix from autogenes .pickle files
#'
#' @param autogenes_pickle_path The path of the .pickle file generated
#'    by autogenes, exactly as returned by the `build_model_autogenes` function
#' @param single_cell_object The matrix with the single-cell data used to build the signature.
#'   Rows are genes, columns are samples. Row and column names need to be set.
#' @param cell_type_annotations The vector of the cell type annotations used to
#'   build the signature. Has to be in the same order as the samples in single_cell_object.
#'
#' @return The signature matrix.
#' @export
#'
#'
extract_signature_autogenes <- function(autogenes_pickle_path,
                                        single_cell_object, cell_type_annotations) {
  sce <- matrix_to_singlecellexperiment(
    as.matrix(single_cell_object),
    cell_type_annotations
  )
  ad <- singlecellexperiment_to_anndata(sce)

  ag <- reticulate::import("autogenes")
  ag$init(ad, celltype_key = "label")

  expr.median <- ag$init(ad, celltype_key = "label")

  median.gene.expr <- expr.median$X

  ag$load(autogenes_pickle_path)
  genes <- ag$adata()$var_names
  selection <- ag$select()
  genes.used <- genes[selection]

  signature <- median.gene.expr[, genes.used]
  signature <- t(signature)
  return(signature)
}


#' Deconvolution Analysis using AutoGeneS.
#' One-step function that performs signature building and deconvolution in one step without saving the signature in between.
#' A signature that has been created with @seealso [build_model_autogenes()] can be supplied as input to this function with the
#' signature parameter, although it is not mandatory.
#'
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param signature OPTIONAL: Path to a .pickle file, created with the build_model method of AutoGeneS.
#' @param ngen (SIGNATURE) Number of generations. The higher, the longer it takes.
#' @param mode (SIGNATURE) In standard mode, the number of genes of a selection is allowed to vary arbitrarily.
#'   In fixed mode, the number of selected genes is fixed (using nfeatures).
#' @param nfeatures (SIGNATURE) Number of genes to be selected in fixed mode.
#' @param weights_signature (SIGNATURE) Weights applied to the objectives. For the optimization, only the sign is
#'   relevant: 1 means to maximize the respective objective, -1 to minimize it and 0 means to ignore
#'   it. The weight supplied here will be the default weight for selection. There must be as many
#'   weights as there are objectives.
#' @param objectives (SIGNATURE) The objectives to maximize or minimize. Must have the same length as weights.
#'   The default objectives (correlation, distance) can be referred to using strings. For custom
#'   objectives, a function has to be passed.
#' @param seed (SIGNATURE) Seed for random number generators.
#' @param population_size (SIGNATURE) Size of every generation (mu parameter).
#' @param offspring_size (SIGNATURE) Number of individuals created in every generation (lambda parameter).
#' @param crossover_pb (SIGNATURE) Crossover probability.
#' @param mutation_pb (SIGNATURE) Mutation probability.
#' @param mutate_flip_pb (SIGNATURE) Mutation flipping probability (fixed mode).
#' @param crossover_thres (SIGNATURE) Crossover threshold (standard mode).
#' @param ind_standard_pb (SIGNATURE) Probability used to generate initial population in standard mode.
#' @param plot_weights (SIGNATURE) Plotting: Weights with which to weight the objective values. For example,
#'   (-1,2) will minimize the first objective and maximize the the second (with higher weight).
#' @param plot_objectives (SIGNATURE) Plotting: The objectives to be plotted. Contains indices of objectives.
#'   The first index refers to the objective that is plotted on the x-axis. For example, (2,1) will
#'   plot the third objective on the x-axis and the second on the y-axis.
#' @param plot (SIGNATURE) Whether to produce a plot at all. This just hands over the reticulate plot and
#'   has some visualization problems. To get a normal plot, use the pickle file, open it in python
#'   and use the plot method there.
#' @param output_dir (SIGNATURE) path to directory where the picke output file will be saved. Default is tempdir().
#' @param model (DECONVOLUTION) Regression model. Available options: NuSVR ("nusvr"), non-negative least
#'   squares("nnls") and linear model ("linear").
#' @param nu (DECONVOLUTION) Nu parameter for NuSVR.
#' @param C (DECONVOLUTION) C parameter for NuSVR.
#' @param normalize_results (DECONVOLUTION) wether to normalize results according to the regression model used.
#'   Default is TRUE
#' @param kernel (DECONVOLUTION) Kernel parameter for NuSVR.
#' @param degree (DECONVOLUTION) Degree parameter for NuSVR.
#' @param gamma (DECONVOLUTION) Gamma parameter for NuSVR.
#' @param coef0 (DECONVOLUTION) Coef0 parameter for NuSVR.
#' @param shrinking (DECONVOLUTION) Shrinking parameter for NuSVR.
#' @param tol (DECONVOLUTION) Tol parameter for NuSVR.
#' @param cache_size (DECONVOLUTION) Cache_size parameter for NuSVR.
#' @param max_iter (DECONVOLUTION) Max_iter parameter for NuSVR.
#' @param weights_deconvolution (DECONVOLUTION) Select Solution: Weights with which to weight the objective values. For example,
#'   (-1,1) will minimize the first objective and maximize the the second (with more weight).
#' @param close_to Select Solution: Select the solution whose objective value is close to a certain
#'   value. Assumes (objective,value). For example, (0,100) will select the solution whose value
#'   for the first objective is closest to 100.
#' @param index Select Solution: If one int is passed, return pareto\[index\] If two ints are passed,
#'   the first is an objective (0 for the first). The second is the nth element if the solutions
#'   have been sorted by the objective in ascending order. For example, (0,1) will return the
#'   solution that has the second-lowest value in the first objective. (1,-1) will return the
#'   solution with the highest value in the second objective.
#' @param verbose Whether to produce an output on the console.
#'
#' @return A list with two elements: 'proportions' is the matrix of cell proportions and
#'   'genes_used' is a vector containing the names of the genes used for the deconvolution, what
#'   is called "solution" by AutoGeneS.
#' @export
#'
deconvolute_autogenes <- function(single_cell_object, bulk_gene_expression, cell_type_annotations,
                                  signature = NULL, ngen = 5000,
                                  mode = c("fixed", "standard"), nfeatures = 500,
                                  weights_signature = list(-1, 1), objectives = list("correlation", "distance"),
                                  seed = 0, population_size = 100, offspring_size = 100,
                                  crossover_pb = 0.7, mutation_pb = 0.3, mutate_flip_pb = 1E-3,
                                  crossover_thres = 1000, ind_standard_pb = 0.1,
                                  plot_weights = NULL, plot_objectives = c(0, 1),
                                   plot = FALSE, output_dir = tempdir(),
                                  model = c("nusvr", "nnls", "linear"), nu = 0.5, C = 0.5,
                                  normalize_results = TRUE, kernel = "linear", degree = 3,
                                  gamma = "scale", coef0 = 0.0, shrinking = TRUE,
                                  tol = 1E-3, cache_size = 200, max_iter = -1,
                                  weights_deconvolution = list(1, 0), index = NULL, close_to = NULL,
                                  verbose = FALSE) {

  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  if (!is.null(signature)) {
    if ("matrix" %in% class(signature)) {
      stop(
        "Parameter 'signature' requires the .pickle file created in the AutoGeneS build model ",
        "method, not a matrix of values."
      )
    }
    if (!file.exists(signature)) {
      stop(
        "The signature parameter has to be a file path to a .pickle file, created with the ",
        "build_model method. This file was not found."
      )
    }
  }
  if (length(mode) > 1) {
    mode <- mode[1]
    message(paste0(mode, " was chosen because multiple values were supplied for \"mode\""))
  }

  if (length(model) > 1) {
    model <- model[1]
    message(paste0(model, " was chosen because multiple values were supplied for \"model\""))
  }

  if (sum(!sapply(list(weights_deconvolution, index, close_to), is.null)) > 1) {
    stop(
      "When selecting a solution in autogenes only one of the parameters 'weights_deconvolution', ",
      "'index' and 'close_to' can be used"
    )
  }

  sce <- matrix_to_singlecellexperiment(single_cell_object, cell_type_annotations)
  ad <- singlecellexperiment_to_anndata(sce)

  ## Start building the signature ##

  ag <- reticulate::import("autogenes")
  ag$init(ad, celltype_key = "label")

  params <- list(
    ngen = as.integer(ngen), weights = weights_signature, objectives = objectives,
    verbose = verbose, seed = as.integer(seed), mode = mode,
    population_size = as.integer(population_size),
    offspring_size = as.integer(offspring_size), crossover_pb = crossover_pb,
    mutation_pb = mutation_pb
  )

  if (mode == "standard") {
    params <- c(params,
                crossover_thres = as.integer(crossover_thres),
                ind_standard_pb = ind_standard_pb
    )
  } else if (mode == "fixed") {
    params <- c(params, nfeatures = as.integer(nfeatures), mutate_flip_pb = mutate_flip_pb)
  } else {
    stop(paste0("Mode ", mode, " not recognized. Please try 'standard' or 'fixed'"))
  }

  do.call(ag$optimize, params)

  if (plot) {
    if (sum(!sapply(list(plot_weights, index, close_to), is.null)) > 1) {
      stop(
        "When selecting a solution in autogenes only one of the parameters 'plot_weights', ",
        "'index' and 'close_to' can be used. It is also possible to use none of them"
      )
    }
    plot_objectives <- as.integer(plot_objectives)
    if (!is.null(plot_weights)) {
      ag$plot(objectives = plot_objectives, weights = as.integer(plot_weights))
    } else if (!is.null(index)) {
      ag$plot(objectives = plot_objectives, index = as.integer(index))
    } else if (!is.null(close_to)) {
      ag$plot(objectives = plot_objectives, close_to = as.integer(close_to))
    } else {
      ag$plot(objectives = plot_objectives)
    }
  }

  ## Start deconvolution ##

  if (!is.null(weights_deconvolution)) {
    ag$select(weights = weights_deconvolution)
  } else if (!is.null(index)) {
    ag$select(index = index)
  } else if (!is.null(close_to)) {
    ag$select(close_to = close_to)
  } else {
    ag$select()
  }

  selection <- ag$selection()
  genes <- ag$adata()$var_names
  genes_used <- genes[selection]
  if (verbose) {
    message(length(genes_used), " genes are used for the deconvolution")
  }

  bulk_data <- data.frame(t(bulk_gene_expression))
  result <- ag$deconvolve(bulk_data,
    model = model, nu = nu, C = C,
    kernel = kernel, degree = as.integer(degree), gamma = gamma, coef0 = coef0,
    shrinking = shrinking, tol = tol, cache_size = as.integer(cache_size),
    verbose = verbose, max_iter = as.integer(max_iter)
  )

  # Addition: this is needed to adjust the format of column names
  names <- ag$adata()$obs_names

  if (is.atomic(names)) {
    colnames(result) <- names
  } else {
    col_names <- names$to_frame()
    col_names_vectorized <- as.vector(col_names[, 1])
    colnames(result) <- col_names_vectorized
  }

  rownames(result) <- rownames(bulk_data)

  if (normalize_results) {
    celltypes <- colnames(result)

    if (model == "nusvr") {
      result[result < 0] <- 0
    }

    result <- t(apply(result, 1, function(row) row / sum(row)))
    if (length(celltypes) == 1) {
      result <- t(result)
      colnames(result) <- celltypes
    }
  }

  return(list(proportions = result, genes_used = genes_used))
}




#' Install AutoGeneS
#'
#' Install AutoGeneS into a previously defined python environment.
#' A python environment can be defined by set_virtualenv() or set_python() command.
#' Alternatively a new environment can be created via create_virtualenv() method.
#'
install_autogenes <- function() {
  reticulate::py_install("autogenes", pip = TRUE)
}

#' Checks if python and the autogenes module are available and installs them if they are not.
#' @param python the python env
#'
autogenes_checkload <- function(python = NULL) {
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
  anndata_checkload()
  if (!reticulate::py_module_available("autogenes")) {
    install_autogenes()
  }
}
