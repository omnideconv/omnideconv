#' Calculates the signature model with AutoGeneS
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes, columns are samples. Row and column names need to be set.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object
#' @param bulk_gene_expression OPTIONAL: A matrix of bulk data. Rows are genes, columns are samples. If the bulk data is supplied, the single cell object loses all gene rows not contained in the bulk data as not to create solution sets with them
#' @param ngen Number of generations. The higher, the longer it takes
#' @param mode In standard mode, the number of genes of a selection is allowed to vary arbitrarily. In fixed mode, the number of selected genes is fixed (using nfeatures)
#' @param nfeatures Number of genes to be selected in fixed mode
#' @param weights Weights applied to the objectives. For the optimization, only the sign is relevant: 1 means to maximize the respective objective, -1 to minimize it and 0 means to ignore it. The weight supplied here will be the default weight for selection. There must be as many weights as there are objectives
#' @param objectives The objectives to maximize or minimize. Must have the same length as weights. The default objectives (correlation, distance) can be referred to using strings. For custom objectives, a function has to be passed
#' @param seed Seed for random number generators
#' @param population_size Size of every generation (mu parameter)
#' @param offspring_size Number of individuals created in every generation (lambda parameter)
#' @param crossover_pb Crossover probability
#' @param mutation_pb Mutation probability
#' @param mutate_flip_pb Mutation flipping probability (fixed mode)
#' @param crossover_thres Crossover threshold (standard mode)
#' @param ind_standard_pb Probability used to generate initial population in standard mode
#' @param plot_weights Plotting: Weights with which to weight the objective values. For example, (-1,2) will minimize the first objective and maximize the the second (with higher weight)
#' @param plot_objectives Plotting: The objectives to be plotted. Contains indices of objectives. The first index refers to the objective that is plotted on the x-axis. For example, (2,1) will plot the third objective on the x-axis and the second on the y-axis
#' @param index Plotting: If one int is passed, return pareto[index] If two ints are passed, the first is an objective (0 for the first). The second is the nth element if the solutions have been sorted by the objective in ascending order. For example, (0,1) will return the solution that has the second-lowest value in the first objective. (1,-1) will return the solution with the highest value in the second objective
#' @param close_to Plotting: Select the solution whose objective value is closest to a certain value. Assumes (objective,value). For example, (0,100) will select the solution whose value for the first objective is closest to 100
#' @param plot Whether to produce a plot at all
#' @param verbose Whether to produce an output on the console
#'
#' @return The path to the pickle file needed for the deconvolution with autogenes.
#' @export
#'
#'
build_model_autogenes <- function(single_cell_object, cell_type_annotations, bulk_gene_expression = NULL,
                                  ngen=2,mode=c("standard","fixed"),nfeatures=NULL,weights=c(-1,1),objectives=c("correlation","distance"),seed=0,population_size=100,
                                  offspring_size=50,crossover_pb=0.7,mutation_pb=0.3,mutate_flip_pb=1E-3,crossover_thres=1000,ind_standard_pb=0.1,plot_weights=NULL,
                                  plot_objectives=c(0,1),index=NULL,close_to=NULL,plot=FALSE,verbose=FALSE){
  if (!is.null(bulk_gene_expression)){
    single_cell_object <- single_cell_object[intersect(rownames(single_cell_object),rownames(bulk_gene_expression)),]
  }
  if (length(mode)>1){
    mode <- mode[1]
  }

  sce <- matrix_to_singlecellexperiment(single_cell_object,cell_type_annotations)
  #TODO: Not export it as global
  ad <<- singlecellexperiment_to_anndata(sce)

  autogenes_checkload()
  if (verbose){
    base::message("import autogenes as ag")
    base::message("ag.init(r.ad,celltype_key='label')")
  }
  reticulate::py_run_string("import autogenes as ag")
  reticulate::py_run_string("ag.init(r.ad,celltype_key='label')")

  command <- glue::glue("ag.optimize(ngen={ngen},weights={create_tuple_string(weights)},objectives={create_tuple_string(objectives,strings=TRUE)},verbose={transform_boolean(verbose)},nfeatures={nfeatures},seed={seed},mode={quote_string(mode)},population_size={population_size},offspring_size={offspring_size},crossover_pb={crossover_pb},mutation_pb={mutation_pb},mutate_flip_pb={mutate_flip_pb},crossover_thres={crossover_thres},ind_standard_pb={ind_standard_pb})", .transformer = null_transformer("None"))
  if (verbose){
    base::message(command)
  }
  reticulate::py_run_string(command)
  if (plot){

    if (sum(!sapply(list(weights,index,close_to),is.null))>1){
      base::stop('When selecting a solution in autogenes only one of the parameters "plot_weights", "index" and "close_to" can be used')
    }
    all_params <- list(glue::glue(",weights={create_tuple_string(plot_weights)}", .transformer = null_transformer("None")),
                       glue::glue(",index={create_tuple_string(index)}", .transformer = null_transformer("None")),
                       glue::glue(",close_to={close_to}", .transformer = null_transformer("None")))
    criterion <- paste0(unlist(all_params[!sapply(list(weights,index,close_to),is.null)])," ")

    command <- glue::glue("ag.plot(objectives={create_tuple_string(plot_objectives)}{criterion})", .transformer = null_transformer("None"))
    if (verbose){
      base::message(command)
    }
    reticulate::py_run_string(command)
  }
  filename <- tempfile(fileext = ".pickle")
  command <- glue::glue('ag.save(r{quote_string(filename)})')
  if (verbose){
    base::message(command)
  }
  reticulate::py_run_string(command)
  return(filename)
}

#' Deconvolution Analysis using AutoGeneS
#'
#' @param bulk_gene_expression Dataframe or matrix of bulk RNA-seq data (genes x individuals)
#' @param signature Path to a .pickle file, created with the build_model method
#' @param model Regression model. Available options: NuSVR ("nusvr"), non-negative least squares("nnls") and linear model ("linear")
#' @param nu Nu parameter for NuSVR
#' @param C C parameter for NuSVR
#' @param kernel Kernel parameter for NuSVR
#' @param degree Degree parameter for NuSVR
#' @param gamma Gamma parameter for NuSVR
#' @param coef0 Coef0 parameter for NuSVR
#' @param shrinking Shrinking parameter for NuSVR
#' @param tol Tol parameter for NuSVR
#' @param cache_size Cache_size parameter for NuSVR
#' @param max_iter Max_iter parameter for NuSVR
#' @param weights Select Solution: Weights with which to weight the objective values. For example, (-1,2) will minimize the first objective and maximize the the second (with more weight)
#' @param index Select Solution: If one int is passed, return pareto[index] If two ints are passed, the first is an objective (0 for the first). The second is the nth element if the solutions have been sorted by the objective in ascending order. For example, (0,1) will return the solution that has the second-lowest value in the first objective. (1,-1) will return the solution with the highest value in the second objective
#' @param close_to Select Solution: Select the solution whose objective value is close to a certain value. Assumes (objective,value). For example, (0,100) will select the solution whose value for the first objective is closest to 100
#' @param verbose Whether the algorithm should print out what it is doing
#'
#' @return cell proportion matrix
#' @export
#'
deconvolute_autogenes <- function(bulk_gene_expression, signature, model=c("nusvr","nnls","linear"),
                                  nu=0.5,C=0.5,kernel="linear",degree=3,gamma="scale",coef0=0.0,
                                  shrinking=TRUE,tol=1E-3,cache_size=200,verbose = FALSE,max_iter=-1,
                                  weights=NULL,index=NULL,close_to=NULL){
  if(length(model)>1){
    model <- model[1]
  }

  if (!file.exists(signature)){
    base::stop("The signature parameter has to be a file path to a .pickle file, created with the build_model method. This file was not found.")
  }

  if (verbose){
    base::message(glue::glue("ag.load(r{quote_string(signature)})", .transformer = null_transformer("None")))
  }
  reticulate::py_run_string(glue::glue("ag.load(r{quote_string(signature)})", .transformer = null_transformer("None")))

  if (sum(!sapply(list(weights,index,close_to),is.null))>1){
    base::stop('When selecting a solution in autogenes only one of the parameters "weights", "index" and "close_to" can be used')
  }
  all_params <- list(glue::glue("weights={create_tuple_string(weights)}", .transformer = null_transformer("None")),
                     glue::glue("index={create_tuple_string(index)}", .transformer = null_transformer("None")),
                     glue::glue("close_to={close_to}", .transformer = null_transformer("None")))
  command = paste0("ag.select(",unlist(all_params[!sapply(list(weights,index,close_to),is.null)]),")")
  if (verbose){
    base::message(command)
  }
  reticulate::py_run_string(command)

  if (verbose){
    base::message("sel = ag.selection()")
    reticulate::py_run_string("sel = ag.selection()")
    base::message("genes = ag.adata().var_names")
    reticulate::py_run_string("genes = ag.adata().var_names")
    base::message("The following genes are used for the deconvolution:")
    #TODO: Change into base::message
    print(reticulate::py$genes[reticulate::py$sel])
  }
  #TODO: Change to not go into global environment
  bulk_data <<- t(bulk_gene_expression)
  command <- glue::glue("coef = ag.deconvolve(r.bulk_data, model={quote_string(model)}, nu={nu}, C={C},
                kernel={quote_string(kernel)}, degree={degree},gamma={quote_string(gamma)}, coef0={coef0},
                shrinking={transform_boolean(shrinking)}, tol={tol},cache_size={cache_size},
                verbose={transform_boolean(verbose)}, max_iter={max_iter})", .transformer = null_transformer("None"))
  if (verbose){
    base::message(command)
  }
  reticulate::py_run_string(command)
  res <- reticulate::py$coef
  if (verbose){
    base::message("")
    base::message("cell_types = ag.adata().obs_names")
  }
  reticulate::py_run_string("cell_types = ag.adata().obs_names")
  colnames(res) <- reticulate::py$cell_types
  rownames(res) <- rownames(bulk_data)
  return(res)
}




#' Install AutoGeneS
#'
#' Install AutoGeneS into a previously defined python environment.
#' A python environment can be defined by set_virtualenv() or set_python() command.
#' Alternatively a new environment can be created via create_virtualenv() method.
#'
install_autogenes <- function() {
  reticulate::py_install("autogenes",pip=TRUE)

}

#' Checks if scaden is installed.
#'
#' If it is available, the python module is imported.
#'
autogenes_checkload <- function(){
  if (!python_available()){
    base::message("Setting up python environment..")
    init_python()
    if (!python_available()){
      base::stop("Could not initiate miniconda python environment. Please set up manually with init_python(python=your/python/version)")
    }
  }
  if (!reticulate::py_module_available("autogenes")){
    install_autogenes()
  }
  reticulate::import("autogenes")
}

#' Creates a string that is evaluated by python as a collection of values
#'
#' @param vector_of_elements A vector of strings or numbers
#' @param strings Whether the input consists of strings that need to be quoted
#'
#' @return The formatted string
create_tuple_string <- function(vector_of_elements,strings=FALSE){
  if(is.null(vector_of_elements)){
    return("None")
  }
  if (strings){
    if(length(vector_of_elements)==1){
      return(quote_string(vector_of_elements))
    }
    return(paste0('("',paste(vector_of_elements,collapse = '","'),'")'))
  }
  if(length(vector_of_elements)==1){
    return(vector_of_elements)
  }
  return(paste0('(',paste(vector_of_elements,collapse = ','),')'))
}

#' Transforms the R TRUE and FALSE into Pythons True and False
#'
#' @param bool The original boolean
#'
#' @return The python readable boolean
transform_boolean <- function(bool){
  if(bool){
    return("True")
  }
  return("False")
}

#' Surrounds a string with quotes
#'
#' @param string The string to be quoted
#'
#' @return The quoted string
quote_string <- function(string){
  return(paste0('"',string,'"'))
}

#' Null transformer for glue: Replaces NULL values with something else
#'
#' @param str The replacement string
#'
#' @return A function used for the transformation
null_transformer <- function(str = "NULL") {
  function(text, envir) {
    out <- glue::identity_transformer(text, envir)
    if (is.null(out)) {
      return(str)
    }
    out
  }
}
