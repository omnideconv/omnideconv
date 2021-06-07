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
#' @return The path to the pickle file needed for the deconvolution with autogenes
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
  ad <- singlecellexperiment_to_anndata(sce)

  autogenes_checkload()
  ag <- reticulate::import("autogenes")
  ag$init(ad,celltype_key="label")

  if (mode=="standard"){
  ag$optimize(ngen=as.integer(ngen),weights=as.integer(weights),objectives=objectives,verbose=verbose,
              seed=as.integer(seed),mode=mode,population_size=as.integer(population_size),
              offspring_size=as.integer(offspring_size),crossover_pb=crossover_pb,mutation_pb=mutation_pb,
              crossover_thres=as.integer(crossover_thres),
              ind_standard_pb=ind_standard_pb)
  } else if (mode=="fixed"){
    ag$optimize(ngen=as.integer(ngen),weights=as.integer(weights),objectives=objectives,verbose=verbose,
                nfeatures=as.integer(nfeatures),seed=as.integer(seed),mode=mode,population_size=as.integer(population_size),
                offspring_size=as.integer(offspring_size),crossover_pb=crossover_pb,mutation_pb=mutation_pb,
                mutate_flip_pb=mutate_flip_pb)
  } else {
    base::stop(paste0("Mode ",mode," not recognized. Please try 'standard' or 'fixed'"))
  }

  if (plot){
    print("plotting")
    if (sum(!sapply(list(weights,index,close_to),is.null))>1){
      base::stop('When selecting a solution in autogenes only one of the parameters "plot_weights", "index" and "close_to" can be used. It is also possible to use none of them')
    }
    plot_objectives <- as.integer(plot_objectives)
    if (!is.null(plot_weights)){
      ag$plot(objectives=plot_objectives,weights=as.integer(plot_weights))
    } else if (!is.null(index)){
      ag$plot(objectives=plot_objectives,index=as.integer(index))
    } else if (!is.null(close_to)){
      ag$plot(objectives=plot_objectives,close_to=as.integer(close_to))
    } else {
      print("befoire")
      ag$plot(objectives=plot_objectives)
      print("after")
    }
  }
  filename <- tempfile(fileext = ".pickle")
  ag$save(filename)
  if (verbose){
    base::message("Successfully saved")
  }
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

  if (sum(!sapply(list(weights,index,close_to),is.null))>1){
    base::stop('When selecting a solution in autogenes only one of the parameters "weights", "index" and "close_to" can be used')
  }
  ag <- reticulate::import("autogenes")
  ag$load(signature)

  if (!is.null(weights)){
    ag$select(weights=weights)
  } else if (!is.null(index)){
    ag$select(index=index)
  } else if (!is.null(close_to)){
    ag$select(close_to=close_to)
  } else {
    ag$select()
  }

  selection <- ag$selection()
  genes <- ag$adata()$var_names
  genes_used <- genes[selection]
  if (verbose){
    base::message(length(genes_used),"genes are used for the deconvolution")
  }

  bulk_data <- t(bulk_gene_expression)
  result <- ag$deconvolve(bulk_data, model=model, nu=nu, C=C,
                        kernel=kernel, degree=as.integer(degree),gamma=gamma, coef0=coef0,
                        shrinking=shrinking, tol=tol,cache_size=as.integer(cache_size),
                        verbose=verbose, max_iter=as.integer(max_iter))

  colnames(result) <- ag$adata()$obs_names
  rownames(result) <- rownames(bulk_data)
  return(list(proportions=result,genes_used=genes_used))
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

#' Checks if python and the autogenes module are available and installs them if they are not.
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
}
