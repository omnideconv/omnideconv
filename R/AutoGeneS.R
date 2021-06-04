#' Calculates the signature model with AutoGeneS
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes, columns are samples. Row and column names need to be set.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples. Necessary for MOMF, defaults to NULL.Row and column names need to be set.
#'
#' @return The path to the pickle file needed for the deconvolution with autogenes.
#' @export
#'
#'
build_model_autogenes <- function(single_cell_object, cell_type_annotations, bulk_gene_expression = NULL,
                                  ngen=2,mode="standard",nfeatures=NULL,weights=c(-1,1),objectives=c("correlation","distance"),seed=0,verbose=FALSE,population_size=100,
                                  offspring_size=50,crossover_pb=0.7,mutation_pb=0.3,mutate_flip_pb=1E-3,crossover_thres=1000,ind_standard_pb=0.1,plot_weights=NULL,
                                  plot_objectives=c(0,1),index=NULL,close_to=NULL,plot=FALSE){
  if (!is.null(bulk_gene_expression)){
    single_cell_object <- single_cell_object[intersect(rownames(single_cell_object),rownames(bulk_gene_expression)),]
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
#' @param signature Signature Matrix (genes x individuals from scRNA-seq)
#' @param single_cell_object scRNA-seq Object (genes x cells)
#' @param method Determines which divergence to use. Options: Kullback-Leibler "KL", Itakura-Saito "IS". Defaults to "KL"
#' @param verbose Whether the algorithm should print out what it is doing.
#' @param ... additional parameters
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

get_buildmodel_parameters_autogenes <- function(ngen=2,mode="standard",nfeatures=NULL,weights=c(-1,1),objectives=c("correlation","distance"),seed=0,verbose=FALSE,population_size=100,
                                      offspring_size=50,crossover_pb=0.7,mutation_pb=0.3,mutate_flip_pb=1E-3,crossover_thres=1000,ind_standard_pb=0.1,plot_weights=NULL,
                                      plot_objectives=c(0,1),index=NULL,close_to=NULL,plot=FALSE){
  string <- paste("-ngen",ngen,"-mode",mode,"-nfeatures",nfeatures,"-weights",paste(weights,collapse = " "),"-objectives",paste(objectives,collapse = " "),"-seed",seed,
                  "-verbose",verbose,"-population_size",population_size,"-offspring_size",offspring_size,"-crossover_pb",crossover_pb,"-mutation_pb",mutation_pb,
                  "-mutate_flip_pb",mutate_flip_pb,"-crossover_thres",crossover_thres,"-ind_standard_pb",ind_standard_pb,"-plot",plot,"-plot_objectives",
                  paste(plot_objectives,collapse = " "))
  if (plot_weights){
    string <- paste(string,"-plot_weights",paste(plot_weights,collapse = " "))
  }
  if (index){
    string <- paste(string,"-index",paste(index,collapse = " "))
  }
  if (close_to){
    string <- paste(string,"-close_to",paste(close_to,collapse = " "))
  }
  return(string)
}

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
transform_boolean <- function(bool){
  if(bool){
    return("True")
  }
  return("False")
}
quote_string <- function(string){
  return(paste0('"',string,'"'))
}
null_transformer <- function(str = "NULL") {
  function(text, envir) {
    out <- glue::identity_transformer(text, envir)
    if (is.null(out)) {
      return(str)
    }
    out
  }
}
