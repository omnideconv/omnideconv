#' Install Scaden
#'
#' Install Scaden into a previously defined python environment.
#' A python environment can be defined by set_virtualenv() or set_python() command.
#' Alternatively a new environment can be created via create_virtualenv() method.
#'
#' @export
#'
#' @examples

install_scaden <- function() {
  reticulate::py_install("scaden",pip=TRUE)

}

#' Builds Scaden Model from scRNA data
#'
#' The model is saved in a defined directory.
#'
#' @param celltype_labels Vector of celltype labels. Order corresponds to rows in single_cell_object matrix.
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes, columns are samples.
#' @param bulk_data Matrix of bulk RNA expression. (rows :genes, columns: samples)
#' @param model_path Path where model directory should be created. (optional)
#' @param samples Bulk simulation: Number of samples to simulate (default: 1000)
#' @param cells Bulk simulation: Number of cells per sample (default: 100)
#' @param dataset_name Bulk simulation: Name of simulated dataset (default scaden)
#' @param var_cutoff Training data processing: Filter out genes with a variance less than the specified cutoff. A low cutoff is recommended,this should only remove genes that are obviously uninformative. (default NULL)
#' @param batch_size Training  of model: Batch size to use for training. (default: 128)
#' @param learning_rate Training of model: Learning rate used for training. (default: 0.0001)
#' @param steps Training of model: Number of training steps (default: 5000)
#' @param verbose If true: Command line output of Scaden is shown (default: false)
#'
#' @return Scaden model
#' @export
#'
#' @examples
scaden_build_model <- function(single_cell_object ,celltype_labels ,bulk_data=NULL ,model_path=NULL, batch_size=128, learning_rate= 0.0001, steps=5000, var_cutoff=NULL, cells=100, samples=1000, dataset_name="scaden", verbose=FALSE){

  if (!verbose){
    if (Sys.info()['sysname']=="Windows"){
      base::message("The windows implementation requires verbose mode. It is now switched on.")
      verbose <- TRUE
    }
  }

  # check if Scaden is installed and loaded.
  scaden_checkload()

  single_cell_object <- t(single_cell_object)
  if (nrow(single_cell_object) != length(celltype_labels)){
    base::stop("Celltype labels must be same length as samples (number of columns) in single_cell_object!")
  }
  gene_labels <- colnames(single_cell_object)


  if (!is.null(bulk_data)){

    training_h5ad <- scaden_simulate(celltype_labels = celltype_labels, gene_labels = gene_labels, single_cell_object = single_cell_object, cells = cells, samples = samples, dataset_name = dataset_name, verbose=verbose)
    processed <- scaden_process(training_h5ad,bulk_data, verbose=verbose, var_cutoff = var_cutoff)
    output_model_path <- scaden_train(processed,model_path=model_path, verbose = verbose, batch_size = batch_size , learning_rate = learning_rate, steps = steps)
    return(output_model_path)
  }
  else{
    base::stop("Scaden needs bulk RNA data for building the model. Please forward a bulk RNA file with build_model(..., bulk_data= your_bulk_data )")
  }
}


#' Performs deconvolution
#'
#' @param model Path to the model directory
#' @param bulk_data Matrix of bulk RNA expression. (genes x samples)
#' @param verbose Defines verbosity of function call (default: false)
#'
#' @return Cell type fractions per sample.
#' @export
#'
#' @examples
deconvolute_scaden <- function(model,bulk_data, verbose=FALSE){
  if (!verbose){
    if (Sys.info()['sysname']=="Windows"){
      base::message("The windows implementation requires verbose mode. It is now switched on.")
      verbose <- TRUE
    }
  }

  scaden_checkload()
  prediction <- scaden_predict(model,bulk_data, verbose = verbose)
  return(prediction)
}


#' Trains a Scaden model
#'
#' Trains a Scaden model on pre-processed training data.
#'
#' @param h5ad_processed pre-processed training data
#' @param batch_size Batch size to use for training. Default: 128
#' @param learning_rate Learning rate used for training. Default: 0.0001
#' @param model_path Path to store the model in
#' @param steps Number of training steps
#' @param verbose Defines verbosity of function call (default: false)
#'
#' @return Scaden model
#' @export
#'
#' @examples
scaden_train <- function(h5ad_processed, batch_size=128, learning_rate= 0.0001, model_path=NULL, steps=5000, verbose=FALSE){

  base::message("Training model")

  # create temporary directory where Scaden input files should be saved at.
  tmp_dir <- tempdir()
  dir.create(tmp_dir,showWarnings = FALSE)

  # the file is only created temporarily, later deleted at unlink(tmp)
  h5ad_processed_tmp <- tempfile(tmpdir = tmp_dir)
  write_anndata(h5ad_processed,h5ad_processed_tmp)

  if (is.null(model_path)){
    currentwd <- base::getwd()
    model_path <- tryCatch(
      {
        base::setwd(tmp_dir)
        model_path <- paste0(tmp_dir,"/model")
        if (dir.exists(model_path)){
          unlink(model_path, recursive = TRUE)
        }
        dir.create(model_path, showWarnings = FALSE)
        model_path

      },
      error=function(cond) {
        base::setwd(currentwd)
        base::stop(cond)
      },
      warning=function(cond) {
        base::warning(cond)
      },
      finally={
        base::setwd(currentwd)
      }
    )
  }

  # Calling Scaden command
  system(paste("scaden train",h5ad_processed_tmp,"--batch_size",batch_size,"--learning_rate",learning_rate,"--steps",steps,"--model_dir",model_path), ignore.stdout = !verbose, ignore.stderr = !verbose)

  model_dirs <- list.dirs(path = model_path, full.names = FALSE, recursive = FALSE)
  model_names <- c('m1024','m256','m512')
  for (model in model_names){
    if (!(model %in% model_dirs)){
      base::stop("Model generation failed!")
    }
  }

  return(model_path)
}



#' Pre-processes training data
#'
#' Training data needs to be processed before a model can be trained.
#'
#' @param h5ad File that should be processed. Must be in AnnData format (.h5ad)
#' @param bulk_data Bulk RNA-seq data. (genes x individuals)
#' @param var_cutoff Filter out genes with a variance less than the specified cutoff. A low cutoff is recommended,this should only remove genes that are obviously uninformative.
#' @param verbose Defines verbosity of function call (default: false)
#'
#' @return processed training file. (.h5ad format)
#' @export
#'
#' @examples
scaden_process <- function(h5ad,bulk_data,var_cutoff=NULL, verbose=FALSE){

  base::message("Processing training data for model creation ...")



  out <- tryCatch(
    {
      tmp_dir <- tempdir()
      dir.create(tmp_dir,showWarnings = FALSE)

      h5ad_tmp <- tempfile(tmpdir = tmp_dir,fileext =".h5ad")
      write_anndata(h5ad,h5ad_tmp)

      bulk_data_tmp <- tempfile(tmpdir = tmp_dir)
      utils::write.table(bulk_data,file = bulk_data_tmp, sep = "\t",row.names = TRUE,col.names = NA,quote = FALSE)

      processed_h5ad <- tempfile(fileext =".h5ad",tmpdir = tmp_dir)


      if(is.null(var_cutoff)){
        system(paste("scaden process",h5ad_tmp,bulk_data_tmp,"--processed_path",processed_h5ad), ignore.stdout = !verbose, ignore.stderr = !verbose)
      }
      else{
        system(paste("scaden process",h5ad_tmp,bulk_data_tmp,"--processed_path",processed_h5ad,"--var_cutoff",var_cutoff), ignore.stdout = !verbose, ignore.stderr = !verbose)
      }
      read_anndata(processed_h5ad)
    },
    error=function(cond){
      base::message("Error preprocessing training data! Make sure training data is not in logarithmic space!")
    },
    warning = function(cond){
      if (verbose){
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
#' @param bulk_data Bulk RNA-seq data. (genes x bulk_RNA samples)
#' @param verbose Defines verbosity of function call (default: false)
#'
#' @return Cell type fractions per sample
#' @export
#'
#' @examples
scaden_predict <- function(model_dir, bulk_data, verbose=FALSE){

  base::message("Predicting cell type proportions")


  current_wd <- base::getwd()

  predictions <- tryCatch(
    {
      tmp_dir <- tempdir()
      dir.create(tmp_dir,showWarnings = FALSE)
      base::setwd(tmp_dir)

      bulk_data_tmp <- tempfile(tmpdir = tmp_dir)
      utils::write.table(bulk_data,file = bulk_data_tmp, sep = "\t",row.names = TRUE,col.names = NA,quote = FALSE)

      system(paste("scaden predict --model_dir",model_dir,bulk_data_tmp), ignore.stdout = !verbose, ignore.stderr = !verbose)

      t(utils::read.table(paste0(tmp_dir,"/scaden_predictions.txt"),sep = "\t",header = TRUE, row.names = 1))

    },
    error=function(cond) {
      base::stop(cond)
    },
    warning=function(cond) {
      base::warning(cond)
    },
    finally={
      base::setwd(current_wd)
    }
  )

  colnames(predictions)[1]<-" "

  return(predictions)

}

#' Simulates example Data provided by Scaden
#'
#' Used for testing.
#'
#' @param example_data_path Path to where example data should be saved. (directory)
#' @param verbose Defines verbosity of function call (default: false)
#'
#' @return List with list$simulated_h5ad =  example training data
#' and list$bulk = example bulk data.
#' @export
#'
#' @examples
scaden_simulate_example <- function(example_data_path=NULL, verbose=FALSE){


  current_wd <- base::getwd()

  if (is.null(example_data_path)){
    tmp_dir <- tempdir()
    dir.create(tmp_dir,showWarnings = FALSE)
    base::setwd(tmp_dir)
  }
  else{
    base::setwd(example_data_path)
  }


  if (!verbose){
    logfile <- tempfile(tmpdir = tmp_dir,fileext = "train.log")
    sink(file = logfile)
  }


  system("mkdir example_data")
  system("scaden example --out example_data/")
  system(paste0("scaden simulate --data ",tmp_dir,"/example_data/ -n 100 --pattern *_counts.txt"))

  simulated_h5ad <- read_anndata(paste0(tmp_dir,"/data.h5ad"))
  bulk <- utils::read.table(paste0(tmp_dir,"/example_data/example_bulk_data.txt"))

  base::setwd(current_wd)
  unlink(tmp_dir)

  output <- list("simulated_h5ad"=simulated_h5ad,"bulk"=bulk)
  return(output)

}

#' Simulates training data from scRNA data
#'
#' @param celltype_labels Vector of celltype labels. Order corresponds to rows in single_cell_object matrix.
#' @param gene_labels Vector of gene labels. Order corresponds to columns in single_cell_object matrix.
#' @param single_cell_object Matrix or dataframe of scRNA data. Rows=cells and columns=genes
#' @param samples Bulk simulation: Number of samples to simulate (default: 1000)
#' @param cells Bulk simulation: Number of cells per sample (default: 100)
#' @param dataset_name Name of dataset
#' @param verbose Defines verbosity of function call (default: false)
#'
#' @return Simulated bulk data of known cell type fractions
#' @export
#'
#' @examples
scaden_simulate <- function(celltype_labels ,gene_labels , single_cell_object, cells=100, samples=1000, dataset_name="scaden" , verbose = FALSE){



    base::message("Simulating training data from single cell experiment: ",samples, " samples of ",cells, " cells")

    current_wd <- base::getwd()


    output <- tryCatch(
      {
        tmp_dir <- tempdir()
        dir.create(tmp_dir,showWarnings = FALSE)
        base::setwd(tmp_dir)
        if (dir.exists(dataset_name)){
          unlink(dataset_name, recursive = TRUE)
        }
        dir.create(dataset_name,showWarnings = FALSE)
        base::setwd(paste0(tmp_dir,"/",dataset_name))

        colnames(single_cell_object)<-gene_labels
        cell_types <- data.frame("Celltype"=celltype_labels)

        utils::write.table(single_cell_object,paste0(tmp_dir,"/",dataset_name,"/",dataset_name,"_counts.txt") ,sep = "\t",row.names = TRUE,col.names = NA,quote = FALSE)
        utils::write.table(cell_types,paste0(tmp_dir,"/",dataset_name,"/",dataset_name,"_celltypes.txt") ,quote = FALSE,row.names = FALSE,col.names = TRUE)
        base::setwd(tmp_dir)

        system(paste("scaden simulate --data",paste0(tmp_dir,"/",dataset_name),"-n",samples,"-c",cells,"--pattern *_counts.txt"), ignore.stdout = !verbose, ignore.stderr = !verbose)


        read_anndata(paste0(tmp_dir,"/data.h5ad"))
      },
      error=function(cond) {
        base::stop(cond)
      },
      warning=function(cond) {
        base::warning(cond)
      },
      finally={
        base::setwd(current_wd)
      }
    )
    return(output)
}

#' Checks if scaden is installed.
#'
#' If it is available, the python module is imported.
#'
#' @export
#'
#' @examples
scaden_checkload <- function(){
  if (python_available()){
    if (reticulate::py_module_available("scaden")){
      reticulate::import("scaden")
    }
    else{
      install_scaden()
      reticulate::import("scaden")
    }
  }
  else{
    base::message("Setting up python environment..")
    init_python()
    if (python_available()){
      install_scaden()
      reticulate::import("scaden")
    }
    else{
      base::stop("Could not initiate miniconda python environment. Please set up manually with init_python(python=your/python/version)")
    }
  }
}
