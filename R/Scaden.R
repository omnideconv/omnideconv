#' Install Scaden
#'
#' Install Scaden into a previously defined python environment.
#' A python environment can be defined by set_virtualenv() or set_python() command.
#' Alternatively a new environment can be created via create_virtualenv() method.
#'
#' @param method
#' @param conda
#'
#' @return
#' @export
#'
#' @examples
install_scaden <- function(method = "auto", conda = "auto") {
  reticulate::py_install("scaden", method = method, conda = conda)
}

#' Builds Scaden Model from scRNA data
#'
#' The model is saved in a defined directory.
#'
#' @param celltype_labels Vector of celltype labels. Order corresponds to rows in count_data matrix.
#' @param gene_labels Vector of gene labels. Order corresponds to columns in count_data matrix.
#' @param count_data Matrix or dataframe of scRNA data. Rows=cells and columns=genes
#' @param bulk_data Matrix of bulk RNA expression. (genes x samples)
#' @param model_path Path where model directory should be created.
#'
#' @return
#' @export
#'
#' @examples
scaden_build_model <- function(celltype_labels ,gene_labels , count_data,bulk_data,model_path){
  # Untraceable error arises if a model is already saved at indicated model path
  if (dir.exists(model_path)){
    base::stop("Either choose another model path or delete existing directory!")
  }
  training_h5ad <- scaden_simulate(celltype_labels,gene_labels,count_data)
  processed <- scaden_process(training_h5ad,bulk_data)
  scaden_train(processed,model_path=model_path)
  return(model_path)
}

#' Builds Scaden model from training data in .h5ad format.
#'
#' @param training_data_h5ad Training data in .h5ad format.
#' @param bulk_data Matrix of bulk RNA expression. (genes x samples)
#' @param model_path Path where model directory should be created.
#'
#' @return
#' @export
#'
#' @examples
scaden_build_model_from_h5ad <- function(training_data_h5ad,bulk_data,model_path){
  # Untraceable error arises if a model is already saved at indicated model path
  if (dir.exists(model_path)){
    base::stop("Either choose another model path or delete existing directory!")
  }
  processed  <- scaden_process(training_data_h5ad,bulk_data)
  scaden_train(processed,model_path = model_path)
  return(model_path)
}

#' Performs deconvolution
#'
#' @param model Path to the model directory
#' @param bulk_data Matrix of bulk RNA expression. (genes x samples)
#'
#' @return
#' @export
#'
#' @examples
scaden_deconvolute <- function(model,bulk_data){
  prediction <- scaden_predict(model,bulk_data)
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
#'
#' @return
#' @export
#'
#' @examples
scaden_train <- function(h5ad_processed, batch_size=128, learning_rate= 0.0001, model_path="model", steps=5000){
  # check if Scaden is installed and loaded.
  scaden_checkload()

  # create temporary directory where Scaden input files should be saved at.
  tmp_dir <- tempdir()
  dir.create(tmp_dir,showWarnings = F)

  # the file is only created temporarily, later deleted at unlink(tmp)
  h5ad_processed_tmp <- tempfile(tmpdir = tmp_dir)
  write_h5ad(h5ad_processed,h5ad_processed_tmp)

  # Calling Scaden command
  system(paste("scaden train",h5ad_processed_tmp,"--batch_size",batch_size,"--learning_rate",learning_rate,"--steps",steps,"--model_dir",model_path))

  # removal and deletion of temporary folder and it's content.
  unlink(tmp_dir)
  return(model_path)
}



#' Pre-processes training data
#'
#' Training data needs to be processed before a model can be trained.
#'
#' @param h5ad File that should be processed. Must be in AnnData format (.h5ad)
#' @param bulk_data Bulk RNA-seq data. (genes x individuals)
#' @param var_cutoff Filter out genes with a variance less than the specified cutoff. A low cutoff is recommended,this should only remove genes that are obviously uninformative.
#'
#' @return
#' @export
#'
#' @examples
scaden_process <- function(h5ad,bulk_data,var_cutoff=NULL){
  # see scaden_train comments for workflow explanation
  scaden_checkload()

  tmp_dir <- tempdir()
  dir.create(tmp_dir,showWarnings = F)

  h5ad_tmp <- tempfile(tmpdir = tmp_dir,fileext =".h5ad")
  h5ad$write(h5ad_tmp)

  bulk_data_tmp <- tempfile(tmpdir = tmp_dir)
  write.table(bulk_data,file = bulk_data_tmp, sep = "\t",row.names = T,col.names = NA,quote = F)

  processed_h5ad <- tempfile(fileext =".h5ad",tmpdir = tmp_dir)

  if(is.null(var_cutoff)){
    system(paste("scaden process",h5ad_tmp,bulk_data_tmp,"--processed_path",processed_h5ad))
  }
  else{
    system(paste("scaden process",h5ad_tmp,bulk_data_tmp,"--processed_path",processed_h5ad,"--var_cutoff",var_cutoff))
  }


  output_h5ad <- read_h5ad(processed_h5ad)
  unlink(tmp_dir)

  return(output_h5ad)
}

#' Predict cell proportions
#'
#' Predicts cell proportions in bulk RNA sample.
#'
#' @param model_dir Directory where model is saved
#' @param bulk_data Bulk RNA-seq data. (genes x bulk_RNA samples)
#'
#' @return
#' @export
#'
#' @examples
scaden_predict <- function(model_dir, bulk_data){
  # see scaden_train comments for workflow explanation
  scaden_checkload()


  current_wd <- getwd()

  tmp_dir <- tempdir()
  dir.create(tmp_dir,showWarnings = F)
  setwd(tmp_dir)

  bulk_data_tmp <- tempfile(tmpdir = tmp_dir)
  write.table(bulk_data,file = bulk_data_tmp, sep = "\t",row.names = T,col.names = NA,quote = F)

  system(paste("scaden predict --model_dir",model_dir,bulk_data_tmp))

  predictions <- read.table(paste0(tmp_dir,"/scaden_predictions.txt"))

  unlink(tmp_dir)
  setwd(current_wd)
  return(predictions)

}

#' Simulates example Data provided by Scaden
#'
#' Used for testing.
#'
#' @param example_data_path Path to where example data should be saved. (directory)
#'
#' @return List with list$simulated_h5ad =  example training data
#' and list$bulk = example bulk data.
#' @export
#'
#' @examples
scaden_simulate_example <- function(example_data_path=NULL){
  # see scaden_train comments for workflow explanation
  scaden_checkload()

  current_wd <- getwd()

  if (is.null(example_data_path)){
    tmp_dir <- tempdir()
    dir.create(tmp_dir,showWarnings = F)
    setwd(tmp_dir)
  }
  else{
    setwd(example_data_path)
  }



  system("mkdir example_data")
  system("scaden example --out example_data/")
  system(paste0("scaden simulate --data ",tmp_dir,"/example_data/ -n 100 --pattern *_counts.txt"))

  simulated_h5ad <- read_h5ad(paste0(tmp_dir,"/data.h5ad"))
  bulk <- read.table(paste0(tmp_dir,"/example_data/example_bulk_data.txt"))

  setwd(current_wd)
  unlink(tmp_dir)

  output <- list("simulated_h5ad"=simulated_h5ad,"bulk"=bulk)
  return(output)

}

#' Simulates training data from scRNA data
#'
#' @param celltype_labels Vector of celltype labels. Order corresponds to rows in count_data matrix.
#' @param gene_labels Vector of gene labels. Order corresponds to columns in count_data matrix.
#' @param count_data Matrix or dataframe of scRNA data. Rows=cells and columns=genes
#'
#' @return
#' @export
#'
#' @examples
scaden_simulate <- function(celltype_labels ,gene_labels , count_data, cells=100, samples=1000){

    scaden_checkload()

    current_wd <- getwd()

    tmp_dir <- tempdir()
    dir.create(tmp_dir,showWarnings = F)
    setwd(tmp_dir)

    colnames(count_data)<-gene_labels
    cell_types <- data.frame("Celltype"=celltype_labels)

    write.table(count_data,paste0(tmp_dir,"/_counts.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
    write.table(cell_types,paste0(tmp_dir,"/_celltypes.txt"),quote = F,row.names = F,col.names = T)

    ftmpdir <- paste0(tmp_dir,"/")

    system(paste("scaden simulate --data",ftmpdir,"-n",samples,"-c",cells,"--pattern *_counts.txt"))

    output <- read_h5ad(paste0(tmp_dir,"/data.h5ad"))

    setwd(current_wd)
    unlink(tmp_dir)
    return(output)


}

#' Checks if Scaden is installed.
#'
#' If it is available, the python module is imported.
#'
#' @return
#' @export
#'
#' @examples
scaden_checkload <- function(){
  if (reticulate::py_module_available("scaden")){
    reticulate::import("scaden")
  }
  else{
    base::stop("python module scaden not available in environment! Run install_scaden() to install it.")
  }
}
