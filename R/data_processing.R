#' Generating a single cell expression set
#'
#' @param single_cell_matrix The single-cell matrix. Rows are genes, columns are samples.
#' @param sample_names A vector of the names of the samples, basically colnames(single_cell_matrix).
#' @param genes A vector of the names of the genes, basically rownames(single_cell_matrix).
#' @param cell_types A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object.
#'
#' @return A Biobase::ExpressionSet of the input data.
#' @export
#' @importClassesFrom Biobase AnnotatedDataFrame
#'
#' @examples
get_single_cell_expression_set <- function(single_cell_matrix, sample_names, genes, cell_types){


  # individual.ids and cell.types should be in the same order as in sampleNames
  sc_pheno <- data.frame(check.names=FALSE, check.rows=FALSE,
                         stringsAsFactors=FALSE,
                         row.names=sample_names,
                         SubjectName=sample_names,
                         cellType=cell_types)
  sc_meta <- data.frame(labelDescription=c("SubjectName",
                                           "cellType"),
                        row.names=c("SubjectName",
                                    "cellType"))
  sc_pdata <- new("AnnotatedDataFrame",
                  data=sc_pheno,
                  varMetadata=sc_meta)
  colnames(single_cell_matrix) <- row.names(sc_pdata)
  rownames(single_cell_matrix) <- genes
  return(Biobase::ExpressionSet(assayData=single_cell_matrix,
                                phenoData=sc_pdata))
}


#' Convert AnnData to SingleCellExperiment
#'
#' @param ad AnnData object
#'
#' @return SingleCellObject
#' @export
#'
#' @examples
anndata_to_singlecellexperiment <- function(ad){
  ad <- ad$transpose()
  X_mat <- ad$X
  rownames(X_mat) <- ad$obs_names
  colnames(X_mat) <- ad$var_names


  assays_list <- list()
  assays_list[['X']]<- X_mat
  assays_list <- c(assays_list,ad$layers)

  meta_list <- list()
  meta_list[ad$uns_keys()]<-ad$uns

  rowdata <- ad$obs
  if (length(ad$obsm) == nrow(rowdata)){
    rowdata$obsm <- ad$obsm
  }


  sce <- SingleCellExperiment::SingleCellExperiment(
    assays      = assays_list,
    rowData     = rowdata,
    colData     = ad$var,
    reducedDims = ad$varm,
    metadata    = meta_list,
    rowPairs    = ad$obsp,
    colPairs    = ad$varp
  )

  return(sce)
}



#' Convert SingleCellExperiment to AnnData object
#'
#'
#' @param sce object of type SingleCellExperiment
#' @param X_name Name of the assay in the SingleCellExperiment to interpret as the X matrix of the AnnData object (default: first in assay list)
#'
#' @return AnnData object
#' @export
#'
#' @examples
singlecellexperiment_to_anndata <- function(sce, X_name = NULL){
  if (is.null(X_name)) {
    if (length(SummarizedExperiment::assays(sce)) == 0) {
      stop("'sce' does not contain any assays")
    }
    X_name <- SummarizedExperiment::assayNames(sce)[1]
    message("Note: using the '", X_name, "' assay as the X matrix")
  }
  X <- SummarizedExperiment::assay(sce,X_name)


  col_data <- SingleCellExperiment::colData(sce)
  var <- NULL
  if (ncol(col_data) > 0) {
    # Manually construct the data.frame to avoid mangling column names
    var <- do.call(
      data.frame,
      c(
        as.list(col_data),
        check.names      = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }

  row_data <- SingleCellExperiment::rowData(sce)
  row_data <- row_data[1]
  obs <- NULL
  if (ncol(row_data) > 0) {
    # Manually construct the data.frame to avoid mangling column names
    obs <- do.call(
      data.frame,
      c(
        as.list(row_data),
        check.names      = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }


  assay_names <- SummarizedExperiment::assayNames(sce)
  assay_names <- assay_names[!assay_names == X_name]
  assays_list <- list()
  for (name in assay_names) {
    assays_list[[name]] = SummarizedExperiment::assay(sce,name)
  }


  meta_list <- S4Vectors::metadata(sce)
  uns_list <- list()
  for (item_name in names(meta_list)) {
    item <- meta_list[[item_name]]
    uns_list[[item_name]] <- item
  }


  #ad$obsp <- as.list(rowPairs(sce))
  #ad$varp <- as.list(colPairs(sce))

  ad <- anndata::AnnData(X=X, obs=obs, var=var, uns=uns_list, varm = as.list(SingleCellExperiment::reducedDims(sce)), obsm = as.list(SingleCellExperiment::rowData(sce)$obsm), layers = assays_list)

  if (!is.null(colnames(sce))) {
    ad$var_names <- colnames(sce)
  }

  if (!is.null(rownames(sce))) {
    ad$obs_names <- rownames(sce)
  }

  ad <- ad$transpose()
  return(ad)
}


# SingleCellExperiment_to_Seurat <- function(sce){
#   seu <- Seurat::as.Seurat(sce)
#   return(seu)
# }
#
# Seurat_to_SingleCellExperiment <- function(seu){
#   sce <- Seurat::as.SingleCellExperiment(seu)
#   return(sce)
# }

#' Convert count matrix to SingleCellExperiment
#'
#' @param matrix Matrix of single cell data (rows: genes, columns: cells), the row names of the matrix should be labeled with gene labels
#' @param cell_labels Annotation vector for celltypes (e.g columns of matrix)
#' @param named_metadata_list A named list containing unstructured metadata information (optional)
#'
#' @return SingleCellExperiment
#' @export
#'
#' @examples
matrix_to_singlecellexperiment <- function(matrix, cell_labels, named_metadata_list=list()){
  gene_labels <- rownames(matrix)
  sce <- SingleCellExperiment::SingleCellExperiment(list(X=matrix),
                              colData= S4Vectors::DataFrame(label=cell_labels),
                              rowData=S4Vectors::DataFrame(length=gene_labels),
                              metadata=named_metadata_list
  )
  return(sce)
}

#' Convert SingleCellExperiment to Matrix + annotation vector for celltypes
#'
#' @param sce SingleCellExperiment
#' @param assay_name name of the assay that should be returned as matrix. (default: first in assay list)
#'
#' @return named list: ..$matrix: matrix object, ..$annotation_vector
#' @export
#'
#' @examples
singlecellexperiment_to_matrix <- function(sce, assay_name=NULL){
  if (is.null(assay_name)) {
    if (length(SummarizedExperiment::assays(sce)) == 0) {
      stop("'sce' does not contain any assays")
    }
    assay_name <- SummarizedExperiment::assayNames(sce)[1]
    message("Note: using the '", assay_name, "' assay as the X matrix")
  }
  X <- SummarizedExperiment::assay(sce,assay_name,withDimnames = T)
  cell_labels <- SingleCellExperiment::colData(sce)$label

  return(list("matrix"=X, "annotation_vector"=cell_labels))
}

#' Check if two SingleCellExperiments are identical
#'
#' @param a SingleCellExperiments
#' @param b SingleCellExperiments
#'
#' @return boolean
#' @export
#'
#' @examples
sces_are_identical <- function(a,b){
  same <- TRUE
  assays_a <- as.list(SummarizedExperiment::assays(a))
  assays_b <- as.list(SummarizedExperiment::assays(b))
  if (length(assays_a) != length(assays_b)){
    same <- FALSE
    message("number of assays not identical")
  }
  else{


    if (!identical(rownames(a),rownames(b))){
      same<-FALSE
      message("rownames not identical")
    }

    if (!identical(colnames(a),colnames(b))){
      same<-FALSE
      message("rownames not identical")
    }

    a_names <- as.list(SummarizedExperiment::assayNames(a))
    b_names <- as.list(SummarizedExperiment::assayNames(b))
    for (i in 1:length(assays_a)) {
      m <- as.numeric(as.matrix(assays_a[[i]]))
      m2 <- as.numeric(as.matrix(assays_b[[i]]))

      if(!identical(m,m2)){
        same<-FALSE

        message(paste("Assay",a_names[[i]],"and assay",b_names[[i]]," not identical"))
      }
    }

    a_col <- base::as.data.frame(SummarizedExperiment::colData(a))
    b_col <- base::as.data.frame(SummarizedExperiment::colData(b))

    if (!identical(a_col,b_col)){
      same <- FALSE
      base::message("Coldata not identical")
    }

    a_row <- base::as.data.frame(SummarizedExperiment::rowData(a))
    b_row <- base::as.data.frame(SummarizedExperiment::rowData(b))

    if (!identical(a_row,b_row)){
      same <- FALSE
      base::message("Rowdata not identical")
    }

  }
  return(same)

}

