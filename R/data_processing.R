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
  sc_pheno <- data.frame(check.names=F, check.rows=F,
                         stringsAsFactors=F,
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


AnnData_to_SingleCellExperiment <- function(ad){
  ad <- ad$transpose()
  X_mat <- ad$X
  rownames(X_mat) <- ad$obs_names
  colnames(X_mat) <- ad$var_names


  assays_list <- list(X = X_mat)
  layer_names <- names(ad$layers)
  if (length(ad$layers)>0){
    for (i in 1:length(ad$layers)) {
      layer_name <- layer_names[i]
      assays_list[[layer_name]] = ad$layers[[i]]
    }
  }

  meta_list <- list()
  uns_keys <- ad$uns_keys()
  if (length(uns_keys)>0){
    for (i in 1:length(uns_keys)) {
      key <- uns_keys[i]
      meta_list[[key]] = ad$uns[[i]]
    }
  }

  rowdata <- ad$obs
  rowdata$obsm <- ad$obsm

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



SingleCellExperiment_to_Anndata <- function(sce, X_name = NULL){
  if (is.null(X_name)) {
    if (length(assays(sce)) == 0) {
      stop("'sce' does not contain any assays")
    }
    X_name <- assayNames(sce)[1]
    message("Note: using the '", X_name, "' assay as the X matrix")
  }
  X <- assay(sce,X_name)


  col_data <- colData(sce)
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

  row_data <- rowData(sce)
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

  assay_names <- assayNames(sce)
  assay_names <- assay_names[!assay_names == X_name]
  assays_list <- list()
  for (name in assay_names) {
    assays_list[[name]] = assay(sce,name)
  }


  meta_list <- metadata(sce)
  uns_list <- list()
  for (item_name in names(meta_list)) {
    item <- meta_list[[item_name]]
    uns_list[[item_name]] <- item
  }


  #ad$obsp <- as.list(rowPairs(sce))
  #ad$varp <- as.list(colPairs(sce))

  ad <- anndata::AnnData(X=X, obs=obs, var=var, uns=uns_list, varm = as.list(reducedDims(sce)), obsm = as.list(rowData(sce)$obsm), layers = assays_list)

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

Matrix_to_SingleCellExperiment <- function(matrix, cell_labels, named_metadata_list=list()){
  gene_labels <- rownames(matrix)
  sce <- SingleCellExperiment(list(X=matrix),
                              colData=DataFrame(label=cell_labels),
                              rowData=DataFrame(length=gene_labels),
                              metadata=named_metadata_list
  )
  return(sce)
}

SingleCellExperiment_to_Matrix <- function(sce, assay_name=NULL){
  if (is.null(assay_name)) {
    if (length(assays(sce)) == 0) {
      stop("'sce' does not contain any assays")
    }
    assay_name <- assayNames(sce)[1]
    message("Note: using the '", assay_name, "' assay as the X matrix")
  }
  X <- assay(sce,assay_name,withDimnames = T)
  cell_labels <- colData(sce)$label

  return(list("matrix"=X, "annotation_vector"=cell_labels))
}



