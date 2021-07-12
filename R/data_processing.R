#' Generating a single cell expression set
#'
#' @param single_cell_matrix The single-cell matrix. Rows are genes, columns are samples.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param genes A vector of the names of the genes, basically rownames(single_cell_matrix).
#' @param cell_types A Vector of the cell type annotations. Has to be in the same order as the
#'   samples in single_cell_object.
#'
#' @return A Biobase::ExpressionSet of the input data.
#' @importClassesFrom Biobase AnnotatedDataFrame
#'
get_single_cell_expression_set <- function(single_cell_matrix, batch_ids, genes, cell_types) {


  # individual.ids and cell.types should be in the same order as in sampleNames
  sc_pheno <- data.frame(
    check.names = FALSE, check.rows = FALSE,
    stringsAsFactors = FALSE,
    row.names = colnames(single_cell_matrix),
    batchId = batch_ids,
    cellType = cell_types
  )
  sc_meta <- data.frame(
    labelDescription = c(
      "batchId",
      "cellType"
    ),
    row.names = c(
      "batchId",
      "cellType"
    )
  )
  sc_pdata <- methods::new("AnnotatedDataFrame",
    data = sc_pheno,
    varMetadata = sc_meta
  )
  colnames(single_cell_matrix) <- row.names(sc_pdata)
  rownames(single_cell_matrix) <- genes
  return(Biobase::ExpressionSet(
    assayData = single_cell_matrix,
    phenoData = sc_pdata
  ))
}

#' Save as .h5ad
#'
#' Transform a single_cell_matrix into an anndata object and saves it
#'
#' @param single_cell_object A matrix or dataframe with the single-cell data. Rows are genes,
#'   columns are samples. Row and column names need to be set. Alternatively a SingleCellExperiment
#'   or an AnnData object can be provided. In that case, note that cell-type labels need to be
#'   indicated either directly providing a vector (cell_type_annotations) or by indicating the
#'   column name that indicates the cell-type labels (cell_type_column_name). (Anndata: obs object,
#'   SingleCellExperiment: colData object)
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order
#' as the samples in single_cell_object
#'
#' @return The path to the saved .h5ad file
#'
save_as_h5ad <- function(single_cell_object, cell_type_annotations) {
  sce <- matrix_to_singlecellexperiment(single_cell_object, cell_type_annotations)
  ad <- singlecellexperiment_to_anndata(sce)
  path <- tempfile(fileext = ".h5ad")
  write_anndata(ad, path)
  return(path)
}


#' Convert AnnData to SingleCellExperiment
#'
#' @param ad AnnData object
#'
#' @return SingleCellObject
#'
anndata_to_singlecellexperiment <- function(ad) {
  ad <- ad$transpose()
  X_mat <- ad$X
  rownames(X_mat) <- ad$obs_names
  colnames(X_mat) <- ad$var_names


  assays_list <- list()
  assays_list[["X"]] <- X_mat
  assays_list <- c(assays_list, lapply(ad$layers$keys(), function(x) ad$layers[x]))
  names(assays_list) <- c("X", ad$layers$keys())

  meta_list <- list()
  meta_list[ad$uns_keys()] <- ad$uns

  rowdata <- ad$obs
  if (length(ad$obsm) == nrow(rowdata)) {
    rowdata$obsm <- ad$obsm
  }


  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    rowData = rowdata,
    colData = ad$var,
    reducedDims = ad$varm,
    metadata = meta_list,
    rowPairs = ad$obsp,
    colPairs = ad$varp
  )

  return(sce)
}



#' Convert SingleCellExperiment to AnnData object
#'
#'
#' @param sce object of type SingleCellExperiment
#' @param X_name Name of the assay in the SingleCellExperiment to interpret as the X matrix of the
#' AnnData object (default: first in assay list)
#'
#' @return AnnData object
#'
singlecellexperiment_to_anndata <- function(sce, X_name = NULL) {
  if (is.null(X_name)) {
    if (length(SummarizedExperiment::assays(sce)) == 0) {
      stop("'sce' does not contain any assays")
    }
    X_name <- SummarizedExperiment::assayNames(sce)[1]
    message("Note: using the '", X_name, "' assay as the X matrix")
  }
  X <- SummarizedExperiment::assay(sce, X_name)


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
    assays_list[[name]] <- SummarizedExperiment::assay(sce, name)
  }


  meta_list <- S4Vectors::metadata(sce)
  uns_list <- list()
  for (item_name in names(meta_list)) {
    item <- meta_list[[item_name]]
    uns_list[[item_name]] <- item
  }


  # ad$obsp <- as.list(rowPairs(sce))
  # ad$varp <- as.list(colPairs(sce))

  ad <- anndata::AnnData(
    X = X, obs = obs, var = var, uns = uns_list, varm = as.list(
      SingleCellExperiment::reducedDims(sce)
    ),
    obsm = as.list(SingleCellExperiment::rowData(sce)$obsm), layers = assays_list
  )

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
#' @param matrix Matrix of single cell data (rows: genes, columns: cells), the row names of the
#'   matrix should be labeled with gene labels
#' @param cell_labels Annotation vector for celltypes (e.g columns of matrix)
#' @param named_metadata_list A named list containing unstructured metadata information (optional)
#'
#' @return SingleCellExperiment
#'
matrix_to_singlecellexperiment <- function(matrix, cell_labels, named_metadata_list = list()) {
  gene_labels <- rownames(matrix)
  sce <- SingleCellExperiment::SingleCellExperiment(list(X = matrix),
    colData = S4Vectors::DataFrame(label = cell_labels),
    rowData = S4Vectors::DataFrame(length = gene_labels),
    metadata = named_metadata_list
  )
  return(sce)
}

#' Convert SingleCellExperiment to Matrix + annotation vector for celltypes
#'
#' @param sce SingleCellExperiment
#' @param assay_name name of the assay that should be returned as matrix. (default: first in assay
#'   list)
#' @param cell_type_column_name name of the column that stores the cell-type labels in the colData
#'   object of the provided SingleCellExperiment
#'
#' @return named list: ..$matrix: matrix object, ..$annotation_vector
#'
singlecellexperiment_to_matrix <- function(sce, assay_name = NULL, cell_type_column_name = NULL) {
  if (is.null(assay_name)) {
    if (length(SummarizedExperiment::assays(sce)) == 0) {
      stop("'sce' does not contain any assays")
    }
    assay_name <- SummarizedExperiment::assayNames(sce)[1]
    message("Note: using the '", assay_name, "' assay as the X matrix")
  }
  X <- SummarizedExperiment::assay(sce, assay_name, withDimnames = T)

  if (is.null(cell_type_column_name)) {
    base::warning("Please provide the column name that contains the cell type labels! (colData
                  object of SingleCellExperiment/ obs object of AnnData)")
    cell_labels <- NULL
  }
  else {
    cell_labels <- SingleCellExperiment::colData(sce)[[cell_type_column_name]]
  }



  return(list("matrix" = X, "annotation_vector" = cell_labels))
}

#' Check if two SingleCellExperiments are identical
#'
#' @param a SingleCellExperiments
#' @param b SingleCellExperiments
#'
#' @return boolean
#'
sces_are_identical <- function(a, b) {
  same <- TRUE
  assays_a <- as.list(SummarizedExperiment::assays(a))
  assays_b <- as.list(SummarizedExperiment::assays(b))
  if (length(assays_a) != length(assays_b)) {
    same <- FALSE
    message("number of assays not identical")
  }
  else {
    if (!identical(rownames(a), rownames(b))) {
      same <- FALSE
      message("rownames not identical")
    }

    if (!identical(colnames(a), colnames(b))) {
      same <- FALSE
      message("rownames not identical")
    }

    a_names <- as.list(SummarizedExperiment::assayNames(a))
    b_names <- as.list(SummarizedExperiment::assayNames(b))
    for (i in 1:length(assays_a)) {
      m <- as.numeric(as.matrix(assays_a[[i]]))
      m2 <- as.numeric(as.matrix(assays_b[[i]]))

      if (!identical(m, m2)) {
        same <- FALSE

        message(paste("Assay", a_names[[i]], "and assay", b_names[[i]], " not identical"))
      }
    }

    a_col <- base::as.data.frame(SummarizedExperiment::colData(a))
    b_col <- base::as.data.frame(SummarizedExperiment::colData(b))

    if (!identical(a_col, b_col)) {
      same <- FALSE
      base::message("Coldata not identical")
    }

    a_row <- base::as.data.frame(SummarizedExperiment::rowData(a))
    b_row <- base::as.data.frame(SummarizedExperiment::rowData(b))

    if (!identical(a_row, b_row)) {
      same <- FALSE
      base::message("Rowdata not identical")
    }
  }
  return(same)
}

### AnnData methods

#' Create anndata object from single cell expression matrix, cell type labels and gene symbols
#'
#' @param x single cell expression matrix, rows = genes, columns = samples
#' @param obs vector of cell type labels
#' @param var vector of gene symbols
#' @param obsm Key-indexed multi-dimensional observations annotation of length #observations.
#' @param varm Key-indexed multi-dimensional variables annotation of length #variables.
#'
#' @return AnnData object
#'
build_anndata <- function(x, obs, var, obsm = NULL, varm = NULL) {
  anndata_check_load()

  x <- as.matrix(x)

  ad <- anndata::AnnData(X = x, obs = obs, var = var, obsm = obsm, varm = varm)

  rownames(ad) <- var
  colnames(ad) <- obs

  return(ad)
}

#' Read anndata object (.h5ad format)
#'
#' @param path path to .h5ad formatted file
#'
#' @return AnnData object
#'
read_anndata <- function(path) {
  anndata_check_load()
  data <- anndata::read_h5ad(path)
  return(data)
}

#' Write anndata object (.h5ad format)
#'
#' @param data AnnData object
#' @param path path where AnnData object should be written to (.h5ad format)
#'
write_anndata <- function(data, path) {
  anndata_check_load()
  data$write_h5ad(path)
}

#' Checks if anndata package is loaded
#'
#' If called and python environment is not set up, this is realized. Else, it checks if the anndata
#' package is loaded, and if not, it does this.
#'
anndata_check_load <- function() {
  if (!python_available()) {
    init_python()
    anndata_check_load()
  }

  if (!reticulate::py_module_available("anndata")) {
    anndata::install_anndata()
  }
}

anndata_is_identical <- function(a, b) {
  same <- TRUE

  if (!identical(as.matrix(a$X), as.matrix(b$X))) {
    base::message("X object is different")
    same <- FALSE
  }
  if (!identical(as.matrix(a$obs), as.matrix(b$obs))) {
    base::message("obs object is different")
    same <- FALSE
  }
  if (!identical(as.matrix(a$var), as.matrix(b$var))) {
    base::message("var object is different")
    same <- FALSE
  }
  if (!all.equal(a$uns, b$uns)) {
    base::message("uns object is different")
    same <- FALSE
  }
  if (!all.equal(a$obsm, b$obsm)) {
    base::message("obsm object is different")
    same <- FALSE
  }
  if (!all.equal(a$varm, b$varm)) {
    base::message("varm object is different")
    same <- FALSE
  }
  if (!all.equal(a$layers, b$layers)) {
    base::message("layers object is different")
    same <- FALSE
  }
  return(same)
}
