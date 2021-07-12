#' No model is build as CPM does not use a signature matrix.
#'
#' Please use the deconvolute method with you single cell and bulk rna seq data to use CPM.
#'
#'
#' @return NULL
#'
#' @export
build_model_cpm <- function() {
  base::message(
    "The deconvolution with CPM is done in only one step. Please just use the",
    "deconvolute method."
  )

  return(NULL)
}

#' CPM Deconvolution
#'
#' This function is to calculate the CPM deconvolution proportions.
#' IMPORTANT: No model is needed. Everything is done inside this method.
#' This method is NOT deterministic, so if it is run multiple times, it will create different
#' outputs.
#'
#' This function initiate the Cellular Population Mapping (CPM) algorithm - a deconvolution
#' algorithm in which single-cell genomics is required in only one or a few samples, where in other
#' samples of the same tissue, only bulk genomics is measured and the underlying fine resolution
#' cellular heterogeneity is inferred.
#' CPM predicts the abundance of cells (and cell types) ranging monotonically from negative to
#' positive levels. Using a relative framework these values correspond to decrease and increase in
#' cell abundance levels, respectively. On the other hand, in an absolute framework lower values
#' (including negatives) correspond to lower abundances and vise versa. These values are comparable
#' between samples.
#'
#'
#' @param bulk_gene_expression  A matrix or dataframe with the bulk data. Rows
#'   are genes, columns are samples.
#' @param single_cell_object A Matrix with the single-cell data. Rows are genes
#'  and columns are samples.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to
#'  be in the same order as the samples in single_cell_object.
#' @param cell_space The cell state space corresponding to the single-cell RNA-seq data. It can be
#'  a vector for a 1-dim space or a 2D matrix for a two space where each column represents a
#'  different dimension. The cell space should incorporate the similarities of cells within cell
#'  types. Similarities between cells from different cell types, based on the cell space, are not
#'  taken into account in CPM.
#'  It is also possible to supply the string "PCA", "UMAP" or "TSNE" which calculates the cell
#'  space using the corresponding method (using the Seurat implementation and default parameters).
#' @param no_cores A number for the amount of cores which will be used for the analysis. The
#'   default (NULL) is total number of cores minus 1.
#' @param neighborhood_size Cell neighborhood size which will be used for the analysis. This should
#'   be lower than the number of cells in the smallest cell type. The default is 10.
#' @param model_size The reference subset size in each iteration of CPM. This should be lower than
#'   the total number of cells. The default is 50.
#' @param min_selection The minimum number of times in which each reference cell is selected.
#'   Increasing this value might have a large effect on the algorithm's running time.
#'   The default is 5.
#' @param calculate_CI A boolean parameter indicating whether the calculation of confidence
#'   intervals is needed. The default is FALSE.
#' @param verbose Whether the algorithm should print out what it is doing.
#' @return A list including:
#' \item{predicted}{CPM predicted cell abundance matrix. Each row represents a sample and
#'   each column a single cell.}
#' \item{cellTypePredictions}{CPM predicted cell-type abundance matrix. Each row represents a sample
#'   and each column a single cell-type.}
#' \item{confIntervals}{A matrix containing the confidence interval for each cell and sample. Each
#'   row represents a sample and each column a single cell. This is calculated
#'   if calculate_CI = TRUE.}
#' \item{numOfRuns}{The number of deconvolution repeats preformed by CPM. }
#' @export
deconvolute_cpm <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                            cell_space = "PCA", no_cores = NULL, neighborhood_size = 10,
                            model_size = 50, min_selection = 5, calculate_CI = FALSE,
                            verbose = FALSE) {
  if (is.null(single_cell_object) || is.null(cell_type_annotations)) {
    base::stop(
      "Single cell object or cell type annotations not provided. Call as: ",
      "deconvolute(bulk_gene_expression, NULL, \"cpm\", single_cell_object, ",
      "cell_type_annotations)"
    )
  }
  if ("character" %in% class(cell_space) && length(cell_space) == 1) {
    cell_space <- calculate_cell_embedding(single_cell_object, cell_type_annotations, cell_space)
  }


  return(scBio::CPM(single_cell_object, cell_type_annotations, bulk_gene_expression, cell_space,
    no_cores = no_cores, neighborhoodSize = neighborhood_size, modelSize = model_size,
    minSelection = min_selection, quantifyTypes = TRUE, typeTransformation = TRUE,
    calculateCI = calculate_CI
  ))
}




#' Calculation of the cell_space parameter needed by CPM
#'
#' @param single_cell_object A Matrix with the single-cell data. Rows are genes
#'  and columns are samples.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to
#'  be in the same order as the samples in single_cell_object.
#' @param method Either "PCA", "UMAP" or "TSNE"
#'
#' @return A matrix with two dimensions. The rows are the cells, the columns are the two dimensions
#'  calculated by Seurat
calculate_cell_embedding <- function(single_cell_object, cell_type_annotations,
                                     method = c("PCA", "UMAP", "TSNE")) {
  if (length(method) > 1) {
    method <- method[[1]]
  }
  method <- base::tolower(method)
  sce <- matrix_to_singlecellexperiment(single_cell_object, cell_type_annotations)
  seurat <- Seurat::as.Seurat(sce, counts = "X", data = NULL)

  seurat <- Seurat::NormalizeData(seurat)
  all.genes <- rownames(seurat)
  seurat <- Seurat::ScaleData(seurat, features = all.genes)
  seurat <- Seurat::FindVariableFeatures(seurat)
  seurat <- Seurat::RunPCA(seurat, features = Seurat::VariableFeatures(object = seurat))
  if (method == "pca") {
    return(seurat@reductions$pca@cell.embeddings[, 1:2])
  } else if (method == "umap") {
    seurat <- Seurat::RunUMAP(seurat, dims = 1:10)
    return(seurat@reductions$umap@cell.embeddings)
  } else if (method == "tsne") {
    seurat <- Seurat::RunTSNE(seurat)
    return(seurat@reductions$tsne@cell.embeddings)
  }
  base::stop("Method ", method, " not recognized")
}
