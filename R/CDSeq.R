#' No model is build as CDSeq does both steps in one.
#'
#' Please use the deconvolute method with your single cell and bulk rna seq data to use CDSeq.
#'
#'
#' @return NULL
#'
#' @export
build_model_cdseq <- function() {
  base::message(
    "The deconvolution with CDSeq is done in only one step. Please just use the ",
    "deconvolute method."
  )

  return(NULL)
}

#' CDSeq Deconvolution
#'
#' This function is to calculate the CDSeq deconvolution proportions.
#' IMPORTANT: No model is needed. Everything is done inside this method.
#' IMPORTANT: The result does not necessarily contain all cell types from the input single cell
#' data. It assigns cell types to clusters found in the bulk data. See CDSeq_cell_type_assignment_df
#' for more information.
#'
#'
#' @param bulk_gene_expression  A matrix or dataframe with the bulk data. Rows are genes, columns
#'   are samples.
#' @param single_cell_object A Matrix with the single-cell data. Rows are genes and columns are
#'   samples.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param batch_ids A vector of the ids of the samples or individuals.
#' @param cdseq_gep_sample_specific CDSeq-estimated sample-specific cell type gene expression, in
#'   the form of read counts. It is a 3 dimension array, i.e. gene by sample by cell type. The
#'   element cdseq_gep_sample_specific\[i,j,k\] represents the reads mapped to gene i from cell type
#'   k in sample j.
#' @param batch_correction perform Harmony batch correction if it is 1.
#' @param harmony_iter Maximum number of rounds to run Harmony. One round of Harmony involves one
#'   clustering and one correction step.
#' @param harmony_cluster Maximum number of rounds to run clustering at each round of Harmony.
#' @param nb_size size parameter for negative binomial distribution, check rnbinom for details.
#' @param nb_mu mu parameter for negative binomial distribution, check rnbinom for details.
#' @param breaksList parameter for pheatmap controling the color scale. See pheatmap function for
#'   details.
#' @param corr_threshold if the correlation between CDSeq-estimated GEPs and the scRNAseq GEP is
#'   below this value, then it is considered the two cell types are not matching.
#' @param pseudo_cell_count an integer indicating how many pseudo cells will be generated from
#'   CDSeq-estimated cell-type-specific gene expression profiles. Default values is 1.
#' @param seurat_count_threshold this parameter will be passed to Seurat subset function
#'   (subset = nCount_RNA > seurat_count_threshold) for filtering out single cells whose total
#'   counts is less this threshold.
#' @param seurat_scale_factor this parameter will be passed to scale.factor in Seurat function
#'   NormalizeData.
#' @param seurat_norm_method this parameter will be passed to normalization.method in Seurat
#'   function NormalizeData.
#' @param seurat_select_method this parameter will be passed to selection.method in Seurat
#'   function FindVariableFeatures
#' @param seurat_nfeatures this parameter will be passed to nfeatures in Seurat function
#'   FindVariableFeatures.
#' @param seurat_npcs this parameter will be passed to npcs in Seurat function RunPCA.
#' @param seurat_dims this parameter will be passed to dims in Seurat function FindNeighbors.
#' @param seurat_reduction this parameter will be passed to reduction in Seurat function
#'   FindNeighbors.
#' @param seurat_resolution this parameter will be passed to resolution in Seurat function
#'   FindClusters.
#' @param seurat_find_marker this parameter controls if run seurat FindMarker function, default
#'   is FALSE.
#' @param seurat_DE_test this parameter will be passed to test.use in Seurat function
#'   FindAllMarkers.
#' @param seurat_DE_logfc this parameter will be passed to logfc.threshold in Seurat function
#'   FindAllMarkers.
#' @param seurat_top_n_markers the number of top DE markers saved from Seurat output.
#' @param sc_pt_size point size of single cell data in umap and tsne plots
#' @param cdseq_pt_size point size of CDSeq-estimated cell types in umap and tsne plots
#' @param plot_umap set 1 to plot umap figure of scRNAseq and CDSeq-estimated cell types,
#'   0 otherwise.
#' @param plot_tsne set 1 to plot tsne figure of scRNAseq and CDSeq-estimated cell types,
#'   0 otherwise.
#' @param plot_per_sample currently disabled for debugging
#' @param fig_save 1 or 0. 1 means save figures to local and 0 means do not save figures to local.
#' @param fig_path the location where the heatmap figure is saved.
#' @param fig_name the name of umap and tsne figures. Umap figure will have the name of
#'   fig_name_umap_date and tsne figure will be named fig_name_tsne_date.
#' @param fig_format "pdf", "jpeg", or "png".
#' @param corr_heatmap_fontsize font size of the correlation heatmap between scRNAseq GEP and
#'   CDSeq-estimated GEPs.
#' @param fig_dpi figure dpi
#' @param verbose Whether to produce an output on the console
#' @return A list including:
#' \item{fig_path}{The same as the input fig_path.}
#' \item{fig_name}{The same as the input fig_name.}
#' \item{cdseq_synth_scRNA}{The synthetic scRNAseq data generated using CDSeq-estiamted GEPs.}
#' \item{cdseq_scRNA_umap}{A ggplot figure of the umap outcome.}
#' \item{cdseq_scRNA_tsne}{A ggplot figure of the tsne outcome.}
#' \item{cdseq_synth_scRNA_seurat}{A Seurat object containing the scRNAseq combined with
#'   CDSeq-estimated cell types. Cell id for CDSeq-estimated cell types start with "CDSeq".}
#' \item{seurat_cluster_purity}{For all cells in a Seurat cluster i,  the ith value in
#'   seurat_cluster_purity is the proportion of the mostly repeated cell annotation from
#'   sc_annotation.
#'   For example, after Seurat clustering, suppose there are 100 cells in cluster 1, out of these
#'   100 cells, 90 cells' annotation in sc_annotation is cell type A, then the fist value in
#'   seurat_cluster_purity is 0.9. This output can be used to assess the agreement between
#'   Seurat clustering and the given sc_annotation.}
#' \item{seurat_unique_clusters}{Unique Seurat cluster numbering. This can be used together with
#'   seurat_cluster_gold_label to match the Seurat clusters with given annotations.}
#' \item{seurat_cluster_gold_label}{The cell type annotations for each unique Seurat cluster based
#'   on sc_annotation.}
#' \item{seurat_markers}{DE genes for each Seurat cluster.}
#' \item{seurat_top_markers}{Top seurat_top_n_markers DE genes for each Seurat cluster.}
#' \item{CDSeq_cell_type_assignment_df}{The cell type assignment for CDSeq-estimated cell types.}
#' \item{CDSeq_cell_type_assignment_confidence}{The cell type assignment confidence matrix, only
#'   available when pseudo_cell_count > 1.}
#' \item{CDSeq_cell_type_assignment_df_all}{The cell type assignment for CDSeq-estimated cell types,
#'   only available when pseudo_cell_count > 1.}
#' \item{cdseq_prop_merged}{CDSeq-estimated cell type proportions with cell type annotations
#'   (annotated using clustering with scRNAseq).}
#' \item{cdseq_gep_sample_specific_merged}{Sample-specific cell-type read counts. It is a 3d array
#'   with dimensions: gene, sample, cell type.}
#' \item{input_list}{The values of the input parameters.}
#' \item{cdseq_sc_comb_umap_df}{The dataframe for umap plot.}
#' \item{cdseq_sc_comb_tsne_df}{The dataframe for tsne plot.}
#' \item{cdseq_prop_merged_byCorr}{CDSeq-estimated cell type proportions with cell type annotations
#'   (annotated using correlation with scRNAseq).}
#' \item{cdseq_gep_merged_byCorr}{CDSeq-estimated cell-type-specific GEPs with cell type annotations
#'   (annotated using correlation with scRNAseq).}
#' \item{cdseq_annotation_byCorr}{CDSeq-estimated cell type annotations (annotated using correlation
#'   with scRNAseq).}
#' @export
deconvolute_cdseq <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                              batch_ids, cdseq_gep_sample_specific = NULL, batch_correction = 1,
                              harmony_iter = 10, harmony_cluster = 20, nb_size = NULL, nb_mu = NULL,
                              corr_threshold = 0, breaksList = seq(0, 1, 0.01), pseudo_cell_count = 1,
                              seurat_count_threshold = 0, seurat_scale_factor = 10000,
                              seurat_norm_method = "LogNormalize", seurat_select_method = "vst",
                              seurat_nfeatures = 1000, seurat_npcs = 30, seurat_dims = 1:30,
                              seurat_reduction = "pca", seurat_resolution = 0.8,
                              seurat_find_marker = FALSE, seurat_DE_test = "wilcox",
                              seurat_DE_logfc = 0.25, seurat_top_n_markers = 10, sc_pt_size = 1,
                              cdseq_pt_size = 3, plot_umap = 0, plot_tsne = 0, plot_per_sample = 0,
                              fig_save = 0, fig_path = getwd(),
                              fig_name = "CDSeqCellTypeAssignSCRNA", fig_format = "jpeg",
                              fig_dpi = 100, corr_heatmap_fontsize = 10, verbose = FALSE) {
  if (is.null(single_cell_object) || is.null(cell_type_annotations) || is.null(batch_ids)) {
    base::stop(
      "Single cell object, cell type annotations or batch_ids not provided. Call as: ",
      "deconvolute(bulk_gene_expression, NULL, \"cdseq\", single_cell_object, ",
      "cell_type_annotations, batch_ids)"
    )
  }

  cdseq_res <- CDSeq::CDSeq(
    bulk_data = bulk_gene_expression,
    cell_type_number = length(unique(cell_type_annotations))
  )

  cdseq_gep <- cdseq_res$estGEP
  cdseq_prop <- cdseq_res$estProp

  cdseq_gep <- cdseq_gep[base::intersect(
    rownames(single_cell_object),
    rownames(cdseq_gep)
  ), ]

  single_cell_object <- single_cell_object[rownames(cdseq_gep), ]

  cell_annotations_as_df <- cbind(colnames(single_cell_object), cell_type_annotations)
  colnames(cell_annotations_as_df) <- c("cell_id", "cell_type")
  cell_annotations_as_df <- data.frame(cell_annotations_as_df)
  return(CDSeq::cellTypeAssignSCRNA(
    cdseq_gep = cdseq_gep,
    cdseq_prop = cdseq_prop,
    sc_gep = single_cell_object,
    sc_annotation = cell_annotations_as_df,
    sc_batch = batch_ids,
    cdseq_gep_sample_specific = cdseq_gep_sample_specific,
    batch_correction = batch_correction,
    harmony_iter = harmony_iter,
    harmony_cluster = harmony_cluster,
    nb_size = nb_size,
    nb_mu = nb_mu,
    corr_threshold = corr_threshold,
    breaksList = breaksList,
    pseudo_cell_count = pseudo_cell_count,
    seurat_count_threshold = seurat_count_threshold,
    seurat_scale_factor = seurat_scale_factor,
    seurat_norm_method = seurat_norm_method,
    seurat_select_method = seurat_select_method,
    seurat_nfeatures = seurat_nfeatures,
    seurat_npcs = seurat_npcs,
    seurat_dims = seurat_dims,
    seurat_reduction = seurat_reduction,
    seurat_resolution = seurat_resolution,
    seurat_find_marker = seurat_find_marker,
    seurat_DE_test = seurat_DE_test,
    seurat_DE_logfc = seurat_DE_logfc,
    seurat_top_n_markers = seurat_top_n_markers,
    sc_pt_size = sc_pt_size,
    cdseq_pt_size = cdseq_pt_size,
    plot_umap = plot_umap,
    plot_tsne = plot_tsne,
    plot_per_sample = plot_per_sample,
    fig_save = fig_save,
    fig_path = fig_path,
    fig_name = fig_name,
    fig_format = fig_format,
    fig_dpi = fig_dpi,
    corr_heatmap_fontsize = corr_heatmap_fontsize,
    verbose = verbose
  ))
}
