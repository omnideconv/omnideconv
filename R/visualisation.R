#' Plot deconvolution results as a barplot
#'
#' @param deconv_result result from deconvolution()
#' @param method_name (optional) title of plot is set to the name of the method used for
#'   deconvolution
#' @param file_name (optional) plot is saved in this file
#' @import dplyr
#' @import ggplot2
#' @return the ggplot object
#'
#' @examples
#' data("single_cell_data_1")
#' data("cell_type_annotations_1")
#' data("batch_ids_1")
#' data("bulk")
#' data("RefData")
#'
#' common_genes <- intersect(rownames(single_cell_data_1), rownames(bulk))[1:2000]
#'
#' single_cell_data <- single_cell_data_1[common_genes, 1:500]
#' cell_type_annotations <- cell_type_annotations_1[1:500]
#' batch_ids <- batch_ids_1[1:500]
#' bulk <- bulk[common_genes, ]
#'
#' model <- build_model(single_cell_data, cell_type_annotations, "bisque", batch_ids)
#' deconvolution <- deconvolute(
#'   bulk, model, "bisque", single_cell_data,
#'   cell_type_annotations, batch_ids
#' )
#' plot_deconv_result(deconvolution, "Bisque")
#' @export
plot_deconv_result <- function(deconv_result, method_name = "", file_name = NULL) {
  plot <- cbind(deconv_result, samples = rownames(deconv_result)) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(!samples, names_to = "cell_type", values_to = "predicted_fraction") %>%
    ggplot2::ggplot(aes(y = samples, x = as.numeric(predicted_fraction), fill = cell_type)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(title = method_name, y = "sample", x = "proportion", fill = "cell type")

  if (ncol(deconv_result) <= 12) {
    plot <- plot + ggplot2::scale_fill_brewer(palette = "Paired")
  }

  if (!is.null(file_name)) {
    ggplot2::ggsave(file_name, plot, width = 7, height = 4)
  }

  return(plot)
}


#' Make a Scatterplot for Benchmarking
#'
#' @param result_list named list containing all deconvolution results that should be considered,
#'   cell type annotations need to contain the same cell types as the ones in ref_data
#' @param ref_data reference cell types which are used as the ground truth
#' @param file_name (optional) plot is saved in this file
#' @import dplyr
#' @import ggplot2
#' @return the ggplot object
#' @export
#'
#' @examples
#' data("single_cell_data_1")
#' data("cell_type_annotations_1")
#' data("batch_ids_1")
#' data("bulk")
#' data("RefData")
#'
#' common_genes <- intersect(rownames(single_cell_data_1), rownames(bulk))[1:2000]
#'
#' single_cell_data <- single_cell_data_1[common_genes, 1:500]
#' cell_type_annotations <- cell_type_annotations_1[1:500]
#' batch_ids <- batch_ids_1[1:500]
#' bulk <- bulk[common_genes, ]
#'
#' RefData <- RefData[, order(colnames(RefData))]
#'
#' sig_bisque <- build_model(
#'   single_cell_data, cell_type_annotations, "bisque",
#'   batch_ids
#' )
#' res_bisque <- deconvolute(
#'   bulk, sig_bisque, "bisque", single_cell_data,
#'   cell_type_annotations, batch_ids
#' )
#'
#' res_scdc <- deconvolute(bulk, NULL, "scdc", batch_ids,
#'   single_cell_object = single_cell_data,
#'   cell_type_annotations = cell_type_annotations
#' )
#'
#' result_list <- list(SCDC = res_scdc, Bisque = res_bisque)
#'
#' # Merging the two T cell props
#' result_list <- lapply(result_list, function(x) {
#'   cbind(x, T = (x[, "CD4 T"] + x[, "CD8 T"]))[, -c(2, 3)]
#' })
#' make_benchmarking_scatterplot(result_list, RefData)
#' # Alternative if you want to save the plot in a file
#' # make_benchmarking_scatterplot(result_list, RefData, "predictionVsGroundtruth.png")
make_benchmarking_scatterplot <- function(result_list, ref_data, file_name = NULL) {
  cell_types <- sort(unique(colnames(result_list[[1]])))
  if (length(colnames(ref_data)) != length(cell_types) ||
    !identical(sort(colnames(ref_data)), cell_types)) {
    stop("Reference and prediction need to include the same cell types")
  }

  li <- lapply(result_list, function(x) cbind(x, sample = rownames(x)))
  li <- lapply(names(li), function(i) cbind(li[[i]], method = rep(i, nrow(li[[i]]))))

  li <- lapply(li, function(x) {
    tidyr::pivot_longer(data.frame(x), !c("sample", "method"),
      names_to = "cell_type", values_to = "predicted_fraction"
    )
  })
  df <- dplyr::bind_rows(li)
  ref_data$sample <- rownames(ref_data)
  plot <- tidyr::pivot_longer(ref_data, !sample,
    names_to = "cell_type",
    values_to = "true_fraction"
  ) %>%
    merge(df, by = c("sample", "cell_type")) %>%
    ggplot2::ggplot(aes(
      x = as.numeric(true_fraction), y = as.numeric(predicted_fraction),
      color = cell_type
    )) +
    geom_point(size = 4) +
    facet_wrap(~method) +
    geom_abline(color = "black") +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "true fraction", y = "predicted fraction", color = "cell type") +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      axis.text = element_text(size = 12), axis.title = element_text(size = 13)
    ) +
    scale_color_manual(
      values = RColorBrewer::brewer.pal(length(cell_types), "Accent"),
      breaks = cell_types,
      labels = cell_types
    )
  if (!is.null(file_name)) {
    ggplot2::ggsave(file_name, plot, width = 6, height = 5)
  }
  return(plot)
}
