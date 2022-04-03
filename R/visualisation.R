#' Plot deconvolution results as a barplot
#'
#' @param result_list A named list containing all deconvolution results
#' @param title (optional) title of the plot
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
#' deconvolution <- deconvolute(
#'   bulk, NULL, "bisque", single_cell_data,
#'   cell_type_annotations, batch_ids
#' )
#' deconvolution <- list(deconvolution)
#' names(deconvolution) <- "Bisque"
#' make_barplot(deconvolution)
#' @export
make_barplot <- function(result_list, title = "", file_name = NULL) {
  if (!"list" %in% class(result_list)) {
    stop("Please supply a list for the parameter 'result_list'")
  }
  if (is.null(names(result_list))) {
    stop("Please supply a NAMED list, names(result_list) returns NULL")
  }
  if (length(unique(sapply(result_list, nrow))) > 1) {
    stop("Please supply a list where every entry contains the same number of samples")
  }
  plots <- lapply(result_list, function(result) {
    cbind(result, samples = rownames(result)) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(!samples, names_to = "cell_type", values_to = "predicted_fraction")
  })
  all <- do.call("rbind", plots)
  all$method <- rep(names(result_list), each = nrow(all) / length(names(result_list)))

  plot <- all %>%
    ggplot2::ggplot(aes(
      y = samples, x = as.numeric(predicted_fraction), fill = cell_type,
      group = method
    )) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::facet_wrap(~method) +
    ggplot2::labs(title = title, y = "sample", x = "proportion", fill = "cell type") +
    ggplot2::theme(legend.position = "bottom")

  if (ncol(result_list[[1]]) <= 12) {
    plot <- plot + ggplot2::scale_fill_brewer(palette = "Paired")
  }

  if (!is.null(file_name)) {
    ggplot2::ggsave(file_name, plot, width = 7, height = 4)
  }

  return(plot)
}


#' Make a Scatterplot for Benchmarking
#'
#' @param result_list A named list containing all deconvolution results that should be considered,
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
#' res_bisque <- deconvolute(
#'   bulk, NULL, "bisque", single_cell_data,
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
  if (!"list" %in% class(result_list)) {
    stop("Please supply a list for the parameter 'result_list'")
  }
  if (is.null(names(result_list))) {
    stop("Please supply a NAMED list, names(result_list) returns NULL")
  }
  if (length(unique(sapply(result_list, nrow))) > 1) {
    stop("Please supply a list where every entry contains the same number of samples")
  }
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


#' Plot Deconvolution results
#'
#' @param deconvolutions A named list of deconvolution results
#' @param plot_method Type of plot to be rendered  ("bar", "jitter", "scatter", "box", "heatmap")
#' @param facet Variable for grouping the plots ("method", "cell_type", "sample")
#' @param palette RColorBrewer palette name (optional), standard = "Set1"
#' @import ggplot2
#' @import tidyr
#' @import magrittr
#' @import RColorBrewer
#' @returns ggplot rendered by plotly for interactivity
#' @export
#'
#' @examples
#' data("single_cell_data_1")
#' data("cell_type_annotations_1")
#' data("batch_ids_1")
#' data("bulk")
#'
#' deconvolution <- omnideconv::deconvolute(
#'   bulk, NULL, "bisque",
#'   single_cell_data_1,
#'   cell_type_annotations_1,
#'   batch_ids_1
#' )
#' deconvolution <- list(deconvolution)
#' names(deconvolution) <- "bisque"
#' omnideconv::plot_deconvolution(deconvolution, "bar", "method", "Spectral")
plot_deconvolution <- function(deconvolutions, plot_method = "bar", facet = "method", palette = "Spectral") {
  # data needs to be a named deconvolution list
  if (is.null(names(deconvolutions))) {
    stop("Please supply a NAMED list, names(deconvolutions) returns NULL")
  }

  # check plot Method
  if (!(plot_method %in% c("bar", "jitter", "scatter", "box", "heatmap"))) {
    stop("plot_method not supported. Please provide one of the following: 'bar', 'jitter', 'scatter', 'box', 'heatmap'")
  }

  # check facet parameter
  if (!(facet %in% c("method", "cell_type", "sample"))) {
    stop("facet not supported. Please provide one of the following: 'method', 'cell_type', 'sample'")
  }

  # check palette name
  if (!(palette %in% rownames(RColorBrewer::brewer.pal.info))) {
    stop("palette not a RColorBrewer palette name")
  }

  # preformat data into a dataframe
  deconvolutions <- lapply(deconvolutions, function(deconvolution) {
    cbind(deconvolution, sample = rownames(deconvolution)) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(!sample, names_to = "cell_type", values_to = "fraction")
  })

  # combine list to one dataframe and add calculation method
  data <- do.call("rbind", deconvolutions)
  data$fraction <- as.numeric(data$fraction) # change datatype of column
  data$method <- rep(names(deconvolutions), each = nrow(data) / length(names(deconvolutions))) # add computation method as column

  # calculate tooltip based on chosen facet
  tooltip <- switch(facet,
    "method" = aes(
      text = paste0("Cell Type: ", cell_type, "\nFraction: ", sprintf("%1.2f%%", 100 * fraction), "\nSample: ", sample)
    ),
    "cell_type" = aes(
      text = paste0("Sample: ", sample, "\nFraction: ", sprintf("%1.2f%%", 100 * fraction), "\nMethod: ", method)
    ),
    "sample" = aes(
      text = paste0("Cell Type: ", cell_type, "\nFraction: ", sprintf("%1.2f%%", 100 * fraction), "\nMethod: ", method)
    )
  )

  # axis information
  axis <- list(
    "method" = list( # x, y, fill/color
      "bar" = list("fraction", "sample", "cell_type"),
      "jitter" = list("fraction", "cell_type", "cell_type"),
      "scatter" = list("fraction", "cell_type", "cell_type"),
      "box" = list("cell_type", "fraction", "cell_type"),
      "heatmap" = list("cell_type", "sample", "fraction")
    ),
    "cell_type" = list(
      "bar" = list("fraction", "sample", "method"),
      "jitter" = list("fraction", "method", "sample"),
      "scatter" = list("fraction", "method", "sample"),
      "box" = list("method", "fraction", "method"),
      "heatmap" = list("sample", "method", "fraction")
    ),
    "sample" = list(
      "bar" = list("fraction", "cell_type", "method"),
      "jitter" = list("fraction", "cell_type", "method"),
      "scatter" = list("fraction", "cell_type", "method"),
      "box" = list("method", "fraction", "method"),
      "heatmap" = list("cell_type", "method", "fraction")
    )
  )

  # extract ggplot aes
  x <- axis[[facet]][[plot_method]][[1]] # x axis
  y <- axis[[facet]][[plot_method]][[2]] # y axis
  col <- axis[[facet]][[plot_method]][[3]] # color / fill, depending on plot type

  aes <- NULL

  if (plot_method %in% c("jitter", "scatter")) {
    aes <- aes_string(x = x, y = y, col = col)
  } else {
    aes <- aes_string(x = x, y = y, fill = col)
  }

  # construct plot
  plot <- ggplot(data, aes)
  plot <- plot + facet_wrap(~ data[[facet]])

  # general theme
  plot <- plot + bbplot::bbc_style() +
    theme(
      legend.title = element_text(size = 16), # legend title font size
      legend.text = element_text(size = 14), # legend element font size
      axis.text.x = element_text(size = 14), # x axis font size
      axis.text.y = element_text(size = 14), # y axis font size
      axis.ticks.x = ggplot2::element_line(colour = "#333333"), # vertical ticks for fractions
      axis.ticks.length = grid::unit(0.26, "cm"), # tick length in cm
      strip.text = element_text(size = 16)
    )

  # add plot content
  if (plot_method == "bar") {
    if (facet != "method") {
      plot <- plot + geom_col(tooltip, position = "dodge") # not stacked
    } else {
      plot <- plot + geom_col(tooltip) +
        theme(panel.grid.major.y = ggplot2::element_blank()) # remove vertical lines
    }
  } else if (plot_method == "jitter") {
    plot <- plot + geom_jitter(tooltip)
  } else if (plot_method == "scatter") {
    plot <- plot + geom_point(tooltip)
  } else if (plot_method == "box") {
    plot <- plot + geom_boxplot(tooltip) +
      coord_flip() # this is mandatory here
  } else if (plot_method == "heatmap") {
    plot <- plot + geom_tile(tooltip) +
      theme(axis.text.x = element_text(angle = 90)) +
      guides(fill = guide_colorbar(barwith = 0.5, barheight = 20))
  }

  # update palette to match number of celltypes
  max_colors <- RColorBrewer::brewer.pal.info[palette, ]$maxcolors # for brewer.pal()
  n_cell_types <- length(unique(data$cell_type)) # number of needed colors
  getPalette <- colorRampPalette(brewer.pal(max_colors, palette)) # function to return custom interpolated palettes

  # add color theme based on plot method
  if (plot_method %in% c("jitter", "scatter")) {
    plot <- plot + ggplot2::scale_colour_manual(values = getPalette(n_cell_types))
  } else if (plot_method == "heatmap") {
    scale_fill_gradient(low = "#FFFFFF", high = RColorBrewer::brewer.pal(max_colors, palette)[1:1])
  } else {
    plot <- plot + ggplot2::scale_fill_manual(values = getPalette(n_cell_types))
  }

  # render
  plotly::ggplotly(plot, tooltip = c("text")) %>%
    plotly::config(
      displaylogo = FALSE, showTips = FALSE, toImageButtonOptions = list(filename = paste0(plot_method, "_plot")),
      modeBarButtonsToRemove = list(
        "hoverCLosestCartesian",
        "hoverCompareCartesian",
        "zoomIn2d", "zoomOut2d",
        "lasso2d", "zoom2d",
        "pan2d", "autoScale2d", "select2d"
      )
    ) %>%
    plotly::layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE))
}
