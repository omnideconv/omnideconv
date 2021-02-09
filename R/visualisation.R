#' Plot deconvolution results as a barplot
#'
#' @param deconv_result result from deconvolution()
#' @param method_name (optional) title of plot is set to the name of the method used for deconvolution
#' @param file_name (optional) plot is saved in this file
#'
#' @return the ggplot object
#' @export
#'
#' @examples plotDeconvResult(immunedeconv2::deconvolute(bulk, immunedeconv2::build_model(single_cell_data, cell_type_annotations, "bisque"), "bisque", single_cell_data, cell_type_annotations), "Bisque")
plotDeconvResult <- function(deconv_result, method_name = "", file_name = NULL){

  plot <- data.table(deconv_result, samples= rownames(deconv_result)) %>%
    tidyr::pivot_longer(!samples, names_to ="cell_type", values_to="predicted_fraction") %>%
    ggplot2::ggplot(aes(y=samples, x=predicted_fraction, fill=cell_type))+geom_bar(stat="identity", position = "stack") +
    labs(title=method_name, y="sample", x="proportion", fill="cell type")

  if(ncol(deconv_result)<= 12){
    plot <- plot + scale_fill_brewer(palette="Paired")
  }

  if(!is.null(file_name)){
    ggsave(file_name, plot, width = 7, height = 4)
  }

  return(plot)
}


#' Make a Scatterplot for Benchmarking
#'
#' @param result_list named list containing all deconvolution results that should be considered,
#' cell type annotations need to contain "T cell", "DC", "Mono", "NK", "B cell" in any way to be included
#'
#' @param file_name (optional) plot is saved in this file
#'
#' @return the ggplot object
#' @export
#'
#' @examples makeBenchmarkingScatterplot(list(Bisque = result_bisque, MOMF = result_momf))
makeBenchmarkingScatterplot <- function(result_list, file_name){

  li <- lapply(result_list, function(x) cbind(x, T_cell = rowSums(x[,grepl("T cell", colnames(x)), drop=FALSE]),
                                              DC = rowSums(x[,grepl("DC", colnames(x)), drop=FALSE]),
                                              Monocyte = rowSums(x[,grepl("Mono", colnames(x)), drop=FALSE]),
                                              NK_cell = rowSums(x[,grepl("NK", colnames(x)), drop=FALSE]),
                                              B_cell = rowSums(x[,grepl("B cell", colnames(x)), drop=FALSE])))
  li <- lapply(li, function(x) cbind(data.table(x), sample = rownames(x)))
  li <- lapply(names(li), function(i) cbind(li[[i]], method = rep(i, nrow(li[[i]]))))

  names(li) <- names(result_list)
  li <- lapply(li, function(x) melt(x, id.vars = c("sample", "method"), variable.name = "cell_type", value.name = "predicted_fraction"))
  df <- dplyr::bind_rows(li)
  df$cell_type <- gsub(" ", "_", df$cell_type)
  load("data/RefData.RData")
  names(RefData) <- c("T_cell", "Monocyte", "B_cell", "DC", "NK_cell")
  RefData$sample <- rownames(RefData)
  plot <- tidyr::pivot_longer(RefData, !sample, names_to="cell_type", values_to="true_fraction") %>% merge(df, by=c("sample", "cell_type")) %>%
    ggplot2::ggplot(aes(x=as.numeric(true_fraction), y=predicted_fraction, color=cell_type))+geom_point(size=4)+facet_wrap(~method)+
    geom_abline(color="black")+scale_y_continuous(breaks=seq(0, 1, 0.25))+scale_x_continuous(breaks=seq(0, 1, 0.25))+
    coord_cartesian(xlim = c(0, 0.75), ylim = c(0, 0.75))+labs(x="true fraction", y="predicted fraction", color="cell type")+
    theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 13), axis.text =
            element_text(size = 12), axis.title = element_text(size = 13))+
    scale_color_manual(values = c("deepskyblue", "springgreen3", "palevioletred1", "red", "blue"), breaks=c("B_cell","Monocyte", "DC","NK_cell","T_cell"),
                       labels=c("B cell","Monocyte", "DC","NK cell","T cell"))
  if(!is.null(file_name)){
    ggplot2::ggsave(file_name, plot, width = 6, height = 5)
  }
  return(plot)
}


