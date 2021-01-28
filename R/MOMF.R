
#' Deconvolution Analysis using MOMF (via Nonnegative Factorization)
#'
#' @param bulk_gene_expression Dataframe or matrix of bulk RNA-seq data (genes x individuals)
#' @param signature Signature Matrix (genes x individuals from scRNA-seq)
#' @param single_cell_object scRNA-seq Object (genes x cells)
#' @param method Determines which divergence to use. Options: Kullback-Leibler "KL", Itakura-Saito "IS". Defaults to "KL"
#' @param verbose Whether the algorithm should print out what it is doing.
#' @param ... additional parameters
#'
#' @return cell proportion matrix
#' @export
#'
#' @examples
deconvolute_MOMF <- function(bulk_gene_expression, signature, single_cell_object, verbose = T, method="KL",...){
  #MOMF needs a list of the single_cell_object with cells x genes and the bulk RNA seq data with individuals x genes
  GList <- list(X1 = t(single_cell_object), X2 = t(bulk_gene_expression))
  if (!verbose){
    sink(tempfile())
    result <- tryCatch(MOMF::momf.fit(DataX = GList, DataPriorU=signature, method=method, ...),
                       finally = sink())
  } else {
    result <- MOMF::momf.fit(DataX = GList, DataPriorU=signature, method=method, ...)
  }
  if (is.null(result$cell.prop)){
    base::stop("Something went wrong. Please switch on verbose mode")
  }
  #return slot in result with cell proportion matrix
  return(result$cell.prop)
}
