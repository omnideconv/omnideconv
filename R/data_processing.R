#' Generating a single cell expression set
#'
#' @param single_cell_matrix The single-cell matrix. Rows are genes, columns are samples.
#' @param sample_names A vector of the names of the samples, basically colnames(single_cell_matrix).
#' @param genes A vector of the names of the genes, basically rownames(single_cell_matrix).
#' @param cell_types A Vector of the cell type annotations. Has to be in the same order as the samples in single_cell_object.
#'
#' @return A Biobase::ExpressionSet of the input data.
#' @export
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
