#' Example RNA-seq dataset from Hoek et al.
#'
#' A dataset of eight RNA-seq samples from patients vaccinated with trivalent inactivated influenza
#' vaccine (TIV)
#'
#' @format matrix with genes x samples
#'
#' @source Hoek et al. (2015), PLOS ONE,
#'   https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118528
"bulk"

#' Small example RNA-seq dataset from Hoek et al.
#'
#' A dataset of eight RNA-seq samples from patients vaccinated with trivalent inactivated influenza
#' vaccine (TIV) including only 100 genes
#'
#' @format matrix with genes x samples
#'
#' @source Hoek et al. (2015), PLOS ONE,
#'   https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118528
"bulk_small"

#' Gold standard measurements with FACS from Hoek et al.
#'
#' A dataset of cell-type fractions of eight RNA-seq samples measured with FACS
#'
#' @format matrix with samples x cell types
#'
#' @source Hoek et al. (2015), PLOS ONE,
#'   https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118528
"RefData"

#' Example sigle-cell RNA-seq dataset from Maynard et al.
#'
#' A subset of TPM-normalized single-cell RNA-seq gene expression, containing 300 cells with types
#' "T cell CD4", "T cell CD8", "T cell dividing", "T cell regulatory", "B cell",
#' "Monocyte conventional", "Monocyte non-conventional", "Macrophage", "NK cell"
#'
#' @format matrix with genes x cells
#'
#' @source Maynard et al. (2020), Cell,
#'   https://www.sciencedirect.com/science/article/pii/S0092867420308825
"single_cell_data"

#' Small example sigle-cell RNA-seq dataset from Maynard et al.
#'
#' A small subset of TPM-normalized single-cell RNA-seq gene expression, containing 50 cells with
#' types "T cell CD4", "T cell CD8", "B cell", "Monocyte conventional", "NK cell" and only 100 genes
#'
#' @format matrix with genes x cells
#'
#' @source Maynard et al. (2020), Cell,
#'   https://www.sciencedirect.com/science/article/pii/S0092867420308825
"single_cell_data_small"

#' Cell type annotations for the example sigle-cell RNA-seq dataset from Maynard et al.
#'
#' cell types are in the same order as cells in the sc RNA-seq matrix
#'
#' @format vector with 300 elements
#'
#' @source Maynard et al. (2020), Cell,
#'   https://www.sciencedirect.com/science/article/pii/S0092867420308825
"cell_type_annotations"

#' Cell type annotations for the small example sigle-cell RNA-seq dataset from Maynard et al.
#'
#' cell types are in the same order as cells in the small sc RNA-seq matrix
#'
#' @format vector with 50 elements
#'
#' @source Maynard et al. (2020), Cell,
#'   https://www.sciencedirect.com/science/article/pii/S0092867420308825
"cell_type_annotations_small"
