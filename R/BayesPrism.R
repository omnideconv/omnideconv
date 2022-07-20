#' No model is built as BayesPrism performs signature building and deconvolution in one step.
#'
#' Please use the deconvolute method with your single cell and bulk rna seq data to use BayesPrism
#'
#'
#' @return NULL.
#'
#' @export
build_model_bayesprism <- function() {
  message(
    "The deconvolution with BayesPrism is done in only one step. Please just use the ",
    "deconvolute method."
  )

  return(NULL)
}


#' Bayesian deconvolution module using BayesPrism
#'
#' IMPORTANT: no model is needed. Everything is done inside this method.
#'
#' Run Bayesian deconvolution to estimate cell type composition and gene expression.
#'
#' @param bulk_gene_expression A matrix of bulk data. Rows are genes, columns are samples.
#'   Row and column names need to be set.
#' @param single_cell_object A matrix with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param cell_type_annotations A vector of the cell type annotations. Has to be in the same order
#'   as the samples in single_cell_object.
#' @param cell_subtype_labels a character or factor vector indicating the cell subtype of each row of
#'  the raw count matrix of scRNA-seq or gene expression profile (GEP). The length needs be equal to
#'  nrow(ref.dat). Default is NULL, which uses the same value of cell.type.labels. Note that TED
#'  computes the posterior sum over the subtypes to get the total fraction / expression of each cell type.
#'   This allows a more fine-grained definition of cell types / cell states.
#' @param tum_key The character in cell.type.labels denoting the tumor cells, e.g. "tumor" or "malignant".
#' @param pseudo_min A numeric value indicating the minimum (non-zero) value of phi. Default=1E-8.
#' @param gibbs_control A list of parameters controlling the Gibbs sampling. Default chain.length=1000,
#'  burn.in=500, thinning=2. A list of parameters controlling the Gibbs sampling. Default chain.length=1000,
#'   burn.in=500, thinning=2. Previous version default is chain.length=400, burn.in=200, thinning=2.
#'   Default chain length has been increased to accommodate spatial transcriptomic data which usually
#'   has lower depth than conventional bulk data, and hence may need longer chain to reach the stationary distribution.
#' @param update_gibbs A logical variable to denote whether run final Gibbs sampling to update theta. Default=TRUE.
#' @param opt_control A list of parameters controlling the optimization by Rcgmin, Default trace=0, maxit= 100000.
#' @param apply_bayes_prism_filtering set to TRUE if you want to run the `cleanup.genes` function of BayesPrism; default is FALSE
#' @param species A character variable to denote if genes are human ("mm") or mouse ("hs").
#' @param exp.cells Genes expressed in number of cells fewer than this will be excluded. Default=1. If the input is GEP, gene
#' will be selected by automatically setting exp.cells is set to min(exp.cells,1). As a result genes expressed in at least 0
#' or 1 cell type will be retained. Only used when `apply_bayes_prism_filtering` is TRUE.
#' @param gene_group a character vector to input gene groups to be removed, must be one or more elements from
#' c("other_Rb","chrM","chrX","chrY","Rb","Mrp","act","hb","MALAT1"). Only used when `apply_bayes_prism_filtering` is TRUE.
#' @param outlier_cut,outlier_fraction Filter genes in X whose expression fraction is greater than outlier.cut
#' (Default=0.01) in more than outlier.fraction (Default=0.1) of bulk data. Typically for dataset with reasonable
#'  quality control, very few genes will be filtered. Removal of outlier genes will ensure that the inference
#'   will not be dominated by outliers, which sometimes may be resulted from poor QC in mapping.
#' @param which_theta A character variable to denote whether to extract results from first or final Gibbs sampling.
#' @param state_or_type A character variable to extract results from cell type or cell state. We caution the
#' interpretation of cell states information when their transcription are highly co-linear.
#' @param n_cores Number of CPU threads used for parallel computing. Default=1
#'
#' @return A list of results is returned including:
#' \itemize{
#'   \item{bp.res}{The result of `run.prism`, a "BayesPrism" S4 object.}
#'   \item{theta}{The result of `get.fraction`, the extracted cell fraction results from the BayesPrism object.}
#'   \item{bp.res$prism}{An S4 object of the class "prism" to represent the input prism object; }
#'   \item{bp.res$posterior.initial.cellState}{An S4 object of the class "jointPost" to represent the posterior mean of cell state fraction and cell state expression outputted by the initial Gibbs sampling using cell state pirors. Contains Z (inferred expression), theta (inferred fraction) and theta.cv (coefficient of variation of posterior of theta)}
#'   \item{bp.res$posterior.initial.cellType}{An S4 object of the class "jointPost" to represent the posterior sum of cell states from each cell type (posterior.initial.cellState); }
#'   \item{bp.res$reference.update}{An S4 obejct of the class "reference" to represent the updated profile ψ;}
#'   \item{bp.res$posterior.theta_f:}{An S4 object of the class "thetaPost" to represent the updated cell type fraction. Contains theta (inferred fraction) and theta.cv (coefficient of variation of posterior of theta); }
#'   \item{bpres$control_param}{A list storing the gibbs.control, opt.control and update.gibbs arguments.}
#' }
#' @export
#'
#' @import snowfall
deconvolute_bayesprism <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                                   cell_subtype_labels = NULL, tum_key = NULL, apply_bayes_prism_filtering = FALSE,
                                   species = "hs", exp.cells = 1, pseudo_min = 1E-8,
                                   gene_group = c("other_Rb", "chrM", "chrX", "chrY", "Rb", "Mrp", "act", "hb", "MALAT1"),
                                   outlier_cut = 0.01, outlier_fraction = 0.1, update_gibbs = TRUE,
                                   gibbs_control = list(chain.length = 1000, burn.in = 500, thinning = 2),
                                   opt_control = list(trace = 0, maxit = 100000), n_cores = 1,
                                   which_theta = "final", state_or_type = "type") {
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }

  ## BayesPrism expects the bulk and single-cell matrices in a transposed format; genes are in columns and cells in rows
  bulk_gene_expression <- t(bulk_gene_expression)
  single_cell_object <- t(single_cell_object)

  # possibility to run the custom gene filtering function
  if (apply_bayes_prism_filtering) {
    message("Cleaning up genes using BayesPrism ...")
    single_cell_object <- BayesPrism::cleanup.genes(
      input = single_cell_object,
      input.type = "count.matrix",
      species = species,
      gene.group = gene_group,
      exp.cells = exp.cells
    )
  }

  # construct BayesPrism object
  myPrism <- BayesPrism::new.prism(
    reference = single_cell_object,
    mixture = bulk_gene_expression,
    input.type = "count.matrix",
    cell.type.labels = cell_type_annotations,
    cell.state.labels = cell_subtype_labels,
    key = tum_key,
    outlier.cut = outlier_cut,
    outlier.fraction = outlier_fraction,
    pseudo.min = pseudo_min
  )

  # run deconvolution
  bp.res <- BayesPrism::run.prism(
    prism = myPrism,
    n.cores = n_cores,
    update.gibbs = update_gibbs,
    gibbs.control = gibbs_control,
    opt.control = opt_control
  )

  # extract cell type fractions from result object
  theta <- BayesPrism::get.fraction(
    bp = bp.res,
    which.theta = which_theta,
    state.or.type = state_or_type
  )

  return(list(
    theta = theta,
    bp.res = bp.res
  ))
}
