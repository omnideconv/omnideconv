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
#' @param alpha One positive numerical parameter or a numeircal vector of length equal nrow(input.phi),
#'  denoting the dirichlet hyper-parameter. Default=1, which represents a uniform prior over the
#'  simplex of theta. For sparser priors, use 0<alpha<1. Note that alpha usually does not affect
#'  the results, due to the dominating likelihood term resulted from the high sequencing depth of bulk RNA-seq.
#' @param sigma One positive numerical parameter or a numeircal vector of length equal number of genes
#' (for gene-specific prior), denoting the prior of the standard deviation of log fold change between
#' the true expression and the reference.Default=2, which represents a weak gene-wise prior. User may
#' provide their own sigma based on prior knowledge, such as differential expression analysis.
#' @param outlier_cut,outlier_fraction  Filter genes in X whose expression fraction is greater than
#' outlier.cut (Default=0.01. previous version used 0.05) in more than outlier.fraction (Default=0.1)
#'  of bulk data. Typically for dataset with reasonalble quality contol in mapping, very few genes
#'  will be filtered. Removal of outlier genes will ensure that the inference will not be dominated
#'  by outliers, which sometimes may be resulted from poor QC in mapping.
#' @param gibbs_control A list of parameters controling the Gibbs sampling. Default chain.length=1000,
#'  burn.in=500, thinning=2. A list of parameters controling the Gibbs sampling. Default chain.length=1000,
#'   burn.in=500, thinning=2. Previous version default is chain.length=400, burn.in=200, thinning=2.
#'   Default chain length has been increased to accomondate spatial transcriptomic data which usually
#'   has lower depth than conventional bulk data, and hence may need longer chain to reach the stationary distribution.
#' @param opt_control A list of parameters controling the optimization by Rcgmin, Default trace=0, maxit= 100000.
#' @param n_cores Number of CPU threads used for parallel computing. Default=1
#' @param n_cores_2g Number of CPU threads used for parallel computing for the final Gibbs sampling.
#'  Default=NULL (same as ncores). Recommended to set to a number smaller than n.cores, if deconvolving
#'  large number of mixtures, such as Visium data, or the number of cell types is large, to avoid memory overflow.
#' @param first_gibbs_only A logical parameter denoting if to only run the first gibbs sampling,
#' i.e. the initial estimates of theta and Z. Default: FALSE
#' @param seed A numerical number specifying the random seed number to generate identical results
#' between different runs. Default: NULL(ignore reproducibility).
#'
#' @return A list of results is returned including:
#'   \item{para}{All input data and parameters.}
#'   \item{res}{All output of TED. }
#'   \item{res$first.gibbs.res$gibbs.theta}{Initial estimates of fraction for all cell subtypes in each bulk sample.}
#'   \item{res$first.gibbs.res$Znkg}{Initial estimates of the mean of posterior read count for each cell subtypes  in each bulk sample.}
#'   \item{res$first.gibbs.res$theta.merged}{Initial estimates of fraction summed across cell types in each bulk sample.}
#'   \item{res$first.gibbs.res$Znkg.merged}{Initial estimates of the mean of posterior read count summed across cell types in each bulk sample.}
#'   \item{res$Zkg.tum}{Mean of posterior of gene expression of tumor in each patient.}
#'   \item{res$Zkg.tum.norm}{Depth normalized Zkg.tum (A pseudo count is added, such that the zero-valued genes have the same value as the min(phi.input)). Refered to as the psi.tum in the TED paper)}
#'   \item{res$Zkg.tum.vst}{Variance stablized transformed value of Zkg.tum. If vst transformation is not feasible, return NULL. }
#'   \item{res$phi.env}{Batch effect corrected expression profiles of stromal cells (refered to as the psi.str in the TED paper)}
#'   \item{res$final.gibbs.theta}{Updated theta after batch correction and tumor expression estimates. This is the final deconvolution result.}
#'   \item{res$cor.mat}{The correlation matrix of the estimated tumor expression profiles across bulk RNA-seq samples.}
#' @export
#'
deconvolute_bayesprism <- function(bulk_gene_expression, single_cell_object, cell_type_annotations,
                                   cell_subtype_labels=NULL, tum_key=NULL, pseudo_min=1E-8, alpha=1,
                                   sigma=2, outlier_cut=0.01, outlier_fraction=0.1,
                                   gibbs_control=list(chain.length=1000,burn.in=500,thinning=2),
                                   opt_control=list(trace=0, maxit= 100000),n_cores=1,
                                   n_cores_2g=NULL, first_gibbs_only=FALSE, seed=NULL) {

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

  # package is required to run bayesprism, also with only a single core
  require(snowfall)

  return(TED::run.Ted(ref.dat = single_cell_object,
                      X = bulk_gene_expression,
                      cell.type.labels = cell_type_annotations,
                      cell.subtype.labels = cell_subtype_labels,
                      tum.key = tum_key,
                      input.type = 'scRNA',
                      pseudo.min = pseudo_min,
                      alpha = alpha,
                      sigma = sigma,
                      outlier.cut = outlier_cut,
                      outlier.fraction = outlier_fraction,
                      gibbs.control = gibbs_control,
                      opt.control = opt_control,
                      n.cores = n_cores,
                      n.cores.2g = n_cores_2g,
                      pdf.name = NULL,
                      first.gibbs.only = first_gibbs_only,
                      seed = seed))
}
