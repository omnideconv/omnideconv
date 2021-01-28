#' Calculates the decomposition using the dwls algorithm
#'
#' Generates a reference profile based on single-cell data. Learns a transformation of bulk expression based  on  observed  single-cell  proportions  and  performs  NNLS  regression  on  these  transformed values to estimate cell proportions.
#'
#' @param bulk_gene_expression An Expression Set containing bulk data.
#' @param signature The Signature matrix.
#' @param dwls_submethod Three alternative methods in DWLS: OLS, SVR, and DampenedWLS.

#'
#' @return A list. Slot bulk.props contains a matrix of cell type proportion estimates with cell types as rows and individuals as columns.
#' @export
#'
#' @examples

deconvolute_dwls = function(bulk_gene_expression, signature, dwls_submethod = c("OLS","SVR","DampenedWLS")){

  if (length(dwls_submethod)>1){
    dwls_submethod <- "OLS"
  }

  message("\nRunning DWLS deconvolution module\n")

  # trim data
  Genes<-intersect(rownames(signature),rownames(bulk_gene_expression))
  bulk<-bulk_gene_expression[Genes,]
  sig<-signature[Genes,]
  if (class(bulk)[[1]]=="numeric"||class(sig)[[1]]=="numeric"){
    base::stop("Either bulk data or signature matrix just contains one row!")
  }

  # perform reconvolution in different sub_methods
  res <- NULL

  if(dwls_submethod == "OLS"){
    solutionsOLS<-NULL
    for (i in 1:ncol(bulk)){
      bulk_i<-bulk[,i]
      sol<-solveOLS(sig,bulk_i)
      sol<-round(sol,5)
      solutionsOLS<-cbind(solutionsOLS,sol)
    }
    colnames(solutionsOLS)<-colnames(bulk)
    res<-solutionsOLS
  }

  if(dwls_submethod == "SVR"){
    solutionsSVR<-NULL
    for (i in 1:ncol(bulk)){
      bulk_i<-bulk[,i]
      sol<-solveSVR(sig,bulk_i)
      sol<-round(sol,5)
      solutionsSVR<-cbind(solutionsSVR,sol)
    }
    colnames(solutionsSVR)<-colnames(bulk)
    res<-solutionsSVR
  }

  if(dwls_submethod == "DampenedWLS"){
    solutionsDampenedWLS<-NULL
    for (i in 1:ncol(bulk)){
      bulk_i<-bulk[,i]
      sol<-solveDampenedWLS(sig,bulk_i)
      sol<-round(sol,5)
      solutionsDampenedWLS<-cbind(solutionsDampenedWLS,sol)
    }
    colnames(solutionsDampenedWLS)<-colnames(bulk)
    res <- solutionsDampenedWLS
  }
  message("Deconvolution sucessful!")
  return (t(res))
}

