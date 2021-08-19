library(BisqueRNA)
library(MOMF)
library(Seurat)
library(MuSiC)
library(SCDC)
library(scBio)
library(CDSeq)
library(DWLS)

bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv", row.names = 1))
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- readr::read_lines("small_test_data/cell_annotations_small.txt")
batch_ids_small_from_file <- readr::read_lines("small_test_data/batch_ids_small.txt")
marker_genes <- readr::read_lines("small_test_data/marker_genes_small.txt")

markers_small <- list(marker_genes[1:9], marker_genes[10:14], marker_genes[15:20])
names(markers_small) <- sort(unique(cell_annotations_small))



matrixRightFormat <- sc_object_small
# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(
  check.names = F, check.rows = F,
  stringsAsFactors = F,
  row.names = colnames(sc_object_small),
  batchId = batch_ids_small_from_file,
  cellType = cell_annotations_small
)
sc.meta <- data.frame(
  labelDescription = c(
    "batchId",
    "cellType"
  ),
  row.names = c(
    "batchId",
    "cellType"
  )
)
sc.pdata <- new("AnnotatedDataFrame",
  data = sc.pheno,
  varMetadata = sc.meta
)
colnames(matrixRightFormat) <- row.names(sc.pdata)
rownames(matrixRightFormat) <- rownames(sc_object_small)
sc.eset <- Biobase::ExpressionSet(
  assayData = matrixRightFormat,
  phenoData = sc.pdata
)
sc.eset <- Biobase::ExpressionSet(
  assayData = Biobase::exprs(sc.eset),
  phenoData = sc.eset@phenoData
)
single_cell_expression_set <- sc.eset
bulk_expression_set <- Biobase::ExpressionSet(assayData = bulk_small)

sc.eset <- BisqueRNA:::CountsToCPM(sc.eset)
sc.eset <- BisqueRNA:::FilterZeroVarianceGenes(sc.eset, FALSE)


model_bisque <- BisqueRNA::GenerateSCReference(sc.eset, "cellType")
utils::write.csv(model_bisque, "test_models/bisque_model_small.csv")
result_bisque <- t(BisqueRNA::ReferenceBasedDecomposition(
  Biobase::ExpressionSet(assayData = bulk_small), single_cell_expression_set,
  use.overlap = FALSE, subject.names = "batchId"
)$bulk.props)
result_bisque <- result_bisque[, order(colnames(result_bisque))]
utils::write.csv(result_bisque, "test_results/bisque_result_small.csv")

# MOMF
model_momf <- MOMF::momf.computeRef(sc_object_small, cell_annotations_small)
utils::write.csv(model_momf, "test_models/momf_model_small.csv")
GList <- list(X1 = t(sc_object_small), X2 = t(bulk_small))
result_momf <- MOMF::momf.fit(DataX = GList, DataPriorU = model_momf)$cell.prop
result_momf <- result_momf[, order(colnames(result_momf))]
utils::write.csv(result_momf, "test_results/momf_result_small.csv")

# DWLS
cell_annotations_small_temp <- gsub(" ", "_", cell_annotations_small)
model_dwls <- DWLS::buildSignatureMatrixMAST(sc_object_small, cell_annotations_small_temp, tempdir())
tr <- trimDataOwn(model_dwls, bulk_small)

allCounts_DWLS <- NULL
allCounts_OLS <- NULL
allCounts_SVR <- NULL
for (i in 1:ncol(tr$bulk)) {
  bulk_i <- tr$bulk[, i]
  solOLS <- DWLS::solveOLS(tr$sig, bulk_i)
  solDWLS <- DWLS::solveDampenedWLS(tr$sig, bulk_i)
  solSVR <- DWLS::solveSVR(tr$sig, bulk_i)

  allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
  allCounts_OLS <- cbind(allCounts_OLS, solOLS)
  allCounts_SVR <- cbind(allCounts_SVR, solSVR)
}

colnames(allCounts_DWLS) <- colnames(tr$bulk)
colnames(allCounts_OLS) <- colnames(tr$bulk)
colnames(allCounts_SVR) <- colnames(tr$bulk)

colnames(model_dwls) <- gsub("_", " ", colnames(model_dwls))
utils::write.csv(model_dwls, "test_models/dwls_model_small.csv")
allCounts_DWLS <- t(allCounts_DWLS)
allCounts_OLS <- t(allCounts_OLS)
allCounts_SVR <- t(allCounts_SVR)
allCounts_DWLS <- allCounts_DWLS[, order(colnames(allCounts_DWLS))]
allCounts_OLS <- allCounts_OLS[, order(colnames(allCounts_OLS))]
allCounts_SVR <- allCounts_SVR[, order(colnames(allCounts_SVR))]
colnames(allCounts_DWLS) <- gsub("_", " ", colnames(allCounts_DWLS))
colnames(allCounts_OLS) <- gsub("_", " ", colnames(allCounts_OLS))
colnames(allCounts_SVR) <- gsub("_", " ", colnames(allCounts_SVR))
utils::write.csv(allCounts_DWLS, "test_results/dwls_dwls_result_small.csv")
utils::write.csv(allCounts_OLS, "test_results/dwls_ols_result_small.csv")
utils::write.csv(allCounts_SVR, "test_results/dwls_svr_result_small.csv")




## CibersortX
single_cell <- rbind(cell_annotations_small, sc_object_small)
rownames(single_cell) <- c("GeneSymbol", rownames(sc_object_small))
single_cell <- data.frame("GeneSymbol" = rownames(single_cell), single_cell)

write.table(single_cell, "sample_file_for_cibersort.txt",
  sep = "\t", quote = F, row.names = F,
  col.names = F
)
write.table(data.frame("Gene" = rownames(bulk_small), bulk_small), "mixture_file_for_cibersort.txt",
  sep = "\t", quote = F, row.names = F
)

root <- getwd() # Alternatively rstudioapi::getSourceEditorContext()$path can be used
email <- "konstantin.pelz@tum.de"
token <- "27308ae0ef1458d381becac46ca7e480"
signature_command <- paste0(
  "docker run -v ", root, ":/src/data -v ", root,
  "/test_models:/src/outdir cibersortx/fractions --username ", email,
  " --token ", token,
  " --single_cell TRUE --refsample sample_file_for_cibersort.txt"
)
system(signature_command)
file.rename(
  paste0(
    root, "/test_models/CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx",
    "_sample_file_for_cibersort_inferred_refsample.bm.K999.txt"
  ),
  paste0(root, "/test_models/cibersortx_model_small.tsv")
)
file.remove(paste0(
  root, "/test_models/CIBERSORTx_sample_file_for_cibersort_inferred",
  "_phenoclasses.CIBERSORTx_sample_file_for_cibersort_inferred",
  "_refsample.bm.K999.pdf"
))
file.remove(paste0(
  root, "/test_models/CIBERSORTx_sample_file_for_cibersort_inferred",
  "_phenoclasses.txt"
))
file.remove(paste0(
  root, "/test_models/CIBERSORTx_sample_file_for_cibersort_inferred",
  "_refsample.txt"
))
file.remove(paste0(root, "/test_models/CIBERSORTx_cell_type_sourceGEP.txt"))

result_command <- paste0(
  "docker run -v ", root, ":/src/data -v ", root,
  "/test_results:/src/outdir cibersortx/fractions --username ", email,
  " --token ", token,
  " --single_cell TRUE --mixture mixture_file_for_cibersort.txt",
  " --sigmatrix cibersortx_model_small.tsv"
)
file.copy(
  paste0(root, "/test_models/cibersortx_model_small.tsv"),
  paste0(root, "/cibersortx_model_small.tsv")
)
system(result_command)
file.rename(
  paste0(root, "/test_results/CIBERSORTx_Results.txt"),
  paste0(root, "/test_results/cibersortx_result_small.tsv")
)

file.remove(paste0(root, "/cibersortx_model_small.tsv"))
file.remove(paste0(root, "/sample_file_for_cibersort.txt"))
file.remove(paste0(root, "/mixture_file_for_cibersort.txt"))


## AutoGeneS
# First the data was saved as a h5ad file with the 'save_as_h5ad' method.
# Then, the following python script was used:
#
# import pandas as pd
# import autogenes as ag
# import scanpy as sc
#
# adata = sc.read_h5ad("test.h5ad")
# ag.init(adata,celltype_key='label')
# ag.optimize()
# bulk_data = pd.read_csv('bulk_small.csv',index_col=0).transpose()
# coef = ag.deconvolve(bulk_data)
# res = pd.DataFrame(coef,columns=ag.adata().obs_names,index=bulk_data.index)
# res.to_csv("autogenes_result_small.csv",sep="\t")

## MuSiC
result_music <- music_prop(
  bulk.eset = bulk_expression_set, sc.eset = single_cell_expression_set,
  clusters = "cellType", samples = "batchId"
)$Est.prop.weighted
result_music <- result_music[, order(colnames(result_music))]
utils::write.csv(result_music, "test_results/music_result_small.csv")

## SCDC
result_scdc <- SCDC_prop(
  bulk.eset = bulk_expression_set, sc.eset = single_cell_expression_set,
  ct.varname = "cellType", sample = "batchId",
  ct.sub = unique(single_cell_expression_set@phenoData@data[, "cellType"])
)$prop.est.mvw
result_scdc <- result_scdc[, order(colnames(result_scdc))]
utils::write.csv(result_scdc, "test_results/scdc_result_small.csv")

eset_one <- getESET(single_cell_data,
  fdata = rownames(single_cell_data),
  pdata = cbind(
    cellname = cell_type_annotations,
    subjects = batch_ids
  )
)
eset_two <- getESET(sc_object_small,
  fdata = rownames(sc_object_small),
  pdata = cbind(
    cellname = cell_annotations_small,
    subjects = batch_ids_small_from_file
  )
)

eset_list <- list(one = eset_one, two = eset_two)

result_scdc_ensemble <- SCDC_ENSEMBLE(
  bulk.eset = bulk_expression_set, sc.eset.list = eset_list,
  ct.varname = "cellname", sample = "subjects",
  ct.sub = Reduce(intersect, sapply(eset_list, function(x) {
    unique(x@phenoData@data[, "cellname"])
  }))
)
props <- wt_prop(result_scdc_ensemble$w_table, result_scdc_ensemble$prop.only)
props <- props[, order(colnames(props))]
utils::write.csv(props, "test_results/scdc_result_ensemble.csv")


## BSEQ-sc

signature_bseqsc <-
  bseqsc_basis(sc_object_small, markers_small, cell_annotations_small, batch_ids_small_from_file)
utils::write.csv(signature_bseqsc, "test_models/bseq_model_small.csv")
# bseqsc_config("C:/Users/Konstantin/Downloads/CIBERSORT.R")
# bseqsc_props <- bseqsc_proportions(bulk_small, reference = signature_bseqsc)
# bseqsc_props <- t(bseqsc_props$coefficients)
# bseqsc_props <- bseqsc_props[, order(colnames(bseqsc_props))]
# utils::write.csv(bseqsc_props, "test_results/bseqsc_result_small.csv")


## CDSeq

# cdseq_res <- CDSeq::CDSeq(
#  bulk_data = bulk,
#  cell_type_number = length(unique(cell_type_annotations))
# )

# cdseq_gep <- cdseq_res$estGEP
# cdseq_prop <- cdseq_res$estProp
# cdseq_prop
# cell_annotations_as_df <- cbind(colnames(single_cell_data), cell_type_annotations)
# colnames(cell_annotations_as_df) <- c("cell_id", "cell_type")
# cell_annotations_as_df <- data.frame(cell_annotations_as_df)
# props <- CDSeq::cellTypeAssignSCRNA(
#  cdseq_gep = cdseq_gep, # CDSeq-estimated cell-type-specific GEPs
#  cdseq_prop = cdseq_prop, # CDSeq-estimated cell type proportions
#  sc_gep = single_cell_data, # PBMC single cell data
#  sc_annotation = cell_annotations_as_df, # PBMC single data annotations
#  sc_batch = batch_ids,
#  plot_umap = 0,
#  plot_tsne = 0,
#  verbose = TRUE
# )
# props$input_list$cdseq_prop
# props$cdseq_prop_merged




## DWLS Stuff

# trim bulk and single-cell data to contain the same genes
# CHANGED FROM THE ORIGINAL FUNCTION TO WORK WITH OUR DATA
trimDataOwn <- function(Signature, bulkData) {
  Genes <- intersect(rownames(Signature), rownames(bulkData))
  B <- bulkData[Genes, ]
  S <- Signature[Genes, ]
  return(list("sig" = S, "bulk" = B))
}
