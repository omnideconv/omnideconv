library(BisqueRNA)
library(MOMF)
library(reshape)
library(quadprog)
library(e1071)
library(Seurat)
library(ROCR)
library(varhandle)
library(MAST)
library(scBio)

bulk_small <- as.matrix(utils::read.csv("small_test_data/bulk_small.csv", row.names = 1))
sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- utils::read.csv("small_test_data/cell_annotations_small.csv",
  row.names = 1
)$x


matrixRightFormat <- sc_object_small
# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(
  check.names = F, check.rows = F,
  stringsAsFactors = F,
  row.names = colnames(sc_object_small),
  SubjectName = colnames(sc_object_small),
  cellType = cell_annotations_small
)
sc.meta <- data.frame(
  labelDescription = c(
    "SubjectName",
    "cellType"
  ),
  row.names = c(
    "SubjectName",
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
  Biobase::ExpressionSet(assayData = bulk_small), sc.eset,
  use.overlap = FALSE
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
model_dwls <- buildSignatureMatrixMAST(sc_object_small, cell_annotations_small_temp, tempdir())
tr <- trimData(model_dwls, bulk_small)

allCounts_DWLS <- NULL
allCounts_OLS <- NULL
allCounts_SVR <- NULL
for (i in 1:ncol(tr$bulk)) {
  bulk_i <- tr$bulk[, i]
  solOLS <- solveOLS(tr$sig, bulk_i)
  solDWLS <- solveDampenedWLS(tr$sig, bulk_i)
  solSVR <- solveSVR(tr$sig, bulk_i)

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

root <- dirname(getwd()) # Alternatively rstudioapi::getSourceEditorContext()$path can be used
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
  clusters = "cellType", samples = "SubjectName"
)$Est.prop.weighted
result_music <- result_music[, order(colnames(result_music))]
utils::write.csv(result_music, "test_results/music_result_small.csv")



## CPM
# cellSpace <- scores(pcaMethods::pca(t(sc_object_small), method="svd", nPcs=2))
# Cell space creation with seurat instead
# result_cpm <- CPM(sc_object_small, cell_annotations_small, bulk_small, cellSpace,
#                  quantifyTypes = TRUE, typeTransformation = TRUE)$cellTypePredictions
# result_cpm <- result_cpm[, order(colnames(result_cpm))]
# utils::write.csv(result_cpm, "test_results/cpm_result_small.csv")


## DWLS Stuff

# trim bulk and single-cell data to contain the same genes
# CHANGED FROM THE ORIGINAL FUNCTION TO WORK WITH OUR DATA
trimData <- function(Signature, bulkData) {
  Genes <- intersect(rownames(Signature), rownames(bulkData))
  B <- bulkData[Genes, ]
  S <- Signature[Genes, ]
  return(list("sig" = S, "bulk" = B))
}


# solve using OLS, constrained such that cell type numbers>0
solveOLS <- function(S, B) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")
  solution <- solve.QP(D / sc, d / sc, A, bzero)$solution
  names(solution) <- colnames(S)
  # print(round(solution/sum(solution),5))
  return(solution / sum(solution))
}

# return cell number, not proportion
# do not print output
# CHANGED FROM THE ORIGINAL FUNCTION TO WORK
solveOLSInternal <- function(S, B) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))


  sc <- norm(D, "2")
  solution <- solve.QP(D / sc, d / sc, A, bzero)$solution
  # solution<-solve.QP(D, d,A,bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

# solve using WLS with weights dampened by a certain dampening constant
solveDampenedWLS <- function(S, B) {
  # first solve OLS, use this solution to find a starting point for the weights
  solution <- solveOLSInternal(S, B)
  # now use dampened WLS, iterate weights until convergence
  iterations <- 0
  changes <- c()
  # find dampening constant for weights using cross-validation
  j <- findDampeningConstant(S, B, solution)
  change <- 1
  while (change > .01 & iterations < 1000) {
    newsolution <- solveDampenedWLSj(S, B, solution, j)
    # decrease step size for convergence
    solutionAverage <- rowMeans(cbind(newsolution, matrix(solution,
      nrow = length(solution),
      ncol = 4
    )))
    change <- norm(as.matrix(solutionAverage - solution))
    solution <- solutionAverage
    iterations <- iterations + 1
    changes <- c(changes, change)
  }
  # print(round(solution/sum(solution),5))
  return(solution / sum(solution))
}

# solve WLS given a dampening constant
solveDampenedWLSj <- function(S, B, goldStandard, j) {
  multiplier <- 1 * 2^(j - 1)
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsDampened <- wsScaled
  wsDampened[which(wsScaled > multiplier)] <- multiplier
  W <- diag(wsDampened)
  D <- t(S) %*% W %*% S
  d <- t(S) %*% W %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")
  solution <- solve.QP(D / sc, d / sc, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

# find a dampening constant for the weights using cross-validation
findDampeningConstant <- function(S, B, goldStandard) {
  solutionsSd <- NULL
  # goldStandard is used to define the weights
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsScaledMinusInf <- wsScaled
  # ignore infinite weights
  if (max(wsScaled) == "Inf") {
    wsScaledMinusInf <- wsScaled[-which(wsScaled == "Inf")]
  }
  # try multiple values of the dampening constant (multiplier)
  # for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))) {
    multiplier <- 1 * 2^(j - 1)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    solutions <- NULL
    seeds <- c(1:100)
    for (i in 1:100) {
      set.seed(seeds[i]) # make nondeterministic
      subset <- sample(length(ws), size = length(ws) * 0.5) # randomly select half of gene set
      # solve dampened weighted least squares for subset
      fit <- lm(B[subset] ~ -1 + S[subset, , drop = FALSE], weights = wsDampened[subset])
      sol <- fit$coef * sum(goldStandard) / sum(fit$coef)
      solutions <- cbind(solutions, sol)
    }
    solutionsSd <- cbind(solutionsSd, apply(solutions, 1, sd))
  }
  # choose dampening constant that results in least cross-validation variance
  j <- which.min(colMeans(solutionsSd^2))
  return(j)
}

solveSVR <- function(S, B) {
  # scaling
  ub <- max(c(as.vector(S), B)) # upper bound
  lb <- min(c(as.vector(S), B)) # lower bound
  Bs <- ((B - lb) / ub) * 2 - 1
  Ss <- ((S - lb) / ub) * 2 - 1

  # perform SVR
  model <- svm(Ss, Bs, nu = 0.5, scale = TRUE, type = "nu-regression", kernel = "linear", cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef < 0)] <- 0
  coef <- as.vector(coef)
  names(coef) <- colnames(S)
  print(round(coef / sum(coef), 5))
  return(coef / sum(coef))
}

## alternative differential expression method using MAST

# functions for DE

Mean.in.log2space <- function(x, pseudo.count) {
  return(log2(mean(2^(x) - pseudo.count) + pseudo.count))
}

stat.log2 <- function(data.m, group.v, pseudo.count) {
  # data.m=data.used.log2
  log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) {
    Mean.in.log2space(x, pseudo.count)
  })
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group", log2.mean.r[1, ], sep = "")
  log2.mean.r <- log2.mean.r[-1, ]
  log2.mean.r <- as.data.frame(log2.mean.r)
  log2.mean.r <- varhandle::unfactor(log2.mean.r) # from varhandle
  log2.mean.r[, 1] <- as.numeric(log2.mean.r[, 1])
  log2.mean.r[, 2] <- as.numeric(log2.mean.r[, 2])
  log2_foldchange <- log2.mean.r$mean.group1 - log2.mean.r$mean.group0
  results <- data.frame(cbind(log2.mean.r$mean.group0, log2.mean.r$mean.group1, log2_foldchange))
  colnames(results) <- c("log2.mean.group0", "log2.mean.group1", "log2_fc")
  rownames(results) <- rownames(log2.mean.r)
  return(results)
}

v.auc <- function(data.v, group.v) {
  prediction.use <- prediction(data.v, group.v, 0:1)
  perf.use <- performance(prediction.use, "auc")
  auc.use <- round(perf.use@y.values[[1]], 3)
  return(auc.use)
}
m.auc <- function(data.m, group.v) {
  AUC <- apply(data.m, 1, function(x) v.auc(x, group.v))
  AUC[is.na(AUC)] <- 0.5
  return(AUC)
}

# perform DE analysis using MAST
DEAnalysisMAST <- function(scdata, id, path) {
  pseudo.count <- 0.1
  data.used.log2 <- log2(scdata + pseudo.count)
  colnames(data.used.log2) <- make.unique(colnames(data.used.log2))
  diff.cutoff <- 0.5
  for (i in unique(id)) {
    cells.symbol.list2 <- colnames(data.used.log2)[which(id == i)]
    cells.coord.list2 <- match(cells.symbol.list2, colnames(data.used.log2))
    cells.symbol.list1 <- colnames(data.used.log2)[which(id != i)]
    cells.coord.list1 <- match(cells.symbol.list1, colnames(data.used.log2))
    data.used.log2.ordered <- cbind(
      data.used.log2[, cells.coord.list1],
      data.used.log2[, cells.coord.list2]
    )
    group.v <- c(rep(0, length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    # ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))

    DE <- bigtable[bigtable$log2_fc > diff.cutoff, ]
    dim(DE)
    if (dim(DE)[1] > 1) {
      data.1 <- data.used.log2[, cells.coord.list1, drop = FALSE]
      data.2 <- data.used.log2[, cells.coord.list2, drop = FALSE]
      genes.list <- rownames(DE)
      log2fold_change <- cbind(genes.list, DE$log2_fc)
      colnames(log2fold_change) <- c("gene.name", "log2fold_change")
      counts <- as.data.frame(cbind(data.1[genes.list, ], data.2[genes.list, ]))
      groups <- c(
        rep("Cluster_Other", length(cells.coord.list1)),
        rep(i, length(cells.coord.list2))
      )
      groups <- as.character(groups)
      data_for_MIST <- as.data.frame(cbind(
        rep(rownames(counts), dim(counts)[2]), melt(counts),
        rep(groups, each = dim(counts)[1]),
        rep(1, dim(counts)[1] * dim(counts)[2])
      ))
      colnames(data_for_MIST) <- c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta <- data_for_MIST
      vbeta.fa <- FromFlatDF(vbeta,
        idvars = c("Subject.ID"),
        primerid = "Gene", measurement = "Et", ncells = "Number.of.Cells",
        geneid = "Gene", cellvars = c("Number.of.Cells", "Population"),
        phenovars = c("Population"), id = "vbeta all"
      )
      vbeta.1 <- subset(vbeta.fa, Number.of.Cells == 1)
      # .3 MAST
      # head(colData(vbeta.1)) ??
      zlm.output <- zlm(~Population, vbeta.1, method = "bayesglm", ebayes = TRUE)
      show(zlm.output)
      coefAndCI <- summary(zlm.output, logFC = TRUE)
      zlm.lr <- lrTest(zlm.output, "Population")
      zlm.lr_pvalue <- melt(zlm.lr[, , "Pr(>Chisq)"])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == "hurdle"), ]



      lrTest.table <- merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <- c(
        "Gene", "test.type", "p_value",
        paste("log2.mean.", "Cluster_Other", sep = ""),
        paste("log2.mean.", i, sep = ""), "log2fold_change", "Auc"
      )
      cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)), ]

      # . 4 save results
      utils::write.csv(cluster_lrTest.table, file = paste(path, "/", i, "_lrTest.csv", sep = ""))
      save(cluster_lrTest.table, file = paste(path, "/", i, "_MIST.RData", sep = ""))
    }
  }
}

# build signature matrix using genes identified by DEAnalysisMAST()
buildSignatureMatrixMAST <- function(scdata, id, path, diff.cutoff = 0.5, pval.cutoff = 0.01) {
  # compute differentially expressed genes for each cell type
  DEAnalysisMAST(scdata, id, path)

  # for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the
  # estimated fold-change is greater than 0.5
  numberofGenes <- c()
  for (i in unique(id)) {
    if (file.exists(paste(path, "/", i, "_MIST.RData", sep = ""))) {
      load(file = paste(path, "/", i, "_MIST.RData", sep = ""))
      pvalue_adjusted <- p.adjust(cluster_lrTest.table[, 3],
        method = "fdr",
        n = length(cluster_lrTest.table[, 3])
      )
      cluster_lrTest.table <- cbind(cluster_lrTest.table, pvalue_adjusted)
      DEGenes <- cluster_lrTest.table$Gene[intersect(
        which(pvalue_adjusted < pval.cutoff),
        which(cluster_lrTest.table$log2fold_change > diff.cutoff)
      )]
      nonMir <- grep("MIR|Mir", DEGenes, invert = T) # because Mir gene is usually not accurate
      assign(
        paste("cluster_lrTest.table.", i, sep = ""),
        cluster_lrTest.table[which(cluster_lrTest.table$Gene %in% DEGenes[nonMir]), ]
      )
      numberofGenes <- c(numberofGenes, length(DEGenes[nonMir]))
    }
  }

  # need to reduce number of genes
  # for each subset, order significant genes by decreasing fold change, choose between 50 and
  # 200 genes for each, iterate and choose matrix with lowest condition number
  conditionNumbers <- c()
  for (G in 50:200) {
    Genes <- c()
    j <- 1
    for (i in unique(id)) {
      if (numberofGenes[j] > 0) {
        temp <- paste("cluster_lrTest.table.", i, sep = "")
        temp <- as.name(temp)
        temp <- eval(parse(text = temp))
        temp <- temp[order(temp$log2fold_change, decreasing = TRUE), ]
        Genes <- c(Genes, varhandle::unfactor(temp$Gene[1:min(G, numberofGenes[j])]))
      }
      j <- j + 1
    }
    Genes <- unique(Genes)
    # make signature matrix
    ExprSubset <- scdata[Genes, ]
    Sig <- NULL
    for (i in unique(id)) {
      Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
    }
    colnames(Sig) <- unique(id)
    conditionNumbers <- c(conditionNumbers, kappa(Sig))
  }
  G <- which.min(conditionNumbers) + min(49, numberofGenes - 1)
  Genes <- c()
  j <- 1
  for (i in unique(id)) {
    if (numberofGenes[j] > 0) {
      temp <- paste("cluster_lrTest.table.", i, sep = "")
      temp <- as.name(temp)
      temp <- eval(parse(text = temp))
      temp <- temp[order(temp$log2fold_change, decreasing = TRUE), ]
      Genes <- c(Genes, varhandle::unfactor(temp$Gene[1:min(G, numberofGenes[j])]))
    }
    j <- j + 1
  }
  Genes <- unique(Genes)
  ExprSubset <- scdata[Genes, ]
  Sig <- NULL
  for (i in unique(id)) {
    Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
  }
  colnames(Sig) <- unique(id)
  save(Sig, file = paste(path, "/Sig.RData", sep = ""))
  return(Sig)
}
