#' Prepare Design matrix and Cross-subject Variance for MuSiC Deconvolution
#'
#' This function is used for generating cell type specific cross-subject mean and variance for each
#' gene. Cell type specific library size is also calculated.
#'
#' @param bulk_gene_expression A matrix or dataframe with the bulk data. Rows
#'   are genes, columns are samples.
#' @param single_cell_object A Matrix with the single-cell data. Rows are genes
#'  and columns are samples.
#' @param cell_type_annotations A Vector of the cell type annotations. Has to
#'  be in the same order as the samples in single_cell_object.
#' @param non_zero logical, default as TRUE. If true, remove all gene with zero expression.
#' @param markers vector or list of gene names. Default as NULL. If NULL, then use all genes
#'  provided.
#' @param clusters character, the phenoData used as clusters;
#' @param samples character,the phenoData used as samples;
#' @param select_ct vector of cell types. Default as NULL. If NULL, then use all cell types
#'  provided.
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd
#'  column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size
#'  from data.
#' @param ct_cov logical. If TRUE, use the covariance across cell types.
#' @param verbose logical, default as TRUE.
#' @return a list of
#'     * gene by cell type matrix of Design matrix
#'     * subject by celltype matrix of Library size
#'     * vector of average library size for each cell type
#'     * gene by celltype matrix of average relative abudance
#'     * gene by celltype matrix of cross-subject variation
#'
#' @export
build_model_music <- function(bulk_gene_expression, single_cell_object = NULL,
                              cell_type_annotations = NULL, non_zero = TRUE, markers = NULL,
                              clusters = "cellType", samples = "SubjectName", select_ct = NULL,
                              cell_size = NULL, ct_cov = FALSE, verbose = FALSE) {
  if (is.null(single_cell_object) || is.null(cell_type_annotations)) {
    base::stop(
      "Single cell object or cell type annotations not provided. Call as: ",
      "build_model(bulk_gene_expression, signature, \"music\", single_cell_object, ",
      "cell_type_annotations)"
    )
  }
  sc_eset <- omnideconv:::get_single_cell_expression_set(
    single_cell_object, colnames(single_cell_object),
    rownames(single_cell_object), cell_type_annotations
  )
  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)

  bulk_gene <- rownames(bulk_eset)[rowMeans(exprs(bulk_eset)) != 0]
  bulk_eset <- bulk_eset[bulk_gene, , drop = FALSE]
  if (is.null(markers)) {
    sc_markers <- bulk_gene
  } else {
    sc_markers <- intersect(bulk_gene, unlist(markers))
  }
  sc_basis <- music_basis(sc_eset,
    non.zero = non_zero, markers = sc_markers, clusters = clusters,
    samples = samples, select.ct = select_ct, cell_size = cell_size,
    ct.cov = ct_cov, verbose = verbose
  )
  return(sc_basis)
}

#' MuSiC Deconvolution
#'
#' This function is to calculate the MuSiC deconvolution proportions
#'
#' @param bulk_gene_expression  A matrix or dataframe with the bulk data. Rows
#'   are genes, columns are samples.
#' @param signature_data Preprocessed single cell data from the build_model_music method
#' @param markers vector or list of gene names, default as NULL. If NULL, use all genes that
#'  provided by both bulk and single cell dataset.
#' @param clusters character, the phenoData of single cell dataset used as clusters;
#' @param samples character,the phenoData of single cell dataset used as samples;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd
#'  column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from
#'  data;
#' @param ct_cov logical. If TRUE, use the covariance across cell types;
#' @param verbose logical, default as TRUE.
#' @param iter.max numeric, maximum iteration number
#' @param nu regulation parameter, take care of weight when taking recipical
#' @param eps Thredshold of convergence
#' @param centered logic, substract avg of Y and D
#' @param normalize logic, divide Y and D by their standard deviation
#' @return a list with elements:
#'    * Estimates of MuSiC
#'    * Estimates of NNLS
#'    * Weight of MuSiC
#'    * r.squared of MuSiC
#'    * Variance of MuSiC estimates
#' @export
deconvolute_music <- function(bulk_gene_expression, signature_data, markers = NULL,
                              clusters = "cellType", samples = "SubjectName", cell_size = NULL,
                              ct_cov = FALSE, verbose = FALSE, iter.max = 1000, nu = 0.0001,
                              eps = 0.01, centered = FALSE, normalize = FALSE) {
  bulk_eset <- Biobase::ExpressionSet(assayData = bulk_gene_expression)

  bulk_gene <- rownames(bulk_eset)[rowMeans(exprs(bulk_eset)) != 0]
  bulk_eset <- bulk_eset[bulk_gene, , drop = FALSE]
  if (is.null(markers)) {
    sc_markers <- bulk_gene
  } else {
    sc_markers <- intersect(bulk_gene, unlist(markers))
  }
  cm_gene <- intersect(rownames(signature_data$Disgn.mtx), bulk_gene)
  if (is.null(markers)) {
    if (length(cm_gene) < 0.2 * min(length(bulk_gene), nrow(signature_data$Disgn.mtx))) {
      stop("Too few common genes!")
    }
  } else {
    if (length(cm_gene) < 0.2 * length(unlist(markers))) {
      stop("Too few common genes!")
    }
  }
  if (verbose) {
    message(paste("Used", length(cm_gene), "common genes..."))
  }

  m.sc <- match(cm_gene, rownames(signature_data$Disgn.mtx))
  m.bulk <- match(cm_gene, bulk_gene)
  D1 <- signature_data$Disgn.mtx[m.sc, ]
  M.S <- colMeans(signature_data$S, na.rm = T)

  if (!is.null(cell_size)) {
    if (!is.data.frame(cell_size)) {
      stop(
        "cell_size paramter should be a data.frame with 1st column for cell type names and 2nd",
        "column for cell sizes"
      )
    } else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
      stop("Cell type names in cell_size must match clusters")
    } else if (any(is.na(as.numeric(cell_size[, 2])))) {
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]), ]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }

  Yjg <- relative.ab(exprs(bulk_eset)[m.bulk, ])
  N.bulk <- ncol(bulk_eset)
  if (ct_cov) {
    Sigma.ct <- signature_data$Sigma.ct[, m.sc]

    Est.prop.allgene <- NULL
    Est.prop.weighted <- NULL
    Weight.gene <- NULL
    r.squared.full <- NULL
    Var.prop <- NULL

    for (i in 1:N.bulk) {
      if (sum(Yjg[, i] == 0) > 0) {
        D1.temp <- D1[Yjg[, i] != 0, ]
        Yjg.temp <- Yjg[Yjg[, i] != 0, i]
        Sigma.ct.temp <- Sigma.ct[, Yjg[, i] != 0]
        if (verbose) message(paste(colnames(Yjg)[i], "has common genes", sum(Yjg[, i] != 0), "..."))
      } else {
        D1.temp <- D1
        Yjg.temp <- Yjg[, i]
        Sigma.ct.temp <- Sigma.ct
        if (verbose) message(paste(colnames(Yjg)[i], "has common genes", sum(Yjg[, i] != 0), "..."))
      }

      lm.D1.weighted <- music.iter.ct(Yjg.temp, D1.temp, M.S, Sigma.ct.temp,
        iter.max = iter.max,
        nu = nu, eps = eps, centered = centered, normalize = normalize
      )
      Est.prop.allgene <- rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted <- rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp <- rep(NA, nrow(Yjg))
      weight.gene.temp[Yjg[, i] != 0] <- lm.D1.weighted$weight.gene
      Weight.gene <- cbind(Weight.gene, weight.gene.temp)
      r.squared.full <- c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop <- rbind(Var.prop, lm.D1.weighted$var.p)
    }
  } else {
    Sigma <- signature_data$Sigma[m.sc, ]

    valid.ct <- (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) == 0) & (!is.na(M.S))

    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }

    if (verbose) {
      message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    }

    D1 <- D1[, valid.ct]
    M.S <- M.S[valid.ct]
    Sigma <- Sigma[, valid.ct]

    Est.prop.allgene <- NULL
    Est.prop.weighted <- NULL
    Weight.gene <- NULL
    r.squared.full <- NULL
    Var.prop <- NULL
    for (i in 1:N.bulk) {
      if (sum(Yjg[, i] == 0) > 0) {
        D1.temp <- D1[Yjg[, i] != 0, ]
        Yjg.temp <- Yjg[Yjg[, i] != 0, i]
        Sigma.temp <- Sigma[Yjg[, i] != 0, ]
        if (verbose) message(paste(colnames(Yjg)[i], "has common genes", sum(Yjg[, i] != 0), "..."))
      } else {
        D1.temp <- D1
        Yjg.temp <- Yjg[, i]
        Sigma.temp <- Sigma
        if (verbose) message(paste(colnames(Yjg)[i], "has common genes", sum(Yjg[, i] != 0), "..."))
      }

      lm.D1.weighted <- music.iter(Yjg.temp, D1.temp, M.S, Sigma.temp,
        iter.max = iter.max,
        nu = nu, eps = eps, centered = centered, normalize = normalize
      )
      Est.prop.allgene <- rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted <- rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp <- rep(NA, nrow(Yjg))
      weight.gene.temp[Yjg[, i] != 0] <- lm.D1.weighted$weight.gene
      Weight.gene <- cbind(Weight.gene, weight.gene.temp)
      r.squared.full <- c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop <- rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) <- colnames(D1)
  rownames(Est.prop.weighted) <- colnames(Yjg)
  colnames(Est.prop.allgene) <- colnames(D1)
  rownames(Est.prop.allgene) <- colnames(Yjg)
  names(r.squared.full) <- colnames(Yjg)
  colnames(Weight.gene) <- colnames(Yjg)
  rownames(Weight.gene) <- cm_gene
  colnames(Var.prop) <- colnames(D1)
  rownames(Var.prop) <- colnames(Yjg)

  return(list(
    Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene,
    Weight.gene = Weight.gene, r.squared.full = r.squared.full, Var.prop = Var.prop
  ))
}
