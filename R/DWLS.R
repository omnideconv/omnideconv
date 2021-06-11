#' Signature matrix creation with DWLS using genes identified by de_analysis_mast()
#'
#' @param scdata A matrix or dataframe with the single-cell data. Rows are genes, columns are
#'   samples. Row and column names need to be set.
#' @param id A Vector of the cell type annotations. Has to be in the same order as the samples in
#'   single_cell_object
#' @param method The method used to create the signature matrix. Options are "mast" and "seurat"
#' @param path The path where the generated files will be saved. If path=NULL, the generated files
#'   will be discarded.
#' @param verbose Whether to output what DWLS is doing
#' @param diff_cutoff Cutoff to determine the FC-limit. How low can the lowest fold change be to
#'   still be considered differentially expressed?
#' @param pval_cutoff Cutoff to determine the pVal-limit. How high can the highest p-Value be to
#'   still be considered statistically significant?
#'
#' @return The signature matrix. Rows are genes, columns are cell types.
#' @export
#'
build_model_dwls <- function(scdata,
                             id,
                             method = c("mast", "seurat"),
                             path,
                             verbose = FALSE,
                             diff_cutoff = 0.5,
                             pval_cutoff = 0.01) {
  if (length(method) > 1) {
    method <- method[1]
  }

  if (method == "mast") {
    return(build_signature_matrix_using_mast(scdata, id, path, verbose, diff_cutoff, pval_cutoff))
  } else if (method == "seurat") {
    return(build_signature_matrix_using_Seurat(scdata, id, path, verbose, diff_cutoff, pval_cutoff))
  } else {
    base::stop("Could not find method " + method + ". Please try \"Mast\" or \"Seurat\"")
  }
}
#' Calculates the decomposition using the dwls algorithm
#'
#' generates a reference profile based on single-cell data. Learns a transformation of bulk
#' expression based on observed single-cell proportions and performs NNLS regression on these
#' transformed values to estimate cell proportions.
#'
#' @param bulk_gene_expression An Expression Set containing bulk data.
#' @param signature The Signature matrix.
#' @param dwls_submethod Three alternative methods in DWLS: OLS, SVR, and DampenedWLS.
#' @param verbose Whether the algorithm should print out what it is doing.

#'
#' @return A list. Slot bulk.props contains a matrix of cell type proportion estimates with cell
#'   types as rows and individuals as columns.
#' @export
#'

deconvolute_dwls <- function(bulk_gene_expression, signature,
                             dwls_submethod = c("DampenedWLS", "OLS", "SVR"), verbose = FALSE) {
  if (length(dwls_submethod) > 1) {
    dwls_submethod <- dwls_submethod[1]
  }

  if (verbose) {
    base::message("Running DWLS deconvolution module")
  }

  # trim data
  genes <- base::intersect(rownames(signature), rownames(bulk_gene_expression))
  bulk <- bulk_gene_expression[genes, , drop = FALSE]
  sig <- signature[genes, , drop = FALSE]
  if (class(bulk)[[1]] == "numeric" || class(sig)[[1]] == "numeric") {
    base::stop("Either bulk data or signature matrix just contains one row!")
  }

  # perform reconvolution in different sub_methods
  res <- NULL

  if (dwls_submethod == "OLS") {
    solutions_ols <- NULL
    for (i in 1:ncol(bulk)) {
      bulk_i <- bulk[, i]
      sol <- solve_ols(sig, bulk_i, verbose)
      # sol<-round(sol,5)
      solutions_ols <- cbind(solutions_ols, sol)
    }
    colnames(solutions_ols) <- colnames(bulk)
    res <- solutions_ols
  } else if (dwls_submethod == "SVR") {
    solutions_svr <- NULL
    for (i in 1:ncol(bulk)) {
      bulk_i <- bulk[, i]
      sol <- solve_svr(sig, bulk_i)
      # sol<-round(sol,5)
      solutions_svr <- cbind(solutions_svr, sol)
    }
    colnames(solutions_svr) <- colnames(bulk)
    res <- solutions_svr
  } else if (dwls_submethod == "DampenedWLS") {
    solutions_dampened_wls <- NULL
    for (i in 1:ncol(bulk)) {
      bulk_i <- bulk[, i]
      sol <- solve_dampened_wls(sig, bulk_i, verbose)
      # sol<-round(sol,5)
      solutions_dampened_wls <- cbind(solutions_dampened_wls, sol)
    }
    colnames(solutions_dampened_wls) <- colnames(bulk)
    res <- solutions_dampened_wls
  } else {
    base::stop("Submethod " + dwls_submethod + " not found. Please provide a valid one.")
  }

  if (verbose) {
    base::message("Deconvolution sucessful!")
  }
  result <- t(res)
  return(result)
}

#' Solver using OLS, constrained such that cell type numbers>0
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param verbose Whether to print output
#'
#' @return A vector with the cell type proportions for the sample
#'
solve_ols <- function(S, B, verbose = FALSE) {
  solution <- solve_ols_internal(S, B, verbose)
  return(solution / sum(solution))
}

#' Solver using OLS, constrained such that cell type numbers>0
#'
#' Returns absolute numbers, not proportions
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param verbose Whether to print output
#'
#' @return A vector with the cell type numbers for the sample
#'
solve_ols_internal <- function(S, B, verbose = FALSE) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- base::cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")

  solution <- tryCatch(
    {
      quadprog::solve.QP(D / sc, d / sc, A, bzero)$solution
    },
    error = function(cond) {
      base::message(
        "Error solving the quadratic programming problem. Your dataset might be too small. Run ",
        "with verbose=TRUE to get the full error message."
      )
      if (verbose) {
        base::stop(cond)
      } else {
        base::stop()
      }
    },
    warning = function(cond) {
      base::warning(cond)
    }
  )

  names(solution) <- colnames(S)
  return(solution)
}


#' Solver using WLS with weights dampened by a certain dampening constant
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param verbose Whether to print output
#'
#' @return A vector with the cell type numbers for the sample
#'
solve_dampened_wls <- function(S, B, verbose = FALSE) {
  # first solve OLS, use this solution to find a starting point for the weights
  solution <- solve_ols_internal(S, B, verbose)
  # now use dampened WLS, iterate weights until convergence
  iterations <- 0
  changes <- c()
  # find dampening constant for weights using cross-validation
  j <- find_dampening_constant(S, B, solution)
  change <- 1
  while (change > .01 & iterations < 1000) {
    newsolution <- solve_dampened_wlsj(S, B, solution, j, verbose)
    # decrease step size for convergence
    solution_average <-
      base::rowMeans(base::cbind(newsolution, matrix(
        solution,
        nrow = length(solution), ncol = 4
      )))
    change <- norm(as.matrix(solution_average - solution))
    solution <- solution_average
    iterations <- iterations + 1
    changes <- c(changes, change)
  }
  # print(round(solution/sum(solution),5))
  return(solution / sum(solution))
}


#' Solve WLS given a dampening constant
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param gold_standard The average of all the solutions so far
#' @param j The dampening constant
#' @param verbose Whether to print output
#'
#' @return A vector with the cell type numbers for the sample
#'
solve_dampened_wlsj <-
  function(S, B, gold_standard, j, verbose = FALSE) {
    multiplier <- 1 * 2^(j - 1)
    sol <- gold_standard
    ws <- as.vector((1 / (S %*% sol))^2)
    ws_scaled <- ws / min(ws)
    ws_dampened <- ws_scaled
    ws_dampened[base::which(ws_scaled > multiplier)] <- multiplier
    W <- base::diag(ws_dampened)
    D <- t(S) %*% W %*% S
    d <- t(S) %*% W %*% B
    A <- base::cbind(diag(dim(S)[2]))
    bzero <- c(rep(0, dim(S)[2]))
    sc <- norm(D, "2")

    solution <- tryCatch(
      {
        quadprog::solve.QP(D / sc, d / sc, A, bzero)$solution
      },
      error = function(cond) {
        base::message(
          "Error solving the quadratic programming problem. Your dataset might be too small. Run ",
          "with verbose=TRUE to get the full error message."
        )
        if (verbose) {
          base::stop(cond)
        } else {
          base::stop()
        }
      },
      warning = function(cond) {
        base::warning(cond)
      }
    )

    names(solution) <- colnames(S)
    return(solution)
  }

#' Finding a dampening constant for the weights using cross-validation
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#' @param gold_standard The solution found with OLS
#'
#' @return The dampening constant (integer)
find_dampening_constant <- function(S, B, gold_standard) {
  solutions_sd <- NULL
  # gold_standard is used to define the weights
  sol <- gold_standard
  ws <- as.vector((1 / (S %*% sol))^2)
  ws_scaled <- ws / min(ws)
  ws_scaled_minus_inf <- ws_scaled
  # ignore infinite weights
  if (max(ws_scaled) == "Inf") {
    ws_scaled_minus_inf <- ws_scaled[-base::which(ws_scaled == "Inf")]
  }
  # try multiple values of the dampening constant (multiplier)
  # for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:base::ceiling(log2(max(ws_scaled_minus_inf)))) {
    multiplier <- 1 * 2^(j - 1)
    ws_dampened <- ws_scaled
    ws_dampened[which(ws_scaled > multiplier)] <- multiplier
    solutions <- NULL
    seeds <- c(1:100)

    for (i in 1:100) {
      base::set.seed(seeds[i]) # make nondeterministic
      subset <-
        sample(length(ws), size = length(ws) * 0.5) # randomly select half of gene set
      # solve dampened weighted least squares for subset
      fit <- stats::lm(B[subset] ~ -1 + S[subset, , drop = FALSE], weights = ws_dampened[subset])
      sol <- fit$coef * sum(gold_standard) / sum(fit$coef)
      solutions <- base::cbind(solutions, sol)
    }

    solutions_sd <-
      base::cbind(solutions_sd, base::apply(solutions, 1, stats::sd))
  }
  # choose dampening constant that results in least cross-validation variance
  j <- which.min(base::colMeans(solutions_sd^2))
  return(j)
}


#' Solver using SVR
#'
#' @param S Signature matrix
#' @param B One column of the bulk data matrix, so one sample
#'
#' @return A vector with the cell type proportions for the sample
#'
solve_svr <- function(S, B) {
  # scaling
  ub <- max(c(as.vector(S), B)) # upper bound
  lb <- min(c(as.vector(S), B)) # lower bound
  Bs <- ((B - lb) / ub) * 2 - 1
  Ss <- ((S - lb) / ub) * 2 - 1

  # perform SVR
  model <-
    e1071::svm(
      Ss,
      Bs,
      nu = 0.5,
      scale = TRUE,
      type = "nu-regression",
      kernel = "linear",
      cost = 1
    )
  coef <- t(model$coefs) %*% model$SV
  coef[base::which(coef < 0)] <- 0
  coef <- as.vector(coef)
  names(coef) <- colnames(S)
  # print(round(coef/sum(coef),5))
  return(coef / sum(coef))
}


#' Performing DE analysis using Seurat
#'
#' When path = NULL, the generated files in the processes will not be saved and output
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#'
#' @return List with the differentially expressed genes for each cell type
#'
de_analysis <- function(scdata, id, path = NULL) {
  list_names <- unique(id)
  list_de_group <- as.list(rep(0, length(list_names)))

  expr_obj <-
    Seurat::CreateSeuratObject(raw.data = as.data.frame(scdata), project = "DE")
  expr_obj2 <- Seurat::SetIdent(expr_obj, ident.use = as.vector(id))
  # print("Calculating differentially expressed genes:")
  for (i in unique(id)) {
    de_group <-
      Seurat::FindMarkers(
        object = expr_obj2,
        ident.1 = i,
        ident.2 = NULL,
        only.pos = TRUE,
        test.use = "bimod"
      )

    index <- which(list_names == i)
    list_de_group[[index]] <- de_group

    if (!is.null(path)) {
      save(de_group, file = base::paste(path, "/de_", i, ".RData", sep = ""))
    }
  }
  return(list_de_group)
}


#' Building the signature matrix using Seurat
#'
#' When path = NULL, the generated files in the processes will not be saved and output.
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#' @param verbose Whether to output what DWLS is doing
#' @param diff_cutoff The FC cutoff
#' @param pval_cutoff The pValue cutoff
#'
#' @return The computed signature matrix
#'
build_signature_matrix_using_Seurat <- function(scdata,
                                                id,
                                                path = NULL,
                                                verbose = FALSE,
                                                diff_cutoff = 0.5,
                                                pval_cutoff = 0.01) {
  # perform differential expression analysis
  list_de_groups <- de_analysis(scdata, id, path)

  number_of_genes <- c()
  for (i in unique(id)) {
    name_list <- unique(id)
    index <- which(name_list == i)
    de_group <- list_de_groups[[index]]

    de_genes <-
      rownames(de_group)[base::intersect(
        which(de_group$p_val_adj < pval_cutoff),
        which(de_group$avg_logFC > diff_cutoff)
      )]
    non_mir <- base::grep("MIR|Mir", de_genes, invert = TRUE)
    base::assign(
      base::paste("cluster_lr_test_table.", i, sep = ""),
      de_group[which(rownames(de_group) %in% de_genes[non_mir]), ]
    )
    number_of_genes <- c(number_of_genes, length(de_genes[non_mir]))
  }

  # need to reduce number of genes
  # for each subset, order significant genes by decreasing fold change,
  # choose between 50 and 200 genes
  # choose matrix with lowest condition number
  condition_numbers <- c()
  for (g in 50:200) {
    genes <- c()
    j <- 1
    for (i in unique(id)) {
      if (number_of_genes[j] > 0) {
        temp <- base::paste("cluster_lr_test_table.", i, sep = "")
        temp <- get(temp)
        temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
        genes <-
          c(genes, (rownames(temp)[1:min(g, number_of_genes[j])]))
      }
      j <- j + 1
    }
    genes <- unique(genes)
    # make signature matrix
    expr_subset <- scdata[genes, , drop = FALSE]
    sig <- NULL
    for (i in unique(id)) {
      sig <-
        base::cbind(sig, (apply(expr_subset, 1, function(y) {
          mean(y[which(id == i)])
        })))
    }
    colnames(sig) <- unique(id)
    condition_numbers <- c(condition_numbers, kappa(sig))
  }

  # g is optimal gene number
  g <- which.min(condition_numbers) + min(49, number_of_genes - 1)
  genes <- c()
  j <- 1
  for (i in unique(id)) {
    if (number_of_genes[j] > 0) {
      temp <- base::paste("cluster_lr_test_table.", i, sep = "")
      temp <- get(temp)
      temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
      genes <-
        c(genes, (rownames(temp)[1:min(g, number_of_genes[j])]))
    }
    j <- j + 1
  }
  genes <- unique(genes)
  expr_subset <- scdata[genes, , drop = FALSE]
  sig <- NULL
  for (i in unique(id)) {
    sig <-
      base::cbind(sig, (apply(expr_subset, 1, function(y) {
        mean(y[which(id == i)])
      })))
  }

  colnames(sig) <- unique(id)
  rownames(sig) <- genes

  if (!is.null(path)) {
    save(sig, file = base::paste(path, "/sig.RData", sep = ""))
  }

  return(sig)
}


# Functions needed for DE

mean_in_log2space <- function(x, pseudo_count) {
  return(log2(mean(2^(x) - pseudo_count) + pseudo_count))
}

stat_log2 <- function(data_m, group_v, pseudo_count) {
  # data_m=data_used_log2
  log2_mean_r <-
    stats::aggregate(t(data_m), list(as.character(group_v)), function(x) {
      mean_in_log2space(x, pseudo_count)
    })
  log2_mean_r <- t(log2_mean_r)
  colnames(log2_mean_r) <-
    base::paste("mean_group", log2_mean_r[1, ], sep = "")
  log2_mean_r <- log2_mean_r[-1, ]
  log2_mean_r <- as.data.frame(log2_mean_r)
  log2_mean_r <- varhandle::unfactor(log2_mean_r) # from varhandle
  log2_mean_r[, 1] <- as.numeric(log2_mean_r[, 1])
  log2_mean_r[, 2] <- as.numeric(log2_mean_r[, 2])
  log2_foldchange <- log2_mean_r$mean_group1 - log2_mean_r$mean_group0
  results <- data.frame(base::cbind(
    log2_mean_r$mean_group0,
    log2_mean_r$mean_group1,
    log2_foldchange
  ))
  colnames(results) <- c("log2_mean_group0", "log2_mean_group1", "log2_fc")
  rownames(results) <- rownames(log2_mean_r)
  return(results)
}

v_auc <- function(data_v, group_v) {
  prediction_use <- ROCR::prediction(data_v, group_v, 0:1)
  perf_use <- ROCR::performance(prediction_use, "auc")
  auc_use <- round(perf_use@y.values[[1]], 3)
  return(auc_use)
}

m_auc <- function(data_m, group_v) {
  AUC <- apply(data_m, 1, function(x) {
    v_auc(x, group_v)
  })
  AUC[is.na(AUC)] <- 0.5
  return(AUC)
}


#' Performing DE analysis using mast
#'
#' When path = NULL, the generated files in the processes will not be saved and output.
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#' @param verbose Whether to print output
#'
#' @return A list with the cell types and their differentially expressed genes
#'
de_analysis_mast <- function(scdata, id, path, verbose = FALSE) {
  list_names <- unique(id)
  list_lr_test_table <- as.list(rep(0, length(list_names)))

  pseudo_count <- 0.1
  data_used_log2 <- log2(scdata + pseudo_count)
  colnames(data_used_log2) <- make.unique(colnames(data_used_log2))
  diff_cutoff <- 0.5
  for (i in unique(id)) {
    cells_symbol_list2 <- colnames(data_used_log2)[which(id == i)]
    cells_coord_list2 <- match(cells_symbol_list2, colnames(data_used_log2))
    cells_symbol_list1 <- colnames(data_used_log2)[which(id != i)]
    cells_coord_list1 <- match(cells_symbol_list1, colnames(data_used_log2))
    data_used_log2_ordered <-
      base::cbind(data_used_log2[, cells_coord_list1], data_used_log2[, cells_coord_list2])
    group_v <-
      c(rep(0, length(cells_coord_list1)), rep(1, length(cells_coord_list2)))
    # ouput
    log2_stat_result <-
      stat_log2(data_used_log2_ordered, group_v, pseudo_count)
    auc <- m_auc(data_used_log2_ordered, group_v)
    bigtable <- data.frame(base::cbind(log2_stat_result, auc))

    de <- bigtable[bigtable$log2_fc > diff_cutoff, ]
    dim(de)
    if (dim(de)[1] > 1) {
      data_1 <- data_used_log2[, cells_coord_list1, drop = FALSE]
      data_2 <- data_used_log2[, cells_coord_list2, drop = FALSE]
      genes_list <- rownames(de)
      log2fold_change <- base::cbind(genes_list, de$log2_fc)
      colnames(log2fold_change) <- c("gene_name", "log2fold_change")
      counts <- as.data.frame(base::cbind(data_1[genes_list, ], data_2[genes_list, ]))
      groups <- c(
        rep("Cluster_Other", length(cells_coord_list1)),
        rep(i, length(cells_coord_list2))
      )
      groups <- as.character(groups)
      data_for_mist <-
        verbose_wrapper(verbose)(as.data.frame(base::cbind(
          rep(rownames(counts), dim(counts)[2]),
          reshape::melt(counts),
          rep(groups, each = dim(counts)[1]),
          rep(1, dim(counts)[1] * dim(counts)[2])
        )))
      colnames(data_for_mist) <- c(
        "gene",
        "Subject.ID",
        "Et",
        "Population",
        "Number.of.Cells"
      )
      vbeta <- data_for_mist
      vbeta_fa <-
        verbose_wrapper(verbose)(
          MAST::FromFlatDF(
            vbeta,
            idvars = c("Subject.ID"),
            primerid = "gene",
            measurement = "Et",
            ncells = "Number.of.Cells",
            geneid = "gene",
            cellvars = c("Number.of.Cells", "Population"),
            phenovars = c("Population"),
            id = "vbeta all"
          )
        )
      vbeta_1 <- subset(vbeta_fa, Number.of.Cells == 1)
      # .3 mast
      utils::head(SummarizedExperiment::colData(vbeta_1))
      zlm_output <-
        verbose_wrapper(verbose)(MAST::zlm(
          ~Population,
          vbeta_1,
          method = "bayesglm",
          ebayes = TRUE
        ))
      if (verbose) {
        methods::show(zlm_output)
      }
      coef_and_c_i <- summary(zlm_output, logFC = TRUE)
      zlm_lr <-
        verbose_wrapper(verbose)(MAST::lrTest(zlm_output, "Population"))
      zlm_lr_pvalue <- reshape::melt(zlm_lr[, , "Pr(>Chisq)"])
      zlm_lr_pvalue <-
        zlm_lr_pvalue[which(zlm_lr_pvalue$test.type == "hurdle"), ]



      lr_test_table <-
        merge(zlm_lr_pvalue, de, by.x = "primerid", by.y = "row.names")
      colnames(lr_test_table) <-
        c(
          "gene",
          "test.type",
          "p_value",
          base::paste("log2_mean_", "Cluster_Other", sep = ""),
          base::paste("log2_mean_", i, sep = ""),
          "log2fold_change",
          "auc"
        )
      cluster_lr_test_table <-
        lr_test_table[rev(order(lr_test_table$auc)), ]

      # . 4 save results

      index <- which(list_names == i)
      list_lr_test_table[[index]] <- cluster_lr_test_table

      if (!is.null(path)) {
        utils::write.csv(cluster_lr_test_table,
          file = base::paste(path, "/", i, "_lr_test.csv", sep = "")
        )
        save(cluster_lr_test_table,
          file = base::paste(path, "/", i, "_mist.RData", sep = "")
        )
      }
    }
  }

  return(list_lr_test_table)
}


#' #' Building the signature matrix using mast
#'
#' When path = NULL, the generated files in the processes will not be saved and output.
#'
#' @param scdata The single cell data matrix
#' @param id A Vector of the cell type annotations
#' @param path OPTIONAL path for saving generated files
#' @param verbose Whether to print output
#' @param diff_cutoff The FC cutoff
#' @param pval_cutoff The pValue cutoff
#'
#' @return The computed signature matrix
#'
build_signature_matrix_using_mast <- function(scdata,
                                              id,
                                              path,
                                              verbose = FALSE,
                                              diff_cutoff = 0.5,
                                              pval_cutoff = 0.01) {
  # compute differentially expressed genes for each cell type
  list_cluster_table <-
    de_analysis_mast(scdata, id, path, verbose = verbose)

  # for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the
  # estimated fold-change is greater than 0.5
  number_of_genes <- c()
  for (i in unique(id)) {
    name_list <- unique(id)
    index <- which(name_list == i)
    cluster_lr_test_table <- list_cluster_table[[index]]
    pvalue_adjusted <-
      stats::p.adjust(
        cluster_lr_test_table[, 3],
        method = "fdr",
        n = length(cluster_lr_test_table[, 3])
      )
    cluster_lr_test_table <-
      base::cbind(cluster_lr_test_table, pvalue_adjusted)
    de_genes <-
      cluster_lr_test_table$gene[base::intersect(
        which(pvalue_adjusted < pval_cutoff),
        which(cluster_lr_test_table$log2fold_change > diff_cutoff)
      )]

    # because Mir gene is usually not accurate
    non_mir <- base::grep("MIR|Mir", de_genes, invert = TRUE)
    base::assign(
      base::paste("cluster_lr_test_table.", i, sep = ""),
      cluster_lr_test_table[which(cluster_lr_test_table$gene %in% de_genes[non_mir]), ]
    )
    number_of_genes <- c(number_of_genes, length(de_genes[non_mir]))
  }


  # need to reduce number of genes
  # for each subset, order significant genes by decreasing fold change,
  # choose between 50 and 200 genes
  # for each, iterate and choose matrix with lowest condition number
  condition_numbers <- c()
  for (g in 50:200) {
    genes <- c()
    j <- 1
    for (i in unique(id)) {
      if (number_of_genes[j] > 0) {
        temp <- base::paste("cluster_lr_test_table.", i, sep = "")
        temp <- get(temp)
        temp <-
          temp[order(temp$log2fold_change, decreasing = TRUE), ]
        genes <-
          c(genes, varhandle::unfactor(temp$gene[1:min(g, number_of_genes[j])]))
      }
      j <- j + 1
    }
    genes <- unique(genes)
    # make signature matrix
    expr_subset <- scdata[genes, , drop = FALSE]
    sig <- NULL
    for (i in unique(id)) {
      sig <-
        base::cbind(sig, (apply(expr_subset, 1, function(y) {
          mean(y[which(id == i)])
        })))
    }
    colnames(sig) <- unique(id)
    condition_numbers <- c(condition_numbers, kappa(sig))
  }
  # g is optimal gene number
  g <- which.min(condition_numbers) + min(49, number_of_genes - 1)
  genes <- c()
  j <- 1
  for (i in unique(id)) {
    if (number_of_genes[j] > 0) {
      temp <- base::paste("cluster_lr_test_table.", i, sep = "")
      temp <- get(temp)
      temp <- temp[order(temp$log2fold_change, decreasing = TRUE), ]
      genes <-
        c(genes, varhandle::unfactor(temp$gene[1:min(g, number_of_genes[j])]))
    }
    j <- j + 1
  }
  genes <- unique(genes)
  expr_subset <- scdata[genes, , drop = FALSE]
  sig <- NULL
  for (i in unique(id)) {
    sig <-
      base::cbind(sig, (apply(expr_subset, 1, function(y) {
        mean(y[which(id == i)])
      })))
  }

  colnames(sig) <- unique(id)
  rownames(sig) <- genes

  if (!is.null(path)) {
    save(sig, file = base::paste(path, "/sig.RData", sep = ""))
  }


  return(sig)
}
