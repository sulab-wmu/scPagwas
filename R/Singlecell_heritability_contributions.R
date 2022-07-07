
#' get_Pathway_sclm
#'
#' @param Pathway_ld_gwas_data Element in Pagwas
#' @param pca_scCell_mat Element in Pagwas
#' @param data_mat Element in Pagwas
#' @param rawPathway_list Element in Pagwas
#' @param snp_gene_df Element in Pagwas
#' @param split_n number of times to compute the singlecell result
#'
#' @return
get_Pathway_sclm <- function(pa_block,
                             pca_scCell_mat,
                             data_mat,
                             ncores = 1,
                             rawPathway_list,
                             snp_gene_df) {
  pathway <- unique(pa_block$block_info$pathway)
  x <- matrix(pca_scCell_mat[pathway, ], nrow = 1)
  rownames(x) <- pathway

  if (pa_block$n_snps == 0) {
    pa_block$include_in_inference <- F
    pa_block$x <- NULL # to make sure we totally replace previous stuffs
    return(pa_block)
  }

  mg <- intersect(rawPathway_list[[pathway]], rownames(data_mat))
  if (length(mg) == 1) {
    x2 <- matrix(data_mat[mg, ], nrow = 1)
    x2 <- x2 / (x2 + 0.0001)
    rownames(x2) <- mg
  } else {
    x2 <- biganalytics::apply(data_mat[mg, ], 2, function(ge) {
      if (sum(ge) == 0) {
        return(rep(0, length(ge)))
      } else {
        return(ge / sum(ge))
      }
    })
    rownames(x2) <- mg
  }

  x2 <- as(x2, "dgCMatrix")

  if (pa_block$n_snps > 1) {
    x2 <- x2[pa_block$snps$label, ]
    pa_block$n_snps <- nrow(pa_block$snps)

    x <- x[rep(1, pa_block$n_snps), ]
    rownames(x) <- pa_block$snps$rsid
    rownames(snp_gene_df) <- snp_gene_df$rsid
    x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
    x2 <- x2 * x
  } else {
    x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
    rownames(x2) <- pa_block$snps$label
    pa_block$n_snps <- nrow(pa_block$snps)

    x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
    rownames(x) <- pa_block$snps$rsid

    rownames(snp_gene_df) <- snp_gene_df$rsid
    x <- matrix(as.numeric(x) * as.numeric(snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
    x2 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
    x2 <- as(x2, "dgCMatrix")
  }
  pa_block$x <- as(pa_block$ld_matrix_squared %*% x2, "dgCMatrix")

  pa_block$include_in_inference <- T
  noise_per_snp <- pa_block$snps$se**2

  if (!is.null(pa_block$x)) {
    if (pa_block$n_snps > 2) {
      na_elements <- is.na(pa_block$y) | apply(pa_block$x, 1, function(x) {
        any(is.na(x))
      }) | is.na(noise_per_snp)

      results <- scParameter_regression(
        Pagwas_x = pa_block$x[!na_elements, ],
        Pagwas_y = pa_block$y[!na_elements],
        noise_per_snp = noise_per_snp[!na_elements],
        ncores = ncores
      )


      results[is.na(results)] <- 0
      names(results) <- colnames(data_mat)
      # Pathway_cell_regression<-results
    } else {
      results <- NULL
    }
  } else {
    results <- NULL
  }

  return(results)
}

#' scPagwas_perform_score
#' @description Get the scPagwas score for each cells
#' @param Pagwas pagwas
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' storage, you can set a large number to save time.
#' @return
#' @export
#' @examples
#' library(scPagwas)
#' scPagwas_perform_score(Pagwas_data)
scPagwas_perform_score <- function(Pagwas,
                                   remove_outlier = TRUE) {
  options(bigmemory.allow.dimnames = TRUE)
  # Pathway_sclm_results<-
  # fit model
  # names(Pathway_sclm_results) <- Pathway_names
  Pathway_sclm_results <- Pagwas$Pathway_sclm_results # [!sapply(Pagwas$Pathway_sclm_results, is.null)]
  Pathway_names <- colnames(Pathway_sclm_results)
  # Pathway_sclm_results <- #data.matrix(as.data.frame(Pathway_sclm_results))
  ncells <- nrow(Pathway_sclm_results)

  pathway_expr <- lapply(Pathway_names, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], rownames(Pagwas$data_mat))
    if (length(a) == 0) {
      # return(rep(0, ncol(Pagwas$data_mat)))
      return(NULL)
    } else if (length(a) == 1) {
      return(Pagwas$data_mat[a, ])
      # return(NULL)
    } else {
      b <- biganalytics::apply(Pagwas$data_mat[a, ], 2, mean)
      return(b)
    }
  })
  pathway_expr <- pathway_expr[!sapply(pathway_expr, is.null)]
  Pathway_names <- Pathway_names[!sapply(pathway_expr, is.null)]
  # rm(data_mat)
  pathway_expr <- data.matrix(as.data.frame(pathway_expr))

  # pathway_expr<-as(pathway_expr,"dgCMatrix")
  colnames(pathway_expr) <- Pathway_names
  pa_exp_mat <- t(Pagwas$pca_scCell_mat[Pathway_names, ]) * pathway_expr
  # pa_exp_mat <-
  rm(pathway_expr)

  # pa_exp_mat<-as(pa_exp_mat,"dgCMatrix")
  Pagwas$Pathway_single_results <- Pathway_sclm_results[, Pathway_names] * pa_exp_mat

  rownames(Pagwas$Pathway_single_results) <- colnames(Pagwas$pca_scCell_mat)

  message("* Get Pathways'rankPvalue for each celltypes!")
  cl <- unique((Pagwas$Celltype_anno$annotation))

  Pagwas$Pathway_single_results <- t(data.matrix(Pagwas$Pathway_single_results))
  Pathways_rankPvalue <- lapply(cl, function(ss) {
    tt <- Pagwas$Celltype_anno$annotation == ss
    PathwayrankPvalue <- rankPvalue(Pagwas$Pathway_single_results[, tt])
    return(PathwayrankPvalue$pValueHigh)
  })
  Pagwas$scPathways_rankPvalue <- Reduce(function(dtf1, dtf2) cbind(dtf1, dtf2), Pathways_rankPvalue)
  Pagwas$scPathways_rankPvalue <- as.data.frame(Pagwas$scPathways_rankPvalue)
  colnames(Pagwas$scPathways_rankPvalue) <- cl
  rownames(Pagwas$scPathways_rankPvalue) <- colnames(Pathway_sclm_results)
  rm(Pathways_rankPvalue)
  message("* Get scPgwas score for each single cell")
  scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
  gc()
  if (remove_outlier) {
    Pagwas$scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
  }
  names(Pagwas$scPagwas.gPAS.score) <- colnames(Pagwas$pca_scCell_mat)

  return(Pagwas)
}



#' scParameter_regression
#' @description Find parameter estimates for the data.
#' @param Pagwas_x x parameter for lm
#' @param Pagwas_y y parameter for lm
#' @param noise_per_snp noise
#' @export
#' @return

scParameter_regression <- function(Pagwas_x,
                                   Pagwas_y,
                                   noise_per_snp,
                                   ncores = 1) {

  # x_df <-Pagwas_x

  if (bigmemory::is.big.matrix(Pagwas_x)) {
    #  liear_m<- biglm.big.matrix(Pagwas_y ~ offset(noise_per_snp) + Pagwas_x)
    liear_m <- bigstatsr::big_univLinReg(
      X = as_FBM(Pagwas_x[]),
      y.train = Pagwas_y,
      covar.train = bigstatsr::covar_from_df(data.frame(offset(noise_per_snp))),
      ncores = ncores
    )
    parameters <- liear_m$estim
  } else if (is(Pagwas_x, "dgCMatrix")) {
    Pagwas_x <- as_matrix(Pagwas_x)
    liear_m <- bigstatsr::big_univLinReg(
      X = as_FBM(Pagwas_x),
      y.train = Pagwas_y,
      covar.train = bigstatsr::covar_from_df(data.frame(offset(noise_per_snp))),
      ncores = ncores
    )
    parameters <- liear_m$estim
  } else {
    # Pagwas_x<-as_matrix(Pagwas_x)
    liear_m <- bigstatsr::big_univLinReg(
      X = as_FBM(Pagwas_x),
      y.train = Pagwas_y,
      covar.train = bigstatsr::covar_from_df(data.frame(offset(noise_per_snp))),
      ncores = ncores
    )
    parameters <- liear_m$estim
  }

  return(parameters)
}


#' scPagwas_score_filter
#' @description filter the cPagwas_score for outliers.
#' @param scPagwas_score (data.frame)
#'
#' @return

scPagwas_score_filter <- function(scPagwas_score) {
  # scPagwas_score <- scPagwas_score$scPagwas_score
  # remove the NAN!
  if (sum(is.nan(scPagwas_score)) > 0) {
    scPagwas_score[is.nan(scPagwas_score)] <- 0
  }
  # remove the inf values!
  if (Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == Inf)] <- max(scPagwas_score[-which(scPagwas_score == Inf)], na.rm = TRUE)
  }
  if (-Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == -Inf)] <- min(scPagwas_score[-which(scPagwas_score == -Inf)], na.rm = TRUE)
  }
  lower_bound <- quantile(scPagwas_score, 0.01, na.rm = TRUE)
  upper_bound <- quantile(scPagwas_score, 0.99, na.rm = TRUE)

  lower_ind <- which(scPagwas_score < lower_bound)
  upper_ind <- which(scPagwas_score > upper_bound)
  scPagwas_score[lower_ind] <- lower_bound
  scPagwas_score[upper_ind] <- upper_bound

  return(scPagwas_score)
}



#' scGet_gene_heritability_contributions
#'
#' @param Pagwas result of scPagwas
#'
#' @return
#' @export
#' @examples
#' Single_data <- readRDS("E:/RPakage/scPagwas/inst/extdata/scRNAexample.rds")
#' Pagwas$sparse_cor <- scGet_gene_heritability_contributions(
#'   scPagwas_score = Pagwas$scPagwas_score,
#'   data_mat = Seurat::GetAssayData(object = Single_data[, names(Pagwas$scPagwas_score)], slot = "data")
#' )
scGet_gene_heritability_correlation <- function(Pagwas) {
  if (all(names(Pagwas$scPagwas.gPAS.score) == colnames(Pagwas$data_mat))) {
    scPagwas.gPAS.score <- data.matrix(Pagwas$scPagwas.gPAS.score)
    sparse_cor <- corSparse(
      X = t(as_matrix(Pagwas$data_mat)),
      Y = scPagwas.gPAS.score
    )
  } else {
    data_mat <- Pagwas$data_mat[, names(scPagwas.gPAS.score)]
    scPagwas.gPAS.score <- data.matrix(scPagwas.gPAS.score)
    sparse_cor <- corSparse(
      X = t(as_matrix(Pagwas$data_mat)),
      Y = scPagwas.gPAS.score
    )
  }
  rownames(sparse_cor) <- rownames(Pagwas$data_mat)
  colnames(sparse_cor) <- "gene_heritability_correlation"
  sparse_cor[is.nan(sparse_cor)] <- 0
  Pagwas$gene_heritability_correlation <- sparse_cor

  return(Pagwas)
}
