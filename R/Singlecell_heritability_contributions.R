
#' get_Pathway_sclm
#'
#' @param pa_block single pathway blocks.
#' @param pca_scCell_mat Element in Pagwas.
#' @param data_mat Element in Pagwas.
#' @param rawPathway_list Element in Pagwas..
#' @param n.cores cores for regression
#' @param backingpath file address for bk files, no "/" in the end.
#' @param snp_gene_df Element in Pagwas.
#' @param Rns Random bk names
#' @return lm result for pathway in single cell.
get_Pathway_sclm <- function(pa_block,
                             pca_scCell_mat,
                             data_mat,
                             rawPathway_list,
                             snp_gene_df,
                             n.cores=1,
                             backingpath,
                             Rns) {

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

  x2 <- as(x2, "CsparseMatrix")

  if (pa_block$n_snps > 1) {
    x2 <- x2[pa_block$snps$label, ]
    pa_block$n_snps <- nrow(pa_block$snps)

    x <- x[rep(1, pa_block$n_snps), ]
    rownames(x) <- pa_block$snps$rsid

    #x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
    x2 <- x2 * x
  } else {
    x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
    rownames(x2) <- pa_block$snps$label
    pa_block$n_snps <- nrow(pa_block$snps)

    x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
    rownames(x) <- pa_block$snps$rsid


    #x <- matrix(as.numeric(x) * as.numeric(
    #  snp_gene_df[pa_block$snps$rsid, "slope"]
    #), nrow = 1)
    x2 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
    x2 <- as(x2, "CsparseMatrix")
  }
  pa_block$x <- as(pa_block$ld_matrix_squared %*% x2 , "CsparseMatrix")

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
        Rns=Rns,
        backingpath=backingpath,
        n.cores=n.cores
      )
      results[is.na(results)] <- 0
      names(results) <- colnames(data_mat)
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
#' @return Pagwas result list including scPagwas.gPAS.score
#' @export

scPagwas_perform_score <- function(Pagwas,
                                   remove_outlier = TRUE) {
  options(bigmemory.allow.dimnames = TRUE)

  Pathway_sclm_results <- Pagwas$Pathway_sclm_results

  Pathway_names <- colnames(Pathway_sclm_results)

  pathway_expr <- lapply(Pathway_names, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], rownames(Pagwas$data_mat))
    if (length(a) == 0) {
      return(NULL)
    } else if (length(a) == 1) {
      return(Pagwas$data_mat[a, ])
    } else {
      b <- biganalytics::apply(Pagwas$data_mat[a, ], 2, mean)
      return(b)
    }
  })
  a<-!sapply(pathway_expr, is.null)
  pathway_expr <- data.matrix(as.data.frame(pathway_expr[a]))
  Pathway_names <- Pathway_names[a]

  colnames(pathway_expr) <- Pathway_names
  pa_exp_mat <- t(Pagwas$pca_scCell_mat[Pathway_names, ]) * pathway_expr
  rm(pathway_expr)
  Pagwas$Pathway_single_results <- Pathway_sclm_results[, Pathway_names] * pa_exp_mat
  rownames(Pagwas$Pathway_single_results) <- colnames(Pagwas$pca_scCell_mat)
  message("* Get Pathways'rankPvalue for each celltypes!")
  cl <- unique((Pagwas$Celltype_anno$annotation))

  Pagwas$Pathway_single_results <- t(data.matrix(
    Pagwas$Pathway_single_results
  ))

  Pathways_rankPvalue <- lapply(cl, function(ss) {
    print(ss)
    tt <- Pagwas$Celltype_anno$annotation == ss
    PathwayrankPvalue <- scGene_rankP(Pagwas$Pathway_single_results[, tt])
    return(PathwayrankPvalue$pValueHigh)
  })

  Pagwas$scPathways_rankPvalue <- Reduce(
    function(dtf1, dtf2) cbind(dtf1, dtf2),
    Pathways_rankPvalue
  )
  Pagwas$scPathways_rankPvalue <- as.data.frame(Pagwas$scPathways_rankPvalue)
  colnames(Pagwas$scPathways_rankPvalue) <- cl
  rownames(Pagwas$scPathways_rankPvalue) <- rownames(Pagwas$Pathway_single_results)
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
#'
#' @param Pagwas_x x parameter for lm
#' @param Pagwas_y y parameter for lm
#' @param noise_per_snp noise
#' @param n.cores cores for regression
#' @param Rns the names for bk temfiles
#' @param backingpath file address for bk files, no "/" in the end.
#'
#' @export


scParameter_regression <- function(Pagwas_x,
                                   Pagwas_y,
                                   noise_per_snp,
                                   Rns,
                                   n.cores=1,
                                   backingpath) {

    backingpath<- paste0(backingpath,"/",Rns)

       #Pagwas_x<-(Pagwas_x)
      liear_m <- bigstatsr::big_univLinReg(
        X = bigstatsr::as_FBM(as_matrix(Pagwas_x),backingfile = backingpath),
        y.train = Pagwas_y,
        covar.train = bigstatsr::covar_from_df(
          data.frame(offset(noise_per_snp))
        ),
        ncores = n.cores
      )
      parameters <- liear_m$estim
      unlink(paste0(backingpath,".bk"),recursive = TRUE)
   # rm(Pagwas_x)
  return(parameters)
}

#' scPagwas_score_filter
#' @description filter the cPagwas_score for outliers.
#' @param scPagwas_score (data.frame)
#'

scPagwas_score_filter <- function(scPagwas_score) {
  # remove the NAN!
  if (sum(is.nan(scPagwas_score)) > 0) {
    scPagwas_score[is.nan(scPagwas_score)] <- 0
  }
  # remove the inf values!
  if (Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == Inf)] <-
      max(scPagwas_score[-which(scPagwas_score == Inf)], na.rm = TRUE)
  }
  if (-Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == -Inf)] <-
      min(scPagwas_score[-which(scPagwas_score == -Inf)], na.rm = TRUE)
  }
  lower_bound <- stats::quantile(scPagwas_score, 0.01, na.rm = TRUE)
  upper_bound <- stats::quantile(scPagwas_score, 0.99, na.rm = TRUE)

  lower_ind <- which(scPagwas_score < lower_bound)
  upper_ind <- which(scPagwas_score > upper_bound)
  scPagwas_score[lower_ind] <- lower_bound
  scPagwas_score[upper_ind] <- upper_bound

  return(scPagwas_score)
}



#' scGet_gene_heritability_contributions
#'
#' @param scPagwas.gPAS.score score of scPagwas pathway lm
#' @param data_mat matrix for single cell data
#' @return result list including gene heritability correlation
#' @export

scGet_gene_heritability_correlation <- function(scPagwas.gPAS.score,data_mat) {

  if (all(names(scPagwas.gPAS.score) == colnames(data_mat))) {
    scPagwas.gPAS.score <- data.matrix(scPagwas.gPAS.score)
    sparse_cor <- corSparse(
      X = t(as_matrix(data_mat)),
      Y = scPagwas.gPAS.score
    )
  } else {
    scPagwas.gPAS.score <- data.matrix(scPagwas.gPAS.score)
    sparse_cor <- corSparse(
      X = t(as_matrix(data_mat)),
      Y = scPagwas.gPAS.score
    )
  }
  rownames(sparse_cor) <- rownames(data_mat)
  colnames(sparse_cor) <- "gene_heritability_correlation"
  sparse_cor[is.nan(sparse_cor)] <- 0
  #
  return(sparse_cor)
}
