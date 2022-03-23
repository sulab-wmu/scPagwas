
#' scPagwas_perform_score
#' @description Get the scPagwas score for each cells
#' @param Pagwas Pagwas data list
#' @param n.cores cores
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' scPagwas_perform_score(Pagwas)
scPagwas_perform_score <- function(Pagwas,
                                   n.cores = 1,
                                   remove_outlier=TRUE) {
  if (is.null(Pathway_ld_gwas_data)) {
    stop("data has not been precomputed, returning without results")
    # return(Pagwas)
  }
  # fit model
  Pathway_names <- names(Pathway_ld_gwas_data)
  message("Run regression for ", length(Pathway_names), " pathways")
  pb <- txtProgressBar(style = 3)

  Pathway_sclm_results <- lapply(Pathway_names, function(Pathway) {

    Pathway_block <- Pathway_ld_gwas_data[[Pathway]]
    noise_per_snp <- Pathway_block$snps$se**2

    if (!is.null(Pathway_block$x)) {
      if (nrow(Pathway_block$x) > 2) {

        na_elements <- is.na(Pathway_block$y) | apply(Pathway_block$x, 1, function(x) {
          any(is.na(x))
        }) | is.na(noise_per_snp)


        x <-as(Pathway_block$x, "matrix")
        rownames(x) <- Pathway_block$snps$rsid
        results <- scParameter_regression(Pagwas_x = x[!na_elements, ], Pagwas_y = Pathway_block$y[!na_elements], noise_per_snp = noise_per_snp[!na_elements], n.cores = n.cores)
        results[is.na(results)] <- 0
      }else{
          return(NULL)
           }
        setTxtProgressBar(pb, which(Pathway_names == Pathway) / length(Pathway_names))
        return(results)
    } else {
      return(NULL)
    }
  })
  close(pb)

  names(Pathway_sclm_results) <- Pathway_names
  Pathway_sclm_results <- Pathway_sclm_results[!sapply(Pathway_sclm_results, is.null)]
  Pathway_names <- names(Pathway_sclm_results)
  Pathway_sclm_results <- as.data.frame(Pathway_sclm_results)
  Pathway_sclm_results<-as(data.matrix(Pathway_sclm_results), "dgCMatrix")

  #pca_scCell_mat<-as(pca_scCell_mat[Pathway_names, ],"dgCMatrix")

  data_mat<-as(data_mat,"dgCMatrix")

  pathway_expr <- lapply(Pathway_names, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], rownames(data_mat))

    if (length(a) == 0) {
      return(rep(0, ncol(data_mat)))
    } else if (length(a) == 1) {
      return(data_mat[intersect(a, rownames(data_mat)), ])
    } else {

      b <- apply(data_mat[intersect(a, rownames(data_mat)), ], 2, mean)
      return(b)
    }
  })
  #rm(data_mat)
  pathway_expr <- as.data.frame(pathway_expr)
  colnames(pathway_expr) <- Pathway_names
  pathway_expr <- as(data.matrix(pathway_expr), "dgCMatrix")
  pa_exp_mat <- as(t(as(pca_scCell_mat[Pathway_names, ],"matrix")), "dgCMatrix") * pathway_expr

  scPagwas_mat <- Pathway_sclm_results * pa_exp_mat

  rm(pa_exp_mat)
  rm(pathway_expr)

  scs <- rowSums(scPagwas_mat)
  rm(scPagwas_mat)
  SOAR::Store(Pathway_sclm_results)

  scs <- sign(scs) * log10(abs(scs) + 0.0001)

  df <- data.frame(cellid = colnames(pca_scCell_mat), scPagwas_score = scs)
  rownames(df) <- df$cellid
  rm(scs)
  #rm(pca_scCell_mat)
  gc()
  if (remove_outlier) {
    Pagwas$scPagwas_score <- scPagwas_score_filter(scPagwas_score = df$scPagwas_score)
  }
  names(Pagwas$scPagwas_score)<-df$cellid
  Pagwas$gene_heritability_correlation<-scGet_gene_heritability_correlation(
    scPagwas_score=Pagwas$scPagwas_score,
    data_mat=raw_data_mat[,names(Pagwas$scPagwas_score)])

  rm(df)

  return(Pagwas)
}



#' scParameter_regression
#' @description Find parameter estimates for the data.
#' @param Pagwas_x x parameter for lm
#' @param Pagwas_y y parameter for lm
#' @param noise_per_snp noise
#' @param n.cores cores
#'
#' @return

scParameter_regression <- function(Pagwas_x, Pagwas_y, noise_per_snp, n.cores = 1) {

  x_df <- bigstatsr::as_FBM(Pagwas_x, type = "double")

  liear_m <- bigstatsr::big_univLinReg(
    X = x_df,
    y.train = Pagwas_y,
    covar.train = bigstatsr::covar_from_df(data.frame(offset(noise_per_snp))),
    ncores = n.cores
  )


  return(liear_m$estim)
}


#' scPagwas_score_filter
#' @description filter the cPagwas_score for outliers.
#' @param scPagwas_score (data.frame)
#'
#' @return

scPagwas_score_filter <- function(scPagwas_score) {
  #scPagwas_score <- scPagwas_score$scPagwas_score
  # remove the NAN!
  if(sum(is.nan(scPagwas_score))>0){
    scPagwas_score[is.nan(scPagwas_score)]<-0
  }
  # remove the inf values!
  if (Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == Inf)] <- max(scPagwas_score[-which(scPagwas_score == Inf)],na.rm=TRUE)
  }
  if (-Inf %in% scPagwas_score) {
    scPagwas_score[which(scPagwas_score == -Inf)] <- min(scPagwas_score[-which(scPagwas_score == -Inf)],na.rm=TRUE)
  }
  lower_bound <- quantile(scPagwas_score, 0.01,na.rm=TRUE)
  upper_bound <- quantile(scPagwas_score, 0.99,na.rm=TRUE)

  lower_ind <- which(scPagwas_score < lower_bound)
  upper_ind <- which(scPagwas_score > upper_bound)
  scPagwas_score[lower_ind] <- lower_bound
  scPagwas_score[upper_ind] <- upper_bound

  return(scPagwas_score)
}



#' scGet_gene_heritability_contributions
#'
#' @param scPagwas_score result of scPagwas
#' @param data_mat the data matrix from single cell data
#'
#' @return
#' @examples
#' Single_data<-readRDS("E:/RPakage/scPagwas/inst/extdata/scRNAexample.rds")
#' Pagwas$sparse_cor<-scGet_gene_heritability_contributions(
#' scPagwas_score=Pagwas$scPagwas_score,
#' data_mat=Seurat::GetAssayData(object = Single_data[,names(Pagwas$scPagwas_score)], slot = "data"))

scGet_gene_heritability_correlation <- function(scPagwas_score,data_mat){

  if(all(names(scPagwas_score)==colnames(data_mat))){
    scPagwas_score<-data.matrix(scPagwas_score)
    sparse_cor<-corSparse(X=t(as(data_mat,"matrix")), Y=scPagwas_score)
  }else{
    data_mat<-data_mat[,names(scPagwas_score)]
    scPagwas_score<-data.matrix(scPagwas_score)
    sparse_cor<-corSparse(X=t(as(data_mat,"matrix")), Y=scPagwas_score)
  }
  return(sparse_cor)
}
