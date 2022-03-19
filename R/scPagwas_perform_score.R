
#' scPagwas_perform_score
#' @description Get the scPagwas score for each cells
#' @param Pagwas Pagwas data list
#' @param n.cores cores
#' @param split_n (integr)default NULL, When the cell number is too big, there may have memory errors, set split_n=10000 or other number can split the cell data.
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' scPagwas_perform_score(Pagwas)
scPagwas_perform_score <- function(Pagwas,
                                   n.cores = 1,
                                   split_n = NULL,
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

      x <- Pathway_block$x
      rownames(x) <- Pathway_block$snps$rsid

    if(!is.null(split_n)){
      if (ncol(Pathway_block$x) > split_n) {
          # message("There are too much cells, regression progress will split the data!")

          cellsN <- colnames(Pathway_block$x)

          n <- ceiling(ncol(Pathway_block$x) / split_n)
          chunk <- function(x, n) split(x, factor(sort(rank(x) %% n)))
          chunk_list <- chunk(x = cellsN, n)

          results <- unlist(lapply(chunk_list, function(index) {
            sclm_1 <- scParameter_regression(Pagwas_x = x[!na_elements, index], Pagwas_y = Pathway_block$y[!na_elements], noise_per_snp = noise_per_snp[!na_elements], n.cores = n.cores)
            return(sclm_1)
          }))
          rm(x)
          rm(noise_per_snp)

      }else{
          results <- scParameter_regression(Pagwas_x = x[!na_elements, ], Pagwas_y = Pathway_block$y[!na_elements], noise_per_snp = noise_per_snp[!na_elements], n.cores = n.cores)

          }
      }else {
          results <- scParameter_regression(Pagwas_x = x[!na_elements, ], Pagwas_y = Pathway_block$y[!na_elements], noise_per_snp = noise_per_snp[!na_elements], n.cores = n.cores)
           }
       rm(x)
       rm(noise_per_snp)


        results[is.na(results)] <- 0
        setTxtProgressBar(pb, which(Pathway_names == Pathway) / length(Pathway_names))
        return(results)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })
  close(pb)


  names(Pathway_sclm_results) <- Pathway_names
  Pathway_sclm_results <- Pathway_sclm_results[!sapply(Pathway_sclm_results, is.null)]
  Pathway_names <- names(Pathway_sclm_results)
  Pathway_sclm_results <- as.data.frame(lapply(Pathway_sclm_results, function(x) {
    x[is.na(x)] <- 0
    return(x)
  }))

  pca_scCell_mat <- pca_scCell_mat[Pathway_names, ]

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
  pathway_expr <- as.data.frame(pathway_expr)
  colnames(pathway_expr) <- Pathway_names


  pathway_expr <- data.matrix(pathway_expr[colnames(pca_scCell_mat), rownames(pca_scCell_mat)])
  pa_exp_mat <- t(pca_scCell_mat) * pathway_expr
  scPagwas_mat <- data.matrix(Pathway_sclm_results) * pa_exp_mat

  rm(pa_exp_mat)
  rm(pathway_expr)
  scs <- rowSums(scPagwas_mat)
  rm(scPagwas_mat)
  SOAR::Store(Pathway_sclm_results)

  scs <- sign(scs) * log10(abs(scs) + 0.0001)

  df <- data.frame(cellid = colnames(pca_scCell_mat), scPagwas_score = scs)
  rownames(df) <- df$cellid
  rm(scs)

  if (remove_outlier) {
    Pagwas$scPagwas_score <- scPagwas_score_filter(scPagwas_score = df$scPagwas_score)
  }
  names(Pagwas$scPagwas_score)<-df$cellid
  rm(df)
  gc()
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

#' scGet_gene_heritability_contributions
#'
#' @param scPagwas_score result of scPagwas
#' @param data_mat the data matrix from single cell data
#'
#' @return
#' @example
#' Single_data<-readRDS("E:/RPakage/scPagwas/inst/extdata/scRNAexample.rds")
Pagwas$sparse_cor<-scGet_gene_heritability_contributions(
scPagwas_score=Pagwas$scPagwas_score,
data_mat=Seurat::GetAssayData(object = Single_data[,names(Pagwas$scPagwas_score)], slot = "data"))

scGet_gene_heritability_contributions <- function(scPagwas_score,data_mat){

 if(all(names(scPagwas_score)==colnames(data_mat))){
   scPagwas_score<-data.matrix(scPagwas_score)
   data_mat<-t(data_mat)
   sparse_cor<-corSparse(X=data_mat, Y=scPagwas_score)
 }else{
   data_mat<-data_mat[,names(scPagwas_score)]
   scPagwas_score<-data.matrix(scPagwas_score)
   sparse_cor<-corSparse(X=data_mat, Y=scPagwas_score)
 }
return(sparse_cor)
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
    scPagwas_score[is.nan(scPagwas_score)]<-min(scPagwas_score)
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
