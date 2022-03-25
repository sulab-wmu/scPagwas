
#' Singlecell_heritability_contributions
#'
#' @param Pagwas Pagwas format, deault is NULL.
#' @param iters number of bootstrap iterations to perform
#' @param part number of bootstrap iterations to perform,default is 0.5
#' @param n.cores Parallel cores,default is 1. use detectCores() to check the cores in computer.
#' @param SCregression whether to run regression funcitons, which takes a lot of time and memmory.
#' @return
#' @export
#'
#' @examples

Singlecell_heritability_contributions<-function(Pagwas,
                                              iters = 200,
                                              n.cores=1,
                                              part = 0.5,
                                              SCregression=FALSE){
  message("** link single cell and gwas snp!")

  Pathway_ld_gwas_data <- link_scCell_pwpca_block(Pagwas, n.cores = n.cores)
  message("done")

  message("** Get scPagwas_score!")
  Pagwas$scPagwas_score <- scPagwas_perform_score(Pagwas=Pagwas,
                                                  Pathway_ld_gwas_data=Pathway_ld_gwas_data,
                                                  n.cores=n.cores)
  message("done")
  message("** Get gene heritability contributions!")
  Pagwas<-scGet_gene_heritability_correlation(Pagwas=Pagwas)

  message("done")
  # add on heritability values
  if(SCregression){
    Pagwas<-scPagwas_perform_regression(Pagwas,
                                Pathway_ld_gwas_data=Pathway_ld_gwas_data,
                                n.cores = 1)
  }


  return(Pagwas)
}

#' link_scCell_pwpca_block
#' @description link the single cell pca score and expression for each
#' pathway genes for each block(single cell)
#' construct the final matrix for regression in blocks
#' @param Pagwas Pagwas format, deault is NULL.
#' @param n.cores Parallel cores,default is 1. use detectCores() to
#' check the cores in computer.
#' @param scPagwasSession "scPagwasSession"
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)


link_scCell_pwpca_block <- function(Pagwas, n.cores = 1) {

  #Sys.setenv(R_LOCAL_CACHE=scPagwasSession)

  if (is.null(Pagwas$Pathway_ld_gwas_data)) {
    message("* no loaded Pathway_ld_gwas data")
    return(Pagwas)
  }

  if (dim(Pagwas$data_mat)[2] == dim(Pagwas$pca_scCell_mat)[2]) {
    Pagwas$dim_data_mat$col <- Pagwas$dim_pca_scCell_mat$col

    } else {
    #if (length(a) == 0) {
      stop("* please check the colnames of Singlecell data.There may have specific symbol")

  }

  message("* Merging pathway score and expression information about blocks in ", length(Pagwas$Pathway_ld_gwas_data), " pathways")
  pb <- txtProgressBar(style = 3)

  Pathway_ld_gwas_data <- papply(Pagwas$Pathway_ld_gwas_data, function(pa_block) {

    pathway <- unique(pa_block$block_info$pathway)

    x <- Pagwas$pca_scCell_mat[which(Pagwas$dim_pca_scCell_mat$row==pathway), ]

   # if(length(pathway)==1){
      x <- matrix(x, nrow = 1)
      rownames(x)<-pathway
    #}

    if (nrow(pa_block$snps) == 0) {
      pa_block$include_in_inference <- F
      pa_block$x <- NULL # to make sure we totally replace previous stuffs
      return(pa_block)

    }

    mg <- intersect(Pagwas$rawPathway_list[[pathway]],Pagwas$dim_data_mat$row)
    if (length(mg) < 2) {
      return(NULL)
    }

    if (length(mg) > 1) {
      colnorm_sub <- function(X, ind) {
        if (sum(X[, ind]) == 0) {
          return(NA)
        } else {
          return(X[, ind] / sum(X[, ind]))
        }
      }
      x2 <- bigstatsr::big_apply(Pagwas$data_mat[match(mg,Pagwas$dim_data_mat$row), ], a.FUN = colnorm_sub, a.combine = "cbind")
      rownames(x2) <- mg
    }

    if (pa_block$n_snps > 1) {
      x2 <-x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- x[rep(1, pa_block$n_snps), ]
      rownames(x) <- pa_block$snps$rsid
      rownames(Pagwas$snp_gene_df) <- Pagwas$snp_gene_df$rsid

      x <- x * Pagwas$snp_gene_df[pa_block$snps$rsid, "slope"]

      x3 <- x2 * x
      #x3 <- bigstatsr::as_FBM(x3)
    } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
      rownames(x) <- pa_block$snps$rsid

      rownames(Pagwas$snp_gene_df) <- Pagwas$snp_gene_df$rsid
      x <- matrix(as.numeric(x) * as.numeric(Pagwas$snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x3 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)

    }
    #rm(x)
    #rm(x2)
    pa_block$x<- bigstatsr::as_FBM( pa_block$ld_matrix_squared %*% x3)
    pa_block$include_in_inference <- T

   # gc()
    setTxtProgressBar(pb, which(names(Pagwas$Pathway_ld_gwas_data) == pathway) / length(Pagwas$Pathway_ld_gwas_data))

    return(pa_block)
  }, n.cores = n.cores)
  close(pb)

  Pathway_ld_gwas_data <- Pathway_ld_gwas_data[!sapply(Pathway_ld_gwas_data, is.null)]
  gc()
  return(Pathway_ld_gwas_data)
}

#' scPagwas_perform_score
#' @description Get the scPagwas score for each cells
#' @param Pagwas pagwas
#' @param Pathway_ld_gwas_data Pagwas data list
#' @param n.cores cores
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' @return
#'
#' @examples
#' library(scPagwas)
#' scPagwas_perform_score(Pagwas)
scPagwas_perform_score <- function(Pagwas,
                                   Pathway_ld_gwas_data,
                                   n.cores = 1,
                                   remove_outlier=TRUE) {

  if (is.null(Pathway_ld_gwas_data)) {
    stop("data has not been precomputed, returning without results")

  }
  # fit model
  Pathway_names <- names(Pathway_ld_gwas_data)
  message("*Run regression for ", length(Pathway_names), " pathways")
  pb <- txtProgressBar(style = 3)

  Pathway_sclm_results <- lapply(Pathway_names, function(Pathway) {

    Pathway_block <- Pathway_ld_gwas_data[[Pathway]]
    noise_per_snp <- Pathway_block$snps$se**2

    if (!is.null(Pathway_block$x)) {
      if (nrow(Pathway_block$x) > 2) {

        na_elements <- is.na(Pathway_block$y) | apply(Pathway_block$x[], 1, function(x) {
          any(is.na(x))
        }) | is.na(noise_per_snp)
        Pathway_block$x<- bigstatsr::as_FBM(Pathway_block$x[!na_elements, ])
        results <- scParameter_regression(Pagwas_x = Pathway_block$x, Pagwas_y = Pathway_block$y[!na_elements], noise_per_snp = noise_per_snp[!na_elements], n.cores = n.cores)

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
  Pathway_sclm_results<-bigstatsr::as_FBM(data.matrix(Pathway_sclm_results))

  pathway_expr <- lapply(Pathway_names, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], Pagwas$dim_data_mat$row)

    if (length(a) == 0) {
      return(rep(0, ncol(Pagwas$data_mat)))
    } else if (length(a) == 1) {
      return(Pagwas$data_mat[match(a,Pagwas$dim_data_mat$row), ])
    } else {

      b <- apply(Pagwas$data_mat[match(a,Pagwas$dim_data_mat$row), ], 2, mean)
      return(b)
    }
  })
  #rm(data_mat)
  pathway_expr <- as.data.frame(pathway_expr)
  colnames(pathway_expr) <- Pathway_names
  pathway_expr <- bigstatsr::as_FBM(data.matrix(pathway_expr))
  pa_exp_mat <- bigstatsr::as_FBM(t(Pagwas$pca_scCell_mat[match(Pathway_names,Pagwas$dim_pca_scCell_mat$row), ]) * pathway_expr[])

  scPagwas_mat <- Pathway_sclm_results[] * pa_exp_mat[]

  rm(pa_exp_mat)
  rm(pathway_expr)

  scs <- rowSums(scPagwas_mat)
  rm(scPagwas_mat)
  #SOAR::Store(Pathway_sclm_results)

  scs <- sign(scs) * log10(abs(scs) + 0.0001)

  df <- data.frame(cellid = Pagwas$dim_pca_scCell_mat$col, scPagwas_score = scs)
  rownames(df) <- df$cellid
  rm(scs)
  #rm(pca_scCell_mat)
  gc()
  if (remove_outlier) {
    scPagwas_score <- scPagwas_score_filter(scPagwas_score = df$scPagwas_score)
  }
  names(scPagwas_score)<-df$cellid

  rm(df)

  return(scPagwas_score)
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

  #x_df <- bigstatsr::as_FBM(Pagwas_x, type = "double")

  liear_m <- bigstatsr::big_univLinReg(
    X = Pagwas_x,
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
#' @param Pagwas result of scPagwas
#'
#' @return
#' @examples
#' Single_data<-readRDS("E:/RPakage/scPagwas/inst/extdata/scRNAexample.rds")
#' Pagwas$sparse_cor<-scGet_gene_heritability_contributions(
#' scPagwas_score=Pagwas$scPagwas_score,
#' data_mat=Seurat::GetAssayData(object = Single_data[,names(Pagwas$scPagwas_score)], slot = "data"))

scGet_gene_heritability_correlation <- function(Pagwas){

  if(all(names(Pagwas$scPagwas_score)==Pagwas$dim_raw_data_mat$col)){
    scPagwas_score<-data.matrix(Pagwas$scPagwas_score)
    sparse_cor<-corSparse(X=t(Pagwas$raw_data_mat[]), Y=scPagwas_score)
  }else{
    data_mat<-Pagwas$raw_data_mat[,match(names(scPagwas_score),Pagwas$dim_raw_data_mat$col)]
    scPagwas_score<-data.matrix(scPagwas_score)
    sparse_cor<-corSparse(X=t(Pagwas$raw_data_mat[]), Y=scPagwas_score)
  }
  Pagwas$gene_heritability_correlation<-sparse_cor
  rownames(Pagwas$gene_heritability_correlation)<-Pagwas$dim_raw_data_mat$row
  return(Pagwas)
}




#' scPagwas_perform_regression
#' @description Functions for inferring relevant annotations using the polyTest model.
#' @param Pagwas Pagwas data list, default is "NULL"
#' @param n.cores (integr)Parallel cores,default is 1. use detectCores() to check the cores in computer.
#' @param Pathway_ld_gwas_data Pathway_ld_gwas_data
#' @return
#'
#' @examples
#' library(scPagwas)
#' scPagwas_perform_regression(Pagwas, n.cores = n.cores)
scPagwas_perform_regression <- function(Pagwas,
                                        Pathway_ld_gwas_data,
                                        n.cores = 1) {
 #Sys.setenv(R_LOCAL_CACHE=scPagwasSession)
  if (is.null(Pathway_ld_gwas_data)) {
    stop("data has not been precomputed, returning without results")
    #return(NULL)
  }
  message("Start inference")
  # fit model
  vectorized_Pagwas_data <- SCxy2vector(Pathway_ld_gwas_data)


  sclm_results <- scParameter_regression(Pagwas_x=data.matrix(vectorized_Pagwas_data[[2]]),
                                         Pagwas_y=vectorized_Pagwas_data[[1]],
                                         noise_per_snp=vectorized_Pagwas_data[[3]],
                                         n.cores = 1)

  sclm_results[is.na(sclm_results)] <- 0
  names(sclm_results)<-colnames(vectorized_Pagwas_data[[2]])
  rm(vectorized_Pagwas_data)
  Pagwas$sclm_results <- sclm_results
  rm(sclm_results)

  Pagwas$scPathway_heritability_contributions <-
    scGet_Pathway_heritability_contributions(
      Pagwas$pca_scCell_mat[],
      Pagwas$sclm_results
    )

  return(Pagwas)
}


#' scBoot_evaluate
#' @description Bootstrap to evaluate for confidence intervals.
#' @param Pagwas Pagwas object
#' @param bootstrap_iters number of bootstrap iterations to run
#' @param n.cores cores
#'
#' @return

scBoot_evaluate <- function(Pagwas, bootstrap_iters = bootstrap_iters, n.cores = n.cores) {
  # run things in parallel if user specified

  scBoot_evaluate <- papply(1:bootstrap_iters, function(i) {

    part_vector <- xy2vector(SC_Pathway_ld_gwas_data[
      sample(seq_len(length(SC_Pathway_ld_gwas_data)),
             floor(length(SC_Pathway_ld_gwas_data) * 0.25))
    ])
    part_vector$x<-as(part_vector$x,"matrix")
    boot_results <- scParameter_regression(part_vector)

    return(boot_results$parameters)
  }, n.cores = n.cores)

  #SOAR::Store(Pathway_ld_gwas_data)

  df <- as.data.frame(sapply(scBoot_evaluate, function(boot) {
    boot
  }))
  rm(scBoot_evaluate)
  Pagwas$scbootstrap_results <- Get_bootresults_df(
    value_collection = df,
    annotations = names(Pagwas$sclm_results$parameters),
    model_estimates = Pagwas$sclm_results$parameters
  )
  return(Pagwas)
}



#' scGet_Pathway_heritability_contributions
#' @description predicting the importance of a pathway based on its pcascore
#' @param pca_scCell_mat pca score matrix
#' @param parameters parameter fit
#'
#' @return

scGet_Pathway_heritability_contributions <- function(pca_scCell_mat, parameters) {
  if (any(is.na(parameters))) {
    warning("NA pameters found!")
    parameters[is.na(parameters)] <- 0
  }

  # only include parameter for which we have block data

  Pathway_block_info <- as.numeric(pca_scCell_mat %*% parameters[intersect(Pagwas$dim_pca_scCell_mat$col,names(parameters))])
  names(Pathway_block_info) <- Pagwas$dim_pca_scCell_mat$row

  return(Pathway_block_info)
}

#' SCxy2vector
#' @description Take a list of Pagwas - Pathway_ld_gwas_data and vectorize it.
#' @param Pathway_ld_gwas_data the list of block information from Pagwas object
#'
#' @return
#'
SCxy2vector <- function(Pathway_ld_gwas_data = NULL) {
  # use only blocks flagged for inference inclusion
  Pathway_ld_gwas_data <- Pathway_ld_gwas_data[sapply(Pathway_ld_gwas_data, function(block) {
    block$include_in_inference
  })]

  # unpack Pathway_ld_gwas_data
  y <- do.call("c", lapply(Pathway_ld_gwas_data, function(block) {
    block$y
  }))
  x <- do.call("rbind", lapply(Pathway_ld_gwas_data, function(block) {
    block$x[]
  }))

  rownames(x) <- do.call("c", lapply(Pathway_ld_gwas_data, function(block) {
    block$snps$rsid
  }))
  noise_per_snp <- do.call("c", lapply(Pathway_ld_gwas_data, function(block) {
    block$snps$se**2
  }))

  # exclude na elements
  na_elements <- is.na(y) | apply(x, 1, function(x) {
    any(is.na(x))
  }) | is.na(noise_per_snp)
  return(list(
    y = y[!na_elements], x = x[!na_elements, ],
    noise_per_snp = noise_per_snp[!na_elements]
  ))
}



