
#' Singlecell_heritability_contributions
#'
#' @param Pagwas Pagwas format, deault is NULL.
#' @param part number of bootstrap iterations to perform,default is 0.5
#' @param remove_outlier
#' @param n.cores Parallel cores,default is 1. use detectCores() to check the cores in computer.
#' @param SimpleResult TRUE
#' @return
#' @export
#'
#' @examples

Singlecell_heritability_contributions<-function(Pagwas,
                                              n.cores=1,
                                              part = 0.5,
                                              remove_outlier=TRUE,
                                              SimpleResult=F
                                              ){
  message("** link single cell and gwas snp and Get scPagwas_score!")
  Pagwas  <- link_scCell_pwpca_block(Pagwas,
                                     n.cores = n.cores,
                                     remove_outlier=TRUE)
  message("done")

  message("** Get gene heritability contributions!")
  Pagwas <- scGet_gene_heritability_correlation(Pagwas=Pagwas)
  message("done")

  message("** Get Pathway heritability contributions!")
  #Pagwas <- scGet_Pathway_heritability_correlation(Pagwas=Pagwas)

  if(SimpleResult){
    Pagwas[c("VariableFeatures","merge_scexpr","data_mat","rawPathway_list","Pathway_list","pca_scCell_mat","snp_gene_df","Pathway_ld_gwas_data")]<-NULL
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
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)

link_scCell_pwpca_block <- function(Pagwas,
                                    n.cores = 1,
                                    remove_outlier=T) {
  #options(bigmemory.allow.dimnames=TRUE)
  #Sys.setenv(R_LOCAL_CACHE=scPagwasSession)

  if (is.null(Pagwas$Pathway_ld_gwas_data)) {
    message("* no loaded Pathway_ld_gwas data")
    return(Pagwas)
  }

  if (dim(Pagwas$data_mat)[2] != dim(Pagwas$pca_scCell_mat)[2]) {
  stop("* please check the colnames of Singlecell data.There may have specific symbol")
  }

  message("* Merging pathway score and expression information about blocks in ", length(Pagwas$Pathway_ld_gwas_data), " pathways")
  pb <- txtProgressBar(style = 3)

  paths<-names(Pagwas$Pathway_ld_gwas_data)
  Pathway_sclm_results <- papply(Pagwas$Pathway_ld_gwas_data, function(pa_block) {

    pathway <- unique(pa_block$block_info$pathway)
    #message(pathway)
    x <- matrix(Pagwas$pca_scCell_mat[pathway, ],nrow = 1)
    rownames(x)<-pathway
    #message(1)
    #}

    if (pa_block$n_snps == 0) {
      pa_block$include_in_inference <- F
      pa_block$x <- NULL # to make sure we totally replace previous stuffs
      return(pa_block)
    }

    mg <- intersect(Pagwas$rawPathway_list[[pathway]],rownames(Pagwas$data_mat))
    if (length(mg) == 1) {
      x2<-matrix(Pagwas$data_mat[mg, ],nrow=1)
      x2<- x2/(x2+0.0001)
      rownames(x2)<-mg
    }else{
      x2 <-  biganalytics::apply(Pagwas$data_mat[mg, ],2,function(ge){
          if (sum(ge) == 0) {
            return(rep(0,length(ge)))
          }else{
            return(ge / sum(ge))
          }
      })
     rownames(x2)<-mg
    }
    #message(2)
    #if(ncol(x2)>10000){
      x2 <- as(x2,"dgCMatrix")
      #message(3)
    #}
    if (pa_block$n_snps > 1) {
      x2 <-x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- x[rep(1, pa_block$n_snps), ]
      rownames(x) <- pa_block$snps$rsid
      rownames(Pagwas$snp_gene_df) <- Pagwas$snp_gene_df$rsid
      #message(4)
      x <- x * Pagwas$snp_gene_df[pa_block$snps$rsid, "slope"]
      #if(ncol(x2)>10000){
      #x <- as(x,"dgCMatrix")
      #}
      x2 <- x2 * x
      #message(5)
      #x3 <- bigstatsr::as_FBM(x3)
    } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
      rownames(x) <- pa_block$snps$rsid

      rownames(Pagwas$snp_gene_df) <- Pagwas$snp_gene_df$rsid
      x <- matrix(as.numeric(x) * as.numeric(Pagwas$snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x2 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
      #if(ncol(x2)>10000){
      x2 <- as(x2,"dgCMatrix")
      #}
    }
    #rm(x)
    #rm(x2)
    #message(6)
    #if(ncol(x2)>bignumber){
    pa_block$x<- as(pa_block$ld_matrix_squared %*% x2,"dgCMatrix")
    #}else{
    #pa_block$x<- pa_block$ld_matrix_squared %*% x2

   # }
    #message(7)
    pa_block$include_in_inference <- T
    #Pathway_block <- Pathway_ld_gwas_data[[Pathway]]
    noise_per_snp <- pa_block$snps$se**2

    if (!is.null(pa_block$x)) {
      if (pa_block$n_snps > 2) {

        na_elements <- is.na(pa_block$y) | apply(pa_block$x, 1, function(x) {
          any(is.na(x))
        }) | is.na(noise_per_snp)

        results <- scParameter_regression(Pagwas_x = pa_block$x[!na_elements,],
                                          Pagwas_y = pa_block$y[!na_elements],
                                          noise_per_snp = noise_per_snp[!na_elements],
                                          n.cores = n.cores
        )


        results[is.na(results)] <- 0
        #Pathway_cell_regression<-results
      }else{
        results<-NULL
      }

    } else {
      results<-NULL
    }

    setTxtProgressBar(pb, which(paths == pathway) / length(paths))

    return(results)
  }, n.cores = n.cores)
  close(pb)

###############################################

  names(Pathway_sclm_results) <- paths
  Pathway_sclm_results <- Pathway_sclm_results[!sapply(Pathway_sclm_results, is.null)]
  paths <- names(Pathway_sclm_results)
  Pathway_sclm_results <- data.matrix(as.data.frame(Pathway_sclm_results))
  ncells<- nrow(Pathway_sclm_results)
  #if(ncells>bignumber){
  #  Pathway_sclm_results <- bigmemory::as.big.matrix(Pathway_sclm_results,shared = FALSE)
  #}else{
  Pathway_sclm_results<-as(Pathway_sclm_results,"dgCMatrix")
  Pagwas$Pathway_sclm_results<-Pathway_sclm_results
  colnames(Pagwas$Pathway_sclm_results)<-paths
 # }

  message("* Get pathways mean expression in single cell")

  pathway_expr <- lapply(paths, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], rownames(Pagwas$data_mat))

    if (length(a) == 0) {
      #return(rep(0, ncol(Pagwas$data_mat)))
      return(NULL)
    } else if (length(a) == 1) {
      #return(Pagwas$data_mat[a, ])
      return(NULL)
    } else {
      b <- biganalytics::apply(Pagwas$data_mat[a, ], 2, mean)
      return(b)
    }
  })
  #rm(data_mat)
  paths<-paths[!unlist(lapply(pathway_expr,is.null))]
  #pa_exp_mat<-pa_exp_mat[,paths]
  Pathway_sclm_results<-Pathway_sclm_results[,paths]
  pathway_expr<-pathway_expr[!unlist(lapply(pathway_expr,is.null))]
  pathway_expr <- data.matrix(as.data.frame(pathway_expr))

  #if(ncells>bignumber){
  #  pathway_expr <- bigmemory::as.big.matrix(pathway_expr,shared = FALSE)
  # colnames(pathway_expr) <- paths
  # pa_exp_mat <-bigmemory::as.big.matrix(t(Pagwas$pca_scCell_mat[paths, ]) * pathway_expr[],shared = FALSE)
  # rm(pathway_expr)
  # message("* Get scPgwas score for each single cell")
  # scs <- rowSums(Pathway_sclm_results[] * pa_exp_mat[])
  #}else{
  pathway_expr<-as(pathway_expr,"dgCMatrix")
  colnames(pathway_expr) <- paths
  pa_exp_mat <-t(Pagwas$pca_scCell_mat[paths, ]) * pathway_expr
  rm(pathway_expr)

  pa_exp_mat<-as(pa_exp_mat,"dgCMatrix")
  message("* Get scPgwas score for each single cell")
  scs <- rowSums(Pathway_sclm_results * pa_exp_mat)
  #}

  rm(pa_exp_mat)

  df <- data.frame(cellid = colnames(Pagwas$pca_scCell_mat), scPagwas_score = sign(scs) * log10(abs(scs) + 0.0001))
  rownames(df) <- df$cellid
  #rm(pca_scCell_mat)
  gc()
  if (remove_outlier) {
    scPagwas_score <- scPagwas_score_filter(scPagwas_score = df$scPagwas_score)
  }
  names(scPagwas_score)<-df$cellid
  Pagwas$scPagwas_score<-scPagwas_score

  return(Pagwas)
}

#' scPagwas_perform_score
#' @description Get the scPagwas score for each cells
#' @param Pagwas pagwas
#' @param Pathway_ld_gwas_data Pagwas data list
#' @param n.cores cores
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' @param bignumber the threshold of number of cells to change matrix to bigmemory data format, if you have enough
#' storage, you can set a large number to save time.
#' @return
#'
#' @examples
#' library(scPagwas)
#' scPagwas_perform_score(Pagwas)
scPagwas_perform_score <- function(Pagwas,
                                   Pathway_ld_gwas_data,
                                   n.cores = 1,
                                   remove_outlier=TRUE,
                                   bignumber=10000
                                   ) {
  options(bigmemory.allow.dimnames=TRUE)
  if (is.null(Pathway_ld_gwas_data)) {
    stop("data has not been precomputed, returning without results")
  }
  # fit model
  Pathway_names <- names(Pathway_ld_gwas_data)
  message("*Run regression for ", length(Pathway_names), " pathways")
  pb <- txtProgressBar(style = 3)

  Pathway_sclm_results <- lapply(Pathway_ld_gwas_data, function(Pathway_block) {

    #Pathway_block <- Pathway_ld_gwas_data[[Pathway]]
    noise_per_snp <- Pathway_block$snps$se**2

    if (!is.null(Pathway_block$x)) {
      if (Pathway_block$n_snps > 2) {

        na_elements <- is.na(Pathway_block$y) | apply(Pathway_block$x, 1, function(x) {
          any(is.na(x))
        }) | is.na(noise_per_snp)

        results <- scParameter_regression(Pagwas_x = Pathway_block$x[!na_elements,],
                                          Pagwas_y = Pathway_block$y[!na_elements],
                                          noise_per_snp = noise_per_snp[!na_elements],
                                          n.cores = n.cores
                                          )


        results[is.na(results)] <- 0
      }else{
          return(NULL)
      }
      Pathway<-Pathway_block$block_info$pathway[1]
        setTxtProgressBar(pb, which(Pathway_names == Pathway) / length(Pathway_names))
        return(results)
    } else {
      return(NULL)
    }
  })
  close(pb)
  rm(Pathway_ld_gwas_data)

  names(Pathway_sclm_results) <- Pathway_names
  Pathway_sclm_results <- Pathway_sclm_results[!sapply(Pathway_sclm_results, is.null)]
  Pathway_names <- names(Pathway_sclm_results)
  Pathway_sclm_results <- data.matrix(as.data.frame(Pathway_sclm_results))
  ncells<- nrow(Pathway_sclm_results)
  if(ncells>bignumber){
  Pathway_sclm_results <- bigmemory::as.big.matrix(Pathway_sclm_results,shared = FALSE)
  }else{
  Pathway_sclm_results<-as(Pathway_sclm_results,"dgCMatrix")
  }

  message("* Get pathways mean expression in single cell")

  pathway_expr <- lapply(Pathway_names, function(pa) {
    a <- intersect(Pagwas$Pathway_list[[pa]], rownames(Pagwas$data_mat))

    if (length(a) == 0) {
      #return(rep(0, ncol(Pagwas$data_mat)))
      return(NULL)
    } else if (length(a) == 1) {
      #return(Pagwas$data_mat[a, ])
      return(NULL)
    } else {
      b <- biganalytics::apply(Pagwas$data_mat[a, ], 2, mean)
      return(b)
    }
  })
  #rm(data_mat)
  pathway_expr <- data.matrix(as.data.frame(pathway_expr))

  #if(ncells>bignumber){
  #  pathway_expr <- bigmemory::as.big.matrix(pathway_expr,shared = FALSE)
   # colnames(pathway_expr) <- Pathway_names
   # pa_exp_mat <-bigmemory::as.big.matrix(t(Pagwas$pca_scCell_mat[Pathway_names, ]) * pathway_expr[],shared = FALSE)
   # rm(pathway_expr)
   # message("* Get scPgwas score for each single cell")
   # scs <- rowSums(Pathway_sclm_results[] * pa_exp_mat[])
  #}else{
    pathway_expr<-as(pathway_expr,"dgCMatrix")
    colnames(pathway_expr) <- Pathway_names
    pa_exp_mat <-t(Pagwas$pca_scCell_mat[Pathway_names, ]) * pathway_expr
    rm(pathway_expr)

    pa_exp_mat<-as(pa_exp_mat,"dgCMatrix")
    message("* Get scPgwas score for each single cell")
    scs <- rowSums(Pathway_sclm_results * pa_exp_mat)
  #}

  rm(pa_exp_mat)

  df <- data.frame(cellid = colnames(Pagwas$pca_scCell_mat), scPagwas_score = sign(scs) * log10(abs(scs) + 0.0001))
  rownames(df) <- df$cellid
  #rm(pca_scCell_mat)
  gc()
  if (remove_outlier) {
    scPagwas_score <- scPagwas_score_filter(scPagwas_score = df$scPagwas_score)
  }
  names(scPagwas_score)<-df$cellid

  return(scPagwas_score)
}



#' scParameter_regression
#' @description Find parameter estimates for the data.
#' @param Pagwas_x x parameter for lm
#' @param Pagwas_y y parameter for lm
#' @param noise_per_snp noise
#' @param n.cores 1
#'
#' @return

scParameter_regression <- function(Pagwas_x,
                                   Pagwas_y,
                                   noise_per_snp,
                                   n.cores = 1) {

  #x_df <-Pagwas_x

  if(bigmemory::is.big.matrix(Pagwas_x)){
  #  liear_m<- biglm.big.matrix(Pagwas_y ~ offset(noise_per_snp) + Pagwas_x)
  liear_m <- bigstatsr::big_univLinReg(
    X = as_FBM(Pagwas_x[]),
    y.train = Pagwas_y,
    covar.train = bigstatsr::covar_from_df(data.frame(offset(noise_per_snp))),
    ncores = n.cores
  )
  parameters<-  liear_m$estim
  }else if(is(Pagwas_x, 'dgCMatrix')){
    Pagwas_x<-as_matrix(Pagwas_x)
    liear_m <- bigstatsr::big_univLinReg(
      X = as_FBM(Pagwas_x),
      y.train = Pagwas_y,
      covar.train = bigstatsr::covar_from_df(data.frame(offset(noise_per_snp))),
      ncores = n.cores
    )
    parameters<-  liear_m$estim
  }else{
    #Pagwas_x<-as_matrix(Pagwas_x)
    liear_m <- bigstatsr::big_univLinReg(
      X = as_FBM(Pagwas_x),
      y.train = Pagwas_y,
      covar.train = bigstatsr::covar_from_df(data.frame(offset(noise_per_snp))),
      ncores = n.cores
    )
    parameters<-  liear_m$estim

  }

  return(parameters)
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

  if(all(names(Pagwas$scPagwas_score)==colnames(Pagwas$data_mat))){
    scPagwas_score<-data.matrix(Pagwas$scPagwas_score)
    sparse_cor<-corSparse(X=t(as_matrix(Pagwas$data_mat)),
      Y=scPagwas_score)
  }else{
    data_mat<-Pagwas$data_mat[,names(scPagwas_score)]
    scPagwas_score<-data.matrix(scPagwas_score)
    sparse_cor<-corSparse(X=t(as_matrix(Pagwas$data_mat)), Y=scPagwas_score)
  }
  Pagwas$gene_heritability_correlation<-sparse_cor
  rownames(Pagwas$gene_heritability_correlation)<-rownames(Pagwas$data_mat)
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


  Pagwas$sclm_results <- scParameter_regression(Pagwas_x=data.matrix(vectorized_Pagwas_data[[2]]),
                                         Pagwas_y=vectorized_Pagwas_data[[1]],
                                         noise_per_snp=vectorized_Pagwas_data[[3]],
                                         n.cores = n.cores)

  Pagwas$sclm_results[is.na(Pagwas$sclm_results)] <- 0
  names(Pagwas$sclm_results)<-colnames(vectorized_Pagwas_data[[2]])
  rm(vectorized_Pagwas_data)


  if (iters > 0) {
    Pagwas <- scBoot_evaluate(Pagwas,
                            Pathway_ld_gwas_data=Pathway_ld_gwas_data,
                            bootstrap_iters = 200,
                            n.cores = 1,
                            part = 0.5)
    #Pagwas$bootstrap_results$annotation<-c("Intercept",colnames(Pagwas$pca_cell_df))
  }
  return(Pagwas)
}


#' scBoot_evaluate
#' @description Bootstrap to evaluate for confidence intervals.
#' @param Pagwas Pagwas object
#' @param bootstrap_iters number of bootstrap iterations to run
#' @param n.cores cores
#' @param part number of bootstrap iterations to perform,default is 0.5
#'
#' @return

scBoot_evaluate <- function(Pagwas,
                            Pathway_ld_gwas_data,
                            bootstrap_iters,
                            n.cores = 1,
                            part = 0.5) {
  # run things in parallel if user specified
  message(paste0("* starting bootstrap iteration for ", bootstrap_iters, " times"))

  pb <- txtProgressBar(style = 3)
  scBoot_evaluate <- papply(1:bootstrap_iters, function(i) {

    part_vector <- SCxy2vector(Pathway_ld_gwas_data[
      sample(seq_len(length(Pathway_ld_gwas_data)),
             floor(length(Pathway_ld_gwas_data) * part))
    ])

    results <- scParameter_regression(Pagwas_x = part_vector$x,
                                      Pagwas_y = part_vector$y,
                                      noise_per_snp = part_vector$noise_per_snp,
                                      n.cores=n.cores
                                      )

    setTxtProgressBar(pb, i / bootstrap_iters)
    return(results)
  }, n.cores = n.cores)

  #SOAR::Store(Pathway_ld_gwas_data)
  df <- as.data.frame(sapply(scBoot_evaluate, function(boot) {
    boot
  }))
  rm(scBoot_evaluate)
  Pagwas$scbootstrap_results <- Get_bootresults_df(
    value_collection = df,
    annotations = names(Pagwas$sclm_results),
    model_estimates = Pagwas$sclm_results
  )
  return(Pagwas)
}



#' scGet_Pathway_heritability_contributions
#' @description predicting the importance of a pathway based on its pcascore
#' @param pca_scCell_mat pca score matrix
#' @param parameters parameter fit
#'
#' @return

scGet_Pathway_heritability_contributions <- function(pca_scCell_mat, Pathway_sclm_results) {
  if (any(is.na(parameters))) {
    warning("NA pameters found!")
    parameters[is.na(parameters)] <- 0
  }

  # only include parameter for which we have block data

  Pathway_block_info <- as.numeric(pca_scCell_mat %*% parameters[intersect(colnames(pca_scCell_mat),names(parameters))])
  names(Pathway_block_info) <- rownames(pca_scCell_mat)

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
