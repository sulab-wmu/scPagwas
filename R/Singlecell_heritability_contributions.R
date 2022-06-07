
#' scPagwas_perform_score
#' @description Get the scPagwas score for each cells
#' @param Pagwas pagwas
#' @param remove_outlier Whether to remove the outlier for scPagwas score.
#' storage, you can set a large number to save time.
#' @return
#' @export
#' @examples
#' library(scPagwas)
#' scPagwas_perform_score(Pagwas)
scPagwas_perform_score <- function(Pagwas,
                                   remove_outlier = TRUE) {
  options(bigmemory.allow.dimnames = TRUE)
  ###########sclm_score
  a<-matrix(Pagwas$sclm_results,ncol = 1)
  mat<-Pagwas$pca_scCell_mat[colnames(weighted_singlecell_mat),]  * t(data.matrix(weighted_singlecell_mat))
  Pagwas$sclm_allsnpscore<-rowSums(apply(mat,1,function(x){
   x*a
  }),na.rm = T)

  gc()
  if (remove_outlier) {
    Pagwas$sclm_allsnpscore<-scPagwas_score_filter(scPagwas_score = Pagwas$sclm_allsnpscore)
  }
  names(Pagwas$sclm_allsnpscore) <- colnames(Pagwas$pca_scCell_mat)

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
  lower_bound <- quantile(scPagwas_score, 0.05, na.rm = TRUE)
  upper_bound <- quantile(scPagwas_score, 0.95, na.rm = TRUE)

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

    sparse_lmcor <- corSparse(
      X = t(as_matrix(Pagwas$data_mat)),
      Y = data.matrix(Pagwas$sclm_allsnpscore)
    )
  rownames(sparse_lmcor) <- rownames(Pagwas$data_mat)
  colnames(sparse_lmcor) <- "allsnp_gene_heritability_correlation"
  sparse_lmcor[is.nan(sparse_lmcor)] <- 0
  Pagwas$allsnp_gene_heritability_correlation <- sparse_lmcor

  return(Pagwas)
}



#' scPagwas_perform_regression
#' @description Functions for inferring relevant annotations using the polyTest model.
#' @param Pagwas Pagwas data list, default is "NULL"
#' @param Pathway_ld_gwas_data Pathway_ld_gwas_data
#' @return
#'
#' @export
#' @examples
#' library(scPagwas)
#' scPagwas_perform_regression(Pagwas, ncores = ncores)
scPagwas_perform_regression <- function(Pagwas,
                                        Pathway_ld_gwas_data,
                                        ncores = 1) {
  # Sys.setenv(R_LOCAL_CACHE=scPagwasSession)
  if (is.null(Pathway_ld_gwas_data)) {
    stop("data has not been precomputed, returning without results")
    # return(NULL)
  }
  message("Start inference")
  # fit model
  vectorized_Pagwas_data <- SCxy2vector(Pathway_ld_gwas_data)


  Pagwas$sclm_results <- scParameter_regression(
    Pagwas_x = data.matrix(vectorized_Pagwas_data[[2]]),
    Pagwas_y = vectorized_Pagwas_data[[1]],
    noise_per_snp = vectorized_Pagwas_data[[3]],
    ncores = ncores
  )

  Pagwas$sclm_results[is.na(Pagwas$sclm_results)] <- 0
  names(Pagwas$sclm_results) <- colnames(vectorized_Pagwas_data[[2]])
  rm(vectorized_Pagwas_data)


  return(Pagwas)
}


#' scBoot_evaluate
#' @description Bootstrap to evaluate for confidence intervals.
#' @param Pagwas Pagwas object
#' @param bootstrap_iters number of bootstrap iterations to run
#' @param part number of bootstrap iterations to perform,default is 0.5
#'
#'
#' @return

scBoot_evaluate <- function(Pagwas,
                            Pathway_ld_gwas_data,
                            bootstrap_iters,
                            ncores = 1,
                            part = 0.5) {
  # run things in parallel if user specified
  message(paste0("* starting bootstrap iteration for ", bootstrap_iters, " times"))

  pb <- txtProgressBar(style = 3)
  scBoot_evaluate <- lapply(1:bootstrap_iters, function(i) {
    part_vector <- SCxy2vector(Pathway_ld_gwas_data[
      sample(
        seq_len(length(Pathway_ld_gwas_data)),
        floor(length(Pathway_ld_gwas_data) * part)
      )
    ])

    results <- scParameter_regression(
      Pagwas_x = part_vector$x,
      Pagwas_y = part_vector$y,
      noise_per_snp = part_vector$noise_per_snp,
      ncores = ncores
    )

    setTxtProgressBar(pb, i / bootstrap_iters)
    return(results)
  })

   #rm(scPathway_ld_gwas_data)
  df <- as.data.frame(sapply(scBoot_evaluate, function(boot) {
    boot
  }))
  rm(scBoot_evaluate)
  Pagwas$scbootstrap_results <- scGet_bootresults_df(
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
#' @export
#' @return

scGet_Pathway_heritability_contributions <- function(pca_scCell_mat, Pathway_sclm_results) {
  if (any(is.na(parameters))) {
    warning("NA pameters found!")
    parameters[is.na(parameters)] <- 0
  }

  # only include parameter for which we have block data

  Pathway_block_info <- as.numeric(pca_scCell_mat %*% parameters[intersect(colnames(pca_scCell_mat), names(parameters))])
  names(Pathway_block_info) <- rownames(pca_scCell_mat)

  # T_scrna<-AddModuleScore(T_scrna,list(Anti_inflammatory,Pro_inflammatory),name=c("Anti_inflammatory","Pro_inflammatory"))

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

#' scGet_bootresults_df
#' @description Helper function to make a summary table of results from bootstrap data.
#' @param value_collection collection of bootstrapped value estimates
#' @param annotations vector of annotation names
#' @param model_estimates estimates for bias parameter estimates
#'
#' @return

scGet_bootresults_df <- function(value_collection, annotations, model_estimates) {

  # in the case we're calculating single parameter estimates
  if (is.null(dim(value_collection))) {
    value_collection <- matrix(c(value_collection,model_estimates), nrow = 1)
  }

  parameter_estimates <- data.frame(
    annotation = annotations,
    model_estimate = model_estimates,
    bootstrap_error = apply(value_collection, 1, stats::sd)
  ) %>%
    dplyr::mutate(
      bt_value = model_estimate / bootstrap_error,
      bp_value = stats::pnorm(bt_value, lower.tail = F),
      bias_corrected_estimate = 2 * rowMeans(value_collection) - model_estimate
    )

  parameter_estimates$CI_lo <-
    apply(value_collection, 1, function(x) {
      stats::quantile(x, 0.1, na.rm = T)
    })
  parameter_estimates$CI_hi <-
    apply(value_collection, 1, function(x) {
      stats::quantile(x, 0.9, na.rm = T)
    })

  return(parameter_estimates)
}

