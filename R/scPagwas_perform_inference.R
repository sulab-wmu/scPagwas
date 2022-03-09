#' scPagwas_perform_regression
#' @description Functions for inferring relevant annotations using the polyTest model.
#' @param Pagwas Pagwas data list, default is "NULL"
#' @param n.cores (integr)Parallel cores,default is 1. use detectCores() to check the cores in computer.
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' scPagwas_perform_regression(Pagwas, n.cores = n.cores)
scPagwas_perform_regression <- function(Pagwas, n.cores = n.cores) {
  if (is.null(Pagwas$Pathway_ld_gwas_data)) {
    warning("data has not been precomputed, returning without results")
    return(Pagwas)
  }
  message("Start inference")
  # fit model
  vectorized_Pagwas_data <- xy2vector(Pagwas$Pathway_ld_gwas_data)
  Pagwas$sclm_results <- scParameter_regression(vectorized_Pagwas_data)

  Pagwas$sclm_results$model$cellnames <- rownames(Pagwas$sclm_results$model$cellnames)
  Pagwas$sclm_results <- as.data.table(Pagwas$sclm_results$model)

  Pagwas$scPathway_heritability_contributions <-
    scGet_Pathway_heritability_contributions(
      Pagwas$pca_scCell_mat,
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
  pb <- txtProgressBar(style = 3)
  scBoot_evaluate <- papply(1:bootstrap_iters, function(i) {
    # message(paste('starting bootstrap iteration',  i, sep = ': '))

    part_vector <- xy2vector(Pagwas$Pathway_ld_gwas_data[
      sample(1:length(Pagwas$Pathway_ld_gwas_data), floor(length(Pagwas$Pathway_ld_gwas_data) * 0.25))
    ])
    boot_results <- scParameter_regression(part_vector)
    # boot_results <- scpara_names_adjust(Pagwas,lm_results=boot_results)
    # return the important bits
    return(boot_results$parameters)
  }, n.cores = n.cores)

  df <- as.data.frame(sapply(scBoot_evaluate, function(boot) {
    boot
  }))

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
  Pathway_block_info <- as.numeric(as.matrix(pca_scCell_mat) %*% parameters[colnames(pca_scCell_mat)])
  names(Pathway_block_info) <- rownames(pca_scCell_mat)
  return(Pathway_block_info)
}
