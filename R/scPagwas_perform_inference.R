#' scPagwas_perform_regression
#' @description Functions for inferring relevant annotations using the polyTest model.
#' @param Pagwas Pagwas data list, default is "NULL"
#' @param n.cores (integr)Parallel cores,default is 1. use detectCores() to check the cores in computer.
#' @param scPagwasSession "scPagwasSession"
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' scPagwas_perform_regression(Pagwas, n.cores = n.cores)
scPagwas_perform_regression <- function(Pagwas, n.cores = 1,scPagwasSession="scPagwasSession") {
  Sys.setenv(R_LOCAL_CACHE=scPagwasSession)
  if (is.null(SC_Pathway_ld_gwas_data)) {
    warning("data has not been precomputed, returning without results")
    return(Pagwas)
  }
  message("Start inference")
  # fit model
  vectorized_Pagwas_data <- xy2vector(SC_Pathway_ld_gwas_data)


  sclm_results <- scParameter_regression(Pagwas_x=as(vectorized_Pagwas_data[[2]],"matrix"),
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
      pca_scCell_mat,
      Pagwas$sclm_results
    )
  if("Pathway_list" %in% names(Pagwas)){
    Pagwas$Pathway_list<-NULL
  }
  if("rawPathway_list" %in% names(Pagwas)){
    Pagwas$rawPathway_list<-NULL
  }
  if("VariableFeatures" %in% names(Pagwas)){
    Pagwas$VariableFeatures<-NULL
  }
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

  Pathway_block_info <- as.numeric(pca_scCell_mat %*% parameters[colnames(pca_scCell_mat)])
  names(Pathway_block_info) <- rownames(pca_scCell_mat)

  return(Pathway_block_info)
}
