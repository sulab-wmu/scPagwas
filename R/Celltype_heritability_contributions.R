#' Pagwas_perform_regression
#' @description Run regression for each pathway block and gwas beta
#' parameter
#' @param Pathway_ld_gwas_data list of pathway block.
#'
#' @return regression results
#' @export
#' @author Chunyu Deng
#' @aliases Pagwas_perform_regression
#' @keywords Pagwas_perform_regression, Run regression function for
#' vectorized pathway blocks.
Pagwas_perform_regression <- function(Pathway_ld_gwas_data) {
  message("** Start inference")
  # fit model
  vectorized_Pagwas_data <- xy2vector(Pathway_ld_gwas_data)

  lm_results <- Parameter_regression(vectorized_Pagwas_data)

  return(lm_results)
}

#' Parameter_regression
#' @description Find parameter estimates for the data.
#' @param vectorized_Pagwas_data Pagwas data that has been vectorized
#'
#' @return regression for vectorized block data.

Parameter_regression <- function(vectorized_Pagwas_data) {
  lm_results <- list()

  m <- stats::lm(vectorized_Pagwas_data$y ~
  offset(vectorized_Pagwas_data$noise_per_snp) +
    vectorized_Pagwas_data$x)

  lm_results$parameters <- stats::coef(m)

  annotation_names <- c("Intercept", colnames(vectorized_Pagwas_data$x))
  names(lm_results$parameters) <- annotation_names
  lm_results$model <- m

  return(lm_results)
}


#' Boot_evaluate
#' @description Bootstrap to evaluate for confidence intervals.
#' @param Pagwas Pagwas format, deault is NULL.
#' @param bootstrap_iters number of bootstrap iterations to run,
#' default is 200
#' @param part number of bootstrap iterations to perform,default is 0.5
#' @export
#' @return Data frame for bootstrap results
#' @author Chunyu Deng
#' @aliases Boot_evaluate
#' @keywords Boot_evaluate, Integrate the bootstrap results.

Boot_evaluate <- function(Pagwas,
                          bootstrap_iters = 200,
                          part = 0.5) {

  # run things in parallel if user specified
  message(paste0(
    "* starting bootstrap iteration for ",
    bootstrap_iters, " times"
  ))

  pb <- utils::txtProgressBar(style = 3)
  Boot_resultlist <-
    lapply(1:bootstrap_iters, function(i) {
      boot_results <- Parameter_regression(
        xy2vector(Pagwas$Pathway_ld_gwas_data[
          sample(
            seq_len(length(Pagwas$Pathway_ld_gwas_data)),
            floor(length(Pagwas$Pathway_ld_gwas_data) * part)
          )
        ])
      )

      names(boot_results$parameters) <- c(
        "Intercept",
        colnames(Pagwas$pca_cell_df)
      )

      utils::setTxtProgressBar(pb, i / bootstrap_iters)

      # return the important bits
      return(
        list(
          boot_parameters = boot_results$parameters,
          block_heritability = Get_Pathway_heritability_contributions(
            Pagwas$pca_cell_df, boot_results$parameters
          )
        )
      )
    })
  close(pb)

  Pagwas$bootstrap_results <- Get_bootresults_df(
    sapply(Boot_resultlist, function(boot) {
      boot$boot_parameters
    }),
    names(Pagwas$lm_results$parameters),
    Pagwas$lm_results$parameters
  )
  return(Pagwas)
}


#' xy2vector
#' @description Take a list of Pagwas - Pathway_ld_gwas_data and vectorize
#' it.
#' @param Pathway_ld_gwas_data the list of block information from Pagwas
#' object
#'
#' @return vector for x and y values.
#'

xy2vector <- function(Pathway_ld_gwas_data = NULL) {
  # use only blocks flagged for inference inclusion
  Pathway_ld_gwas_data <- Pathway_ld_gwas_data[sapply(
    Pathway_ld_gwas_data,
    function(block) {
      block$include_in_inference
    }
  )]

  # unpack Pathway_ld_gwas_data
  y <- do.call("c", lapply(Pathway_ld_gwas_data, function(block) {
    block$y
  }))
  x <- do.call("rbind", lapply(Pathway_ld_gwas_data, function(block) {
    block$x
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


#' Get_Pathway_heritability_contributions
#' @description Caclulate predicted block values based on block
#' information and model fit.
#' @param pca_cell_df pca score dataframe
#' @param parameters parameter fit
#'
#' @return vector for pathway heritability contributions
#'
Get_Pathway_heritability_contributions <- function(pca_cell_df,
                                                   parameters) {
  if (any(is.na(parameters))) {
    warning("NA pameters found!")
    parameters[is.na(parameters)] <- 0
  }

  # only include parameter for which we have block data
  Pathway_block_info <- as.numeric(data.matrix(pca_cell_df) %*%
    parameters[colnames(pca_cell_df)])
  names(Pathway_block_info) <- rownames(pca_cell_df)


  return(Pathway_block_info)
}



#' Get_bootresults_df
#' @description Helper function to make a summary table of results from
#' bootstrap data.
#' @param value_collection collection of bootstrapped value estimates
#' @param annotations vector of annotation names
#' @param model_estimates estimates for bias parameter estimates
#'
#' @return data frame for bootstrap results.

Get_bootresults_df <- function(value_collection, annotations,
                               model_estimates) {
  bootstrap_estimate <- bootstrap_error <- bt_value <- NULL
  # in the case we're calculating single parameter estimates
  if (is.null(dim(value_collection))) {
    value_collection <- matrix(value_collection, nrow = 1)
  }

  parameter_estimates <- data.frame(
    annotation = annotations,
    bootstrap_estimate = rowMeans(value_collection),
    bootstrap_error = apply(value_collection, 1, stats::sd)
  ) %>%
    dplyr::mutate(
      bt_value = bootstrap_estimate / bootstrap_error,
      bp_value = stats::pnorm(bt_value, lower.tail = F),
      bias_corrected_estimate = 2 * model_estimates - bootstrap_estimate
    )

  parameter_estimates$CI_lo <-
    apply(value_collection, 1, function(x) {
      stats::quantile(x, 0.025, na.rm = T)
    })
  parameter_estimates$CI_hi <-
    apply(value_collection, 1, function(x) {
      stats::quantile(x, 0.975, na.rm = T)
    })

  return(parameter_estimates)
}
