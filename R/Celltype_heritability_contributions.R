#' Celltype_heritability_contributions
#'
#' @param Pagwas Pagwas format, deault is NULL.
#' @param iters number of bootstrap iterations to perform
#' @param part number of bootstrap iterations to perform,default is 0.5
#' @param ncores Parallel cores,default is 1. use detectCores() to check the cores in computer.
#'
#' @return
#' @export
#'
#' @examples
Celltype_heritability_contributions<-function(Pagwas,
                                              iters = 200,
                                              part = 0.5){

  #Pathway_ld_gwas_data <- link_pwpca_block(Pagwas)
  #message("done")
  Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data=Pagwas$Pathway_ld_gwas_data)
  #message("done")
  # add on heritability values
  #message("** Get Pathway heritability contributions!")
  #Pagwas$Pathway_block_heritability <-
  #  Get_Pathway_heritability_contributions(
  #    Pagwas$pca_cell_df,
  #    Pagwas$lm_results$parameters
   # )
  #message("done")
  # Bootstrap error and 95% confidence interval estimates

  if (iters > 0) {
    Pagwas <- Boot_evaluate(Pagwas,
                            bootstrap_iters = iters, part = part)
    Pagwas$bootstrap_results$annotation<-c("Intercept",colnames(Pagwas$pca_cell_df))
  }
  Pagwas$Pathway_ld_gwas_data<-NULL
  if (sum(is.na(Pagwas$lm_results$parameters)) > 1) {
    stop("There is NA in parameters,can not appropriate to continue the calculation,Please check whether the GWAS data is too small!")
  }

return(Pagwas)
}

#' link_pwpca_block
#' @description Link the pca score and expression for each pathway genes
#' for each block
#' Requires rownames that are identitcal to block labels loaded previously.
#' @param Pagwas Pagwas format, deault is NULL.
#'
#' @return
#' @export
#' @examples
#' library(scPagwas)
#' # Pagwas should have inhibit data
#' Pagwas <- link_pwpca_block(Pagwas)

link_pwpca_block <- function(pa_block,
                             pca_cell_df,
                             merge_scexpr,
                             snp_gene_df,
                             rawPathway_list) {

    pathway <- unique(pa_block$block_info$pathway)

    x <-pca_cell_df[pathway, ]
    if(length(pathway)==1){
    x <- matrix(x, nrow = 1)
    rownames(x)<-pathway
    }

    if (nrow(pa_block$snps) == 0) {
      pa_block$include_in_inference <- F
      pa_block$x <- NULL # to make sure we totally replace previous stuffs
      return(pa_block)
      # stop("remove duplicates from pa_block data")
    }
    proper_genes <- rownames(merge_scexpr)
    mg <- intersect(rawPathway_list[[pathway]], proper_genes)
    x2 <- merge_scexpr[mg, ]

    if (length(mg) > 1) {
      x2 <- apply(x2, 2, function(x) (x - min(x)) / (max(x) - min(x)))

    }

    if (pa_block$n_snps > 1) {

      x2 <- x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <-x[rep(1, pa_block$n_snps), ]
      rownames(x) <- pa_block$snps$rsid

      #snp_gene_df <- Pagwas$snp_gene_df
      rownames(snp_gene_df) <- snp_gene_df$rsid
      x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
      x3 <- x2 * x

    } else {

      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label
      pa_block$n_snps <- nrow(pa_block$snps)

      rownames(x) <- pa_block$snps$rsid
      #snp_gene_df <- Pagwas$snp_gene_df
      rownames(snp_gene_df) <- snp_gene_df$rsid

      x <- matrix(as.numeric(x) * as.numeric(snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x3 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
    }

    rm(x)
    rm(x2)

    pa_block$x <- t(pa_block$ld_matrix_squared) %*% x3
    rownames(pa_block$x) <- pa_block$snps$rsid
    colnames(pa_block$x) <- colnames(merge_scexpr)
    rm(x3)

    pa_block$include_in_inference <- T

  return(pa_block)
}


#' Pa_Pagwas_perform_regression
#' @description Run regression
#' @param pa_block Pagwas format, deault is NULL.
#'
#' @return
#'
#' @examples
#' library(scPagwas)
#' Pagwas <- Pagwas_perform_regression(Pagwas, iters = 200)
Pa_Pagwas_perform_regression <- function(pa_block) {
  #Sys.setenv(R_LOCAL_CACHE=scPagwasSession)
  if (is.null(pa_block)) {
    warning("data has not been precomputed, returning without results")
    #return(Pagwas)
  }
  #message("** Start inference")
  # fit model
  if (!is.null(pa_block$x)) {
    if (pa_block$n_snps > 2) {
  vectorized_Pagwas_data <- xy2vector_pa(pa_block)

  lm_results <- pa_parameter_regression(vectorized_Pagwas_data)

    }else{
      lm_results<-NULL
    }

  } else {
    lm_results<-NULL
  }
  return(lm_results)
}


#' Pagwas_perform_regression
#' @description Run regression
#' @param Pathway_ld_gwas_data Pagwas format, deault is NULL.
#'
#' @return
#'
#' @examples
#' library(scPagwas)
#' Pagwas <- Pagwas_perform_regression(Pagwas, iters = 200)
Pagwas_perform_regression <- function(Pathway_ld_gwas_data) {
  #Sys.setenv(R_LOCAL_CACHE=scPagwasSession)
  if (is.null(Pathway_ld_gwas_data)) {
    warning("data has not been precomputed, returning without results")
    #return(Pagwas)
  }
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
#' @return

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


#' pa_parameter_regression
#' @description Find parameter estimates for the data.
#' @param vectorized_Pagwas_data Pagwas data that has been vectorized
#'
#' @return

pa_parameter_regression <- function(vectorized_Pagwas_data) {
  lm_results <- list()

  m <- stats::lm(vectorized_Pagwas_data$y ~
                   offset(vectorized_Pagwas_data$noise_per_snp) +
                   vectorized_Pagwas_data$x)

  lm_results$parameters <- stats::coef(m)

  lm_results$parameters[is.na(lm_results$parameters)] <- 0
  #names(results)<-colnames(data_mat)

  annotation_names <- c("Intercept", colnames(vectorized_Pagwas_data$x))
  names(lm_results$parameters) <- annotation_names
  #lm_results$model <- m

  return(lm_results$parameters[-1])
}


#' Boot_evaluate
#' @description Bootstrap to evaluate for confidence intervals.
#' @param Pagwas Pagwas format, deault is NULL.
#' @param bootstrap_iters number of bootstrap iterations to run,default is 200
#' @param part number of bootstrap iterations to perform,default is 0.5
#'
#' @return

Boot_evaluate <- function(Pagwas,
                          bootstrap_iters = 200,
                          part = 0.5) {

  # run things in parallel if user specified
  message(paste0("* starting bootstrap iteration for ", bootstrap_iters, " times"))

  pb <- txtProgressBar(style = 3)
  Boot_resultlist <-
    lapply(1:bootstrap_iters, function(i) {

      boot_results <- Parameter_regression(
        xy2vector(Pagwas$Pathway_ld_gwas_data[
          sample(seq_len(length(Pagwas$Pathway_ld_gwas_data)), floor(length(Pagwas$Pathway_ld_gwas_data) * part))
        ])
      )


      boot_results <- para_names_adjust(Pagwas, lm_results = boot_results)
      setTxtProgressBar(pb, i / bootstrap_iters)

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


#' para_names_adjust
#' @description Adjustment the names of parameters and pca_cell_df, for the unidentified signals;
#' @param Pagwas Pagwas format
#' @param lm_results Pagwas data list for lm_results, it is the result of lm
#'
#' @return

para_names_adjust <- function(Pagwas, lm_results = Pagwas$lm_results) {
  #pca_cell_df <- Pagwas$pca_cell_df
  if (sum(names(lm_results$parameters) %in% colnames(Pagwas$pca_cell_df)) < ncol(Pagwas$pca_cell_df)) {
    # message("There is blank or '+' within cell names!")
    names(lm_results$parameters) <- stringr::str_replace_all(names(lm_results$parameters), " ", ".")
    names(lm_results$parameters) <- stringr::str_replace_all(names(lm_results$parameters), "\\+", ".")
    names(lm_results$parameters) <- stringr::str_replace_all(names(lm_results$parameters), "-", ".")
  }

  if (sum(names(lm_results$parameters) %in% colnames(Pagwas$pca_cell_df)) < ncol(Pagwas$pca_cell_df)) {
    stop("unidentified signal within cell names, please remove it!")
  }
  return(lm_results)
}




#' xy2vector
#' @description Take a list of Pagwas - Pathway_ld_gwas_data and vectorize it.
#' @param Pathway_ld_gwas_data the list of block information from Pagwas object
#'
#' @return
#'

xy2vector <- function(Pathway_ld_gwas_data = NULL) {
  # use only blocks flagged for inference inclusion
  Pathway_ld_gwas_data <- Pathway_ld_gwas_data[sapply(Pathway_ld_gwas_data, function(block) {
    block$include_in_inference
  })]

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


xy2vector_pa <- function(pa_block = NULL) {
  # use only blocks flagged for inference inclusion

  # unpack Pathway_ld_gwas_data
  y <- pa_block$y
  x <- pa_block$x

  rownames(x) <- pa_block$snps$rsid
  noise_per_snp <- pa_block$snps$se**2

  na_elements <- is.na(y) | apply(x, 1, function(x) {
    any(is.na(x))
  }) | is.na(noise_per_snp)
  return(list(
    y = y[!na_elements], x = x[!na_elements, ],
    noise_per_snp = noise_per_snp[!na_elements]
  ))
}


#' Get_Pathway_heritability_contributions
#' @description Caclulate predicted block values based on block information and model fit.
#' @param pca_cell_df pca score dataframe
#' @param parameters parameter fit
#'
#' @return
#'
Get_Pathway_heritability_contributions <- function(pca_cell_df, parameters) {
  if (any(is.na(parameters))) {
    warning("NA pameters found!")
    parameters[is.na(parameters)] <- 0
  }

  # only include parameter for which we have block data
  Pathway_block_info <- as.numeric(data.matrix(pca_cell_df) %*% parameters[colnames(pca_cell_df)])
  names(Pathway_block_info) <- rownames(pca_cell_df)


  return(Pathway_block_info)
}



#' Get_bootresults_df
#' @description Helper function to make a summary table of results from bootstrap data.
#' @param value_collection collection of bootstrapped value estimates
#' @param annotations vector of annotation names
#' @param model_estimates estimates for bias parameter estimates
#'
#' @return

Get_bootresults_df <- function(value_collection, annotations, model_estimates) {

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



#' Pagwas_perform_regularized_inference
#'
#' @description Run inference with added regularization.
#' If p-values are desired use the other inference function. This for
#' prediction purposes.
#' @param Pagwas Pagwas data
#' @param n_folds folds for regularized inference
#'
#' @return

Pagwas_perform_regularized_regression <- function(Pagwas, n_folds = 10) {
  message("performing cross validation")
  Pagwas$cv_regularized_lm_results <- cv_regularized_parameter_estimator(
    xy2vector(),
    n_folds = n_folds
  )

  # add on block values
  Pagwas$regularized_Pathway_heritability_contributions <- Get_Pathway_heritability_contributions(
    Pagwas$pca_cell_df, Pagwas$cv_regularized_lm_results$parameters
  )

  return(Pagwas)
}

#' cv_regularized_parameter_estimator
#' @description Perform regularization inference.
#' Use CV to find appropriate values of lambda for either feature selection
#' @param vectorized_Pagwas_data Pagwas data used for inference
#' @param n_folds number of folds for cross validation
#' @param ... other arguments to pass to cv.glmnet
#'
#' @return

cv_regularized_parameter_estimator <- function(vectorized_Pagwas_data,
                                               n_folds = 10,
                                               ...) {
  lm_results <- list()
  m <- glmnet::cv.glmnet(
    x = vectorized_Pagwas_data$x,
    y = vectorized_Pagwas_data$y,
    offset = vectorized_Pagwas_data$noise_per_snp,
    foldid = cut(seq_len(length(vectorized_Pagwas_data$y)), breaks = n_folds, labels = F),
    family = "gaussian", ... = ...
  )

  # can choose coefficients with either: lambda.min, or lambda.1se
  lm_results$parameters <- stats::coef(m, s = "lambda.min") %>% as.numeric()
  annotation_names <- c("intercept", colnames(vectorized_Pagwas_data$x))
  names(lm_results$parameters) <- annotation_names
  lm_results$model <- m
  return(lm_results)
}




