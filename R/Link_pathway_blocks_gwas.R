
#' Link_pathway_blocks_gwas
#'
#' @description Link pathway blocks region and gwas summary snps
#' Takes block information, potentially independent LD pathway_blocks
#' or gene blocks,
#' @param Pagwas Pagwas format, deault is NULL.
#' @param chrom_ld LD data for 22 chromosome.
#' @param split_n number of times to compute the singlecell result
#' @param ncores Parallel cores,default is 1. use detectCores()
#' to check the cores in computer.
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(chrom_ld)
#' Pagwas <- Link_pathway_blocks_gwas(Pagwas = Pagwas, chrom_ld = chrom_ld)
Link_pathway_blocks_gwas <- function(Pagwas,
                                     chrom_ld = NULL,
                                     split_n = 1,
                                     singlecell = T,
                                     celltype = T,
                                     ncores = 1) {
  options(bigmemory.allow.dimnames = TRUE)
  Pachrom_block_list <- lapply(Pagwas$pathway_blocks, function(pa_blocks) split(pa_blocks, f = as.vector(pa_blocks$chrom)))

  names(Pachrom_block_list) <- names(Pagwas$pathway_blocks)
  Pagwas$pathway_blocks <- NULL

  chrom_gwas_list <- lapply(split(Pagwas$gwas_data, f = Pagwas$gwas_data$chrom), function(gwas) {
    gwas <- data.table::data.table(gwas)
    data.table::setkey(gwas, pos)
    return(gwas)
  })
  Pagwas$gwas_data <- NULL

  if (!singlecell & celltype) split_n <- 1
#  if (split_n == 1) {
    Pagwas <- Pathway_block_func(
      Pagwas = Pagwas,
      Pachrom_block_list = Pachrom_block_list,
      chrom_gwas_list = chrom_gwas_list,
      singlecell = singlecell,
      celltype = celltype,
      ncores = ncores
    )
#  }
  # if (split_n > 1) {
  #   pa_blocksnplist <- Pachrom_func(
  #     Pagwas = Pagwas,
  #     Pachrom_block_list = Pachrom_block_list,
  #     chrom_gwas_list = chrom_gwas_list
  #   )
  #
  #   if (celltype) {
  #     Pagwas <- celltypes_Pathway_block_func(Pagwas = Pagwas, pa_blocksnplist = pa_blocksnplist)
  #   }
  #   # if(singlecell){}
  #   a <- ncol(Pagwas$pca_scCell_mat)
  #   for (ai in 1:10) {
  #     # print(ai)
  #     if (a %% split_n == 0) break
  #     split_n <- split_n + 1
  #   }
  #
  #   la <- gl(split_n, a / split_n, length = a)
  #
  #   sclm_list <- list()
  #   for (i in 1:split_n) {
  #     message("** start to run the ", i, " split!")
  #     # sclm_list[[i]]<- matrix(nrow=sum(la==i), ncol=nrow(Pagwas$pca_scCell_mat))
  #     sclm_list[[i]] <- scPathway_block_splitfunc(
  #       pa_blocksnplist = pa_blocksnplist,
  #       subpca_scCell_mat = Pagwas$pca_scCell_mat[, which(la == i)],
  #       subdata_mat = Pagwas$data_mat[, which(la == i)],
  #       rawPathway_list = Pagwas$rawPathway_list,
  #       snp_gene_df = Pagwas$snp_gene_df,
  #       ncores = ncores
  #     )
  #   }
  #   # Reduce(function(dtf1, dtf2) cbind(dtf1, dtf2),sclm_list)
  #   Pagwas$Pathway_sclm_results <- as(data.matrix(bigreadr::rbind_df(sclm_list)), "dgCMatrix")
  #   rownames(Pagwas$Pathway_sclm_results) <- colnames(Pagwas$pca_scCell_mat)
  #   rm(pa_blocksnplist)
  # }

  rm(chrom_gwas_list)
  rm(chrom_ld)

  Pagwas$merge_scexpr <- NULL
  gc()
  return(Pagwas)
}

#' make_ld_matrix
#' @description construct the final matrix for regression in blocks
#' @param all_snps The snps that were queried
#' @param ld_data A returned LD matrix with SNP
#'
#' @return
#'
make_ld_matrix <- function(all_snps = snp_data$rsid, ld_data = sub_ld) {
  mat_dim <- length(all_snps)
  ld_matrix <- diag(mat_dim)

  if (mat_dim == 1) {
    return(data.matrix(1))
  } # no need to check psd
  if (mat_dim >= 2) {
    rownames(ld_matrix) <- all_snps
    colnames(ld_matrix) <- all_snps
    for (n in seq_len(nrow(ld_data))) {
      n_x <- as.vector(unlist(ld_data[n, ]))
      ld_matrix[n_x[1], n_x[2]] <- as.numeric(n_x[3])
      ld_matrix[n_x[2], n_x[1]] <- as.numeric(n_x[3])
    }
  }
  return(ld_matrix)
}


#' Title
#'
#' @param Pagwas
#' @param Pachrom_block_list
#' @param ncores
#' @param chrom_gwas_list
#' @param singlecell
#' @param celltype
#'
#' @return
#' @export
#'
#' @examples
Pathway_block_func <- function(Pagwas = NULL,
                               Pachrom_block_list,
                               chrom_gwas_list,
                               ncores = 1,
                               singlecell = T,
                               celltype = T) {
  Pathway_ctlm_results <- list()
  Pathway_ld_gwas_data <- list()
  scPathway_ld_gwas_data <- list()
  weighted_singlecell_mat <- list()
  weighted_celltypes_mat <- list()

  message(paste0("* Start to link gwas and pathway block annotations for ",
    length(Pachrom_block_list), " pathways!"
  ))
  options(bigmemory.allow.dimnames = TRUE)
  pb <- txtProgressBar(style = 3)

  for (pathway in names(Pachrom_block_list)) {
    #print(pathway)
    Pa_chrom_block <- Pachrom_block_list[[pathway]]

    Pa_chrom_data <- lapply(names(Pa_chrom_block), function(chrom) {
      chrom_block <- Pa_chrom_block[[chrom]]
      ld_data <- chrom_ld[[chrom]]

      data.table::setkey(ld_data, SNP_A) # just in case
      if (!(chrom %in% names(chrom_gwas_list))) {
        warning(paste(chrom, " for gwas is missing, could be a problem!", sep = ""))
        return(NULL)
      }
      if (is.null(chrom_gwas_list[[chrom]])) {
        warning(paste(chrom, " data missing, could be a problem", sep = ""))
        return(NULL)
      }

      rsids <- Pagwas$snp_gene_df[which(Pagwas$snp_gene_df$label %in%
        chrom_block$label), c("rsid", "label")]

      c2 <- chrom_gwas_list[[chrom]]
      rsids_gwas <- suppressMessages(dplyr::inner_join(rsids, c2))
      if (is.null(nrow(rsids_gwas))) {
        return(NULL)
      }
      beta_squared <- rsids_gwas$beta^2
      # again, pretty sure this is binarize data table lookup
      sub_ld <- ld_data[.(as.vector(rsids_gwas$rsid)), nomatch = 0L]

      return(list(rsids_gwas, beta_squared, sub_ld, chrom_block))
    })
    #print(1)
    # cbind_df代替 rbind_df(list_df)
    rm(Pa_chrom_block)
    pa_block <- list()
    pa_block$block_info <- bigreadr::rbind_df(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[4]]
    ))
    pa_block$snps <- bigreadr::rbind_df(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[1]]
    ))
    pa_block$y <- unlist(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[2]]
    ))
    #print(2)
    # message(paste0("snp is ",nrow(snp_data)))
    pa_block$ld_data <- bigreadr::rbind_df(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[3]]
    ))

    rsid_x <- intersect(pa_block$snps$rsid, unique(unlist(pa_block$ld_data[, 1:2])))

    pa_block$ld_data <- as.data.frame(pa_block$ld_data[pa_block$ld_data$SNP_A
      %in% rsid_x & pa_block$ld_data$SNP_B
      %in% rsid_x, ])
    rm(Pa_chrom_data)
    # message(paste0("ld is ",nrow(sub_ld)))
    if (nrow(pa_block$ld_data) == 0) {
      pa_block$ld_matrix_squared <- diag(1, nrow = nrow(pa_block$snps))
    } else {
      ld_matrix <- make_ld_matrix(
        all_snps = pa_block$snps$rsid,
        ld_data = pa_block$ld_data
      )
      pa_block$ld_matrix_squared <- ld_matrix * ld_matrix
    }

    pa_block$n_snps <- nrow(pa_block$snps)

    if (singlecell) {
      #Pathway_sclm_results[[pathway]] <- get_Pathway_sclm(
       # Pagwas=Pagwas,
       # pathway=pathway,
      #  path_block = pa_block,
      #  ncores = ncores
      #)

      scPathway_ld_gwas_data[[pathway]]<- link_pwpca_block(
        pa_block = pa_block,
        pca_cell_df = data.matrix(Pagwas$pca_scCell_mat),
        merge_scexpr = Pagwas$data_mat,
        snp_gene_df = Pagwas$snp_gene_df,
        rawPathway_list = Pagwas$rawPathway_list
      )
      weighted_singlecell_mat[[pathway]]<-colSums(as.data.frame(scPathway_ld_gwas_data[[pathway]]$x2))
      scPathway_ld_gwas_data[[pathway]]$x2<-NULL
    }

    if (celltype) {

      Pathway_ld_gwas_data[[pathway]] <-link_pwpca_block(
        pa_block = pa_block,
        pca_cell_df = data.matrix(Pagwas$pca_cell_df),
        merge_scexpr = Pagwas$merge_scexpr,
        snp_gene_df = Pagwas$snp_gene_df,
        rawPathway_list = Pagwas$rawPathway_list
      )
      if(pa_block$n_snps>=10){

         scv<-Pathway_ld_gwas_data[[pathway]]
         scv$noise_per_snp <- pa_block$snps$se**2

         # exclude na elements
         na_elements <- is.na(scv$y) | apply(scv$x, 1, function(x) {
           any(is.na(x))
         }) | is.na(scv$noise_per_snp)
         scv<-list(
           y = scv$y[!na_elements], x = scv$x[!na_elements, ],
           noise_per_snp = scv$noise_per_snp[!na_elements]
         )
         scl<-scParameter_regression(Pagwas_x=scv$x,
                                            Pagwas_y=scv$y,
                                            noise_per_snp=scv$noise_per_snp,
                                            ncores = ncores)
         #scl<-Parameter_regression(scv)
         #scl<-scl$model$coefficients[-1]
         scl[is.na(scl)]<-NA
         Pathway_ctlm_results[[pathway]]<-scl

      }else{
         Pathway_ctlm_results[[pathway]]<-NULL
      }
      weighted_celltypes_mat[[pathway]]<- colSums(as.data.frame(Pathway_ld_gwas_data[[pathway]]$x2))
      Pathway_ld_gwas_data[[pathway]]$x2<-NULL
    }
    setTxtProgressBar(pb, which(names(Pachrom_block_list) == pathway) / length(names(Pachrom_block_list)))
  }

  close(pb)

  if (singlecell) {
    names(scPathway_ld_gwas_data) <- names(Pachrom_block_list)
    weighted_singlecell_mat <- as.data.frame(weighted_singlecell_mat)

    SOAR::Store(scPathway_ld_gwas_data)
    SOAR::Store(weighted_singlecell_mat)

  }

  if (celltype) {
    names(Pathway_ld_gwas_data) <- names(Pachrom_block_list)
    SOAR::Store(Pathway_ld_gwas_data)

    if(length(Pathway_ctlm_results)>1){

    Pathway_ctlm_results <- Pathway_ctlm_results[!sapply(Pathway_ctlm_results, is.null)]
    Pathway_ctlm_results <- data.matrix(as.data.frame(Pathway_ctlm_results))
    rownames(Pathway_ctlm_results) <- colnames(Pagwas$pca_cell_df)
    Pagwas$Pathway_ctlm_results <- as(Pathway_ctlm_results, "dgCMatrix")
    }
    weighted_celltypes_mat <- as.data.frame(weighted_celltypes_mat)
    SOAR::Store(weighted_celltypes_mat)
  }
  return(Pagwas)
}



#' Pathway_block_splitfunc
#'
#' @param Pagwas
#' @param Pachrom_block_list
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples
scPathway_block_splitfunc <- function(pa_blocksnplist,
                                      subpca_scCell_mat,
                                      subdata_mat,
                                      rawPathway_list,
                                      snp_gene_df,
                                      ncores = 1) {
  Pathway_sclm_results <- list()
  #+++++++Pathway_lm_results<-list()
  message("* Start to link gwas and pathway block annotations for pathways!")
  options(bigmemory.allow.dimnames = TRUE)
  pb <- txtProgressBar(style = 3)

  for (pathway in names(pa_blocksnplist)) {
    pa_block <- pa_blocksnplist[[pathway]]
    Pathway_sclm_results[[pathway]] <- subget_Pathway_sclm(
      pathway=pathway,
      path_block = pa_block,
      subdata_mat = subdata_mat,
      rawPathway_list = rawPathway_list,
      snp_gene_df = snp_gene_df,
      ncores = ncores
    )
    setTxtProgressBar(pb, which(names(pa_blocksnplist) == pathway) / length(names(pa_blocksnplist)))
  }

  close(pb)
  Pathway_sclm_results <- Pathway_sclm_results[!sapply(Pathway_sclm_results, is.null)]
  Pathway_sclm_results <- as.data.frame(Pathway_sclm_results)
  rownames(Pathway_sclm_results) <- colnames(subpca_scCell_mat)

  return(Pathway_sclm_results)
}


#' celltypes_Pathway_block_func
#'
#' @param Pagwas
#' @param pa_blocksnplist
#'
#' @return
#' @export
#'
#' @examples
celltypes_Pathway_block_func <- function(Pagwas, pa_blocksnplist) {
  # Pathway_ld_gwas_data<-list()

  message("* Start to run celltypes ")
  options(bigmemory.allow.dimnames = TRUE)
  pb <- txtProgressBar(style = 3)

  Pathway_ld_gwas_data <- lapply(names(pa_blocksnplist), function(pathway) {

    # for (pathway in names(pa_blocksnplist)) {
    pa_block <- pa_blocksnplist[[pathway]]

    pa_block <- link_pwpca_block(
      pa_block = pa_block,
      pca_cell_df = data.matrix(Pagwas$pca_cell_df),
      merge_scexpr = Pagwas$merge_scexpr,
      snp_gene_df = Pagwas$snp_gene_df,
      rawPathway_list = Pagwas$rawPathway_list
    )
    # Pathway_lm_results[[pathway]]<- Pa_Pagwas_perform_regression(pa_block=pa_block)
    # Pathway_ld_gwas_data[[pathway]]<-pa_block
    setTxtProgressBar(pb, which(names(pa_blocksnplist) == pathway) / length(names(pa_blocksnplist)))
    return(pa_block)
  })
  close(pb)
  names(Pathway_ld_gwas_data) <- names(pa_blocksnplist)
  Pagwas$Pathway_ld_gwas_data <- Pathway_ld_gwas_data

  return(Pagwas)
}

#' Title
#'
#' @param Pagwas
#' @param Pachrom_block_list
#' @param chrom_gwas_list
#'
#' @return
#' @export
#'
#' @examples
Pachrom_func <- function(Pagwas, Pachrom_block_list, chrom_gwas_list) {
  message("* Start to link pathway and snps ")

  options(bigmemory.allow.dimnames = TRUE)

  pb <- txtProgressBar(style = 3)

  pa_blocksnplist <- lapply(names(Pachrom_block_list), function(pathway) {
    Pa_chrom_block <- Pachrom_block_list[[pathway]]

    Pa_chrom_data <- lapply(names(Pa_chrom_block), function(chrom) {
      chrom_block <- Pa_chrom_block[[chrom]]
      ld_data <- chrom_ld[[chrom]]

      data.table::setkey(ld_data, SNP_A) # just in case
      if (!(chrom %in% names(chrom_gwas_list))) {
        warning(paste(chrom, " for gwas is missing, could be a problem!", sep = ""))
        return(NULL)
      }
      if (is.null(chrom_gwas_list[[chrom]])) {
        warning(paste(chrom, " data missing, could be a problem", sep = ""))
        return(NULL)
      }

      rsids <- Pagwas$snp_gene_df[which(Pagwas$snp_gene_df$label %in%
        chrom_block$label), c("rsid", "label")]

      c2 <- chrom_gwas_list[[chrom]]
      rsids_gwas <- suppressMessages(dplyr::inner_join(rsids, c2))
      if (is.null(nrow(rsids_gwas))) {
        return(NULL)
      }
      beta_squared <- rsids_gwas$beta^2
      # again, pretty sure this is binarize data table lookup
      sub_ld <- ld_data[.(as.vector(rsids_gwas$rsid)), nomatch = 0L]
      setTxtProgressBar(pb, which(names(Pachrom_block_list) == pathway) / length(names(Pachrom_block_list)))

      return(list(rsids_gwas, beta_squared, sub_ld, chrom_block))
    })
    close(pb)
    # cbind_df代替 rbind_df(list_df)
    rm(Pa_chrom_block)
    pa_block <- list()
    pa_block$block_info <- bigreadr::rbind_df(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[4]]
    ))
    pa_block$snps <- bigreadr::rbind_df(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[1]]
    ))
    pa_block$y <- unlist(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[2]]
    ))
    # message(paste0("snp is ",nrow(snp_data)))
    pa_block$ld_data <- bigreadr::rbind_df(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[3]]
    ))

    rsid_x <- intersect(pa_block$snps$rsid, unique(unlist(pa_block$ld_data[, 1:2])))

    pa_block$ld_data <- as.data.frame(pa_block$ld_data[pa_block$ld_data$SNP_A
      %in% rsid_x & pa_block$ld_data$SNP_B
      %in% rsid_x, ])
    rm(Pa_chrom_data)
    # message(paste0("ld is ",nrow(sub_ld)))
    if (nrow(pa_block$ld_data) == 0) {
      ld_matrix <- diag(1, nrow = nrow(pa_block$snps))
    } else {
      ld_matrix <- make_ld_matrix(
        all_snps = pa_block$snps$rsid,
        ld_data = pa_block$ld_data
      )
    }

    pa_block$n_snps <- nrow(pa_block$snps)
    pa_block$ld_matrix_squared <- ld_matrix * ld_matrix
    return(pa_block)
  })
  names(pa_blocksnplist) <- names(Pachrom_block_list)
  return(pa_blocksnplist)
}
