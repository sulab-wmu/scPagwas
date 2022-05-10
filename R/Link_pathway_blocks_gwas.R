
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
                                     split_n=1,
                                     ncores = 1) {

  options(bigmemory.allow.dimnames=TRUE)
  Pachrom_block_list <- lapply(Pagwas$pathway_blocks, function(pa_blocks) split(pa_blocks, f = as.vector(pa_blocks$chrom)))

  names(Pachrom_block_list) <- names(Pagwas$pathway_blocks)
  Pagwas$pathway_blocks<-NULL

  chrom_gwas_list <- lapply(split(Pagwas$gwas_data, f = Pagwas$gwas_data$chrom), function(gwas) {
      gwas <- data.table::data.table(gwas)
      data.table::setkey(gwas, pos)
      return(gwas)
    }
  )
  Pagwas$gwas_data<-NULL
  #rm(gwas_data)
  # prevent naming issues and indexing issues
  Pathway_sclm_results<-list()
  #+++++++Pathway_lm_results<-list()
  Pathway_ld_gwas_data<-list()

  message(paste0("* Start to link gwas and pathway block annotations for ",
                 length(Pachrom_block_list), " pathways!"))
  options(bigmemory.allow.dimnames=TRUE)
  pb <- txtProgressBar(style = 3)

  #Pathway_ld_gwas_data <- papply(names(Pachrom_block_list), function(pathway) {
  for (pathway in names(Pachrom_block_list)) {

 #   }
    # message(paste(' - starting blocks on pathway: ', pa, sep = ''))
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
    #cbind_df代替 rbind_df(list_df)
    snp_data <- bigreadr::rbind_df(lapply(seq_len(length(Pa_chrom_data)),
                              function(i) Pa_chrom_data[[i]][[1]]))

    # message(paste0("snp is ",nrow(snp_data)))
    sub_ld <- bigreadr::rbind_df(lapply(seq_len(length(Pa_chrom_data)),
                                          function(i) Pa_chrom_data[[i]][[3]]))

    rsid_x <- intersect(snp_data$rsid, unique(unlist(sub_ld[, 1:2])))

    sub_ld <- as.data.frame(sub_ld[sub_ld$SNP_A
                                   %in% rsid_x & sub_ld$SNP_B
                                   %in% rsid_x, ])

    # message(paste0("ld is ",nrow(sub_ld)))
    if (nrow(sub_ld) == 0) {
      ld_matrix <- diag(1, nrow = nrow(snp_data))
    } else {
      ld_matrix <- make_ld_matrix(all_snps = snp_data$rsid, ld_data = sub_ld)
    }
    block_info<- bigreadr::rbind_df(lapply(seq_len(length(Pa_chrom_data)),
                           function(i) Pa_chrom_data[[i]][[4]]))


    pa_block<-list(
      block_info = block_info,
      snps = snp_data,
      ld_data = sub_ld,
      n_snps = nrow(snp_data),
      ld_matrix_squared = ld_matrix * ld_matrix,
      y = unlist(lapply(seq_len(length(Pa_chrom_data)),
                        function(i) Pa_chrom_data[[i]][[2]]))
    )

    ##compute the single cell result
    if(!is.null(split_n) & split_n>1 & class(split_n)== "numeric"){

      a<-ncol(Pagwas$pca_scCell_mat)

      for (ai in 1:10) {
        #print(ai)
        if(a%%split_n==0) break;
        split_n <- split_n+1
      }

      la<-gl(split_n,a/split_n,length=a)
      Pathway_sclm_part<-list()
      for(i in 1:split_n){
        Pathway_sclm_part[[i]] <- get_Pathway_sclm(pa_block=pa_block,
                                                      pca_scCell_mat=Pagwas$pca_scCell_mat[,which(la==i)],
                                                      data_mat=Pagwas$data_mat[,which(la==i)],
                                                      rawPathway_list=Pagwas$rawPathway_list,
                                                      snp_gene_df=Pagwas$snp_gene_df,
                                                      ncores=ncores )

      }


      Pathway_sclm_results[[pathway]] <- unlist(Pathway_sclm_part)[colnames(Pagwas$pca_scCell_mat)]

    }else{
      Pathway_sclm_results[[pathway]]<-get_Pathway_sclm(pa_block=pa_block,
                                             pca_scCell_mat=Pagwas$pca_scCell_mat,
                                             data_mat=Pagwas$data_mat,
                                             rawPathway_list=Pagwas$rawPathway_list,
                                             snp_gene_df=Pagwas$snp_gene_df,
                                             ncores=ncores)
    }

      pa_block<-link_pwpca_block(pa_block=pa_block,
                                 pca_cell_df=data.matrix(Pagwas$pca_cell_df),
                                 merge_scexpr=Pagwas$merge_scexpr,
                                 snp_gene_df=Pagwas$snp_gene_df,
                                 rawPathway_list=Pagwas$rawPathway_list)
    #Pathway_lm_results[[pathway]]<- Pa_Pagwas_perform_regression(pa_block=pa_block)
    Pathway_ld_gwas_data[[pathway]]<-pa_block[c("x","y","snps","include_in_inference")]
    setTxtProgressBar(pb, which(names(Pachrom_block_list) == pathway) / length(names(Pachrom_block_list)))

  }

  close(pb)

  rm(chrom_gwas_list)
  rm(chrom_ld)
  Pathway_sclm_results <- Pathway_sclm_results[!sapply(Pathway_sclm_results, is.null)]
  Pathway_sclm_results <- data.matrix(as.data.frame(Pathway_sclm_results))
  rownames(Pathway_sclm_results)<- colnames(Pagwas$pca_scCell_mat)
  Pagwas$Pathway_sclm_results<-as(Pathway_sclm_results,"dgCMatrix")

  #Pathway_lm_results <- Pathway_lm_results[!sapply(Pathway_lm_results, is.null)]
  #Pathway_lm_results <- data.matrix(as.data.frame(Pathway_lm_results))
  #rownames(Pathway_lm_results)<- colnames(Pagwas$pca_cell_df)
  #Pagwas$Pathway_lm_results<-as(Pathway_lm_results,"dgCMatrix")

  names(Pathway_ld_gwas_data) <- names(Pachrom_block_list)
  Pagwas$Pathway_ld_gwas_data<-Pathway_ld_gwas_data
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
