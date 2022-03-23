
#' Link_pathway_blocks_gwas
#'
#' @description Link pathway blocks region and gwas summary snps
#' Takes block information, potentially independent LD pathway_blocks
#' or gene blocks,
#' @param Pagwas Pagwas format, deault is NULL.
#' @param chrom_ld LD data for 22 chromosome.
#' @param n.cores Parallel cores,default is 1. use detectCores()
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
                                     n.cores = 1) {

  # nested lists of chrom and then individual pathway_blocks
  Pachrom_block_list <- lapply(Pagwas$pathway_blocks, function(pa_blocks) split(pa_blocks, f = as.vector(pa_blocks$chrom)))

  names(Pachrom_block_list) <- names(Pagwas$pathway_blocks)
  Pagwas$pathway_blocks<-NULL
  #rm(pathway_blocks)

  chrom_gwas_list <- lapply(split(Pagwas$gwas_data, f = Pagwas$gwas_data$chrom), function(gwas) {
      gwas <- data.table::data.table(gwas)
      data.table::setkey(gwas, pos)
      return(gwas)
    }
  )
  Pagwas$gwas_data<-NULL
  #rm(gwas_data)
  # prevent naming issues and indexing issues

  message(paste0("* Start to link gwas and pathway block annotations for ",
                 length(Pachrom_block_list), " pathways!"))
  pb <- txtProgressBar(style = 3)

  Pathway_ld_gwas_data <- papply(names(Pachrom_block_list), function(pathway) {

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

      rsids <- snp_gene_df[which(snp_gene_df$label %in%
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

    snp_data <- Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),
                       lapply(seq_len(length(Pa_chrom_data)),
                              function(i) Pa_chrom_data[[i]][[1]]))

    # message(paste0("snp is ",nrow(snp_data)))
    sub_ld <- as.data.frame(Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),
                                   lapply(seq_len(length(Pa_chrom_data)),
                                          function(i) Pa_chrom_data[[i]][[3]])))

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


    setTxtProgressBar(pb, which(names(Pachrom_block_list) == pathway) / length(names(Pachrom_block_list)))

    return(list(
      block_info = Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),
                          lapply(seq_len(length(Pa_chrom_data)),
                                 function(i) Pa_chrom_data[[i]][[4]])),
      snps = snp_data,
      ld_data = sub_ld,
      n_snps = nrow(snp_data),
      ld_matrix_squared = ld_matrix * ld_matrix,
      y = unlist(lapply(seq_len(length(Pa_chrom_data)),
                        function(i) Pa_chrom_data[[i]][[2]]))
    ))
  }, n.cores = n.cores)

  close(pb)

  names(Pathway_ld_gwas_data) <- names(Pachrom_block_list)
  rm(chrom_gwas_list)
  rm(chrom_ld)
  #SOAR::Store(Pathway_ld_gwas_data)
  message("*** Start to store the variables: Pathway_ld_gwas_data")
  SOAR::Store(Pathway_ld_gwas_data)
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
  ld_matrix <- as(ld_matrix, "dgCMatrix")

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

  #storage.mode(ld_matrix)<-"integer"
  return(ld_matrix)
}
