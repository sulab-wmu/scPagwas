
#' Link_pathway_blocks_gwas
#'
#' @description Link pathway blocks region and gwas summary snps
#' Takes block information, potentially independent LD pathway_blocks
#' or gene blocks,
#'
#' @param Pagwas Pagwas format, deault is NULL.
#' @param chrom_ld LD data for 22 chromosome.
#' @param singlecell whether to run singlecell progress
#' @param celltype whether to run celltype progress
#' @param n.cores cores for regression
#' @param backingpath file address for bk files, no "/" in the end.
#' @return Pagwas list of blocks for snp to gene and pathway score.
#' @export
#' @author Chunyu Deng
#' @aliases Link_pathway_blocks_gwas
#' @keywords Link_pathway_blocks_gwas, calculate the svd matrix for
#' singlecell data.

Link_pathway_blocks_gwas <- function(Pagwas,
                                     chrom_ld = NULL,
                                     singlecell = T,
                                     celltype = T,
                                     backingpath,
                                     n.cores=1) {
  options(bigmemory.allow.dimnames = TRUE)
  pos <- NULL

  Pachrom_block_list <- lapply(
    Pagwas$pathway_blocks,
    function(pa_blocks) {
      split(pa_blocks,
        f = as.vector(pa_blocks$chrom)
      )
    }
  )

  names(Pachrom_block_list) <- names(Pagwas$pathway_blocks)
  Pagwas$pathway_blocks <- NULL

  chrom_gwas_list <- lapply(split(Pagwas$gwas_data,
    f = Pagwas$gwas_data$chrom
  ), function(gwas) {
    gwas <- data.table::data.table(gwas)
    data.table::setkey(gwas, pos)
    return(gwas)
  })
  Pagwas$gwas_data <- NULL


  Pagwas <- Pathway_block_func(
    Pagwas = Pagwas,
    Pachrom_block_list = Pachrom_block_list,
    chrom_ld = chrom_ld,
    chrom_gwas_list = chrom_gwas_list,
    singlecell = singlecell,
    celltype = celltype,
    backingpath=backingpath,
    n.cores=n.cores
  )


  rm(chrom_gwas_list)
  rm(chrom_ld)


  gc()
  return(Pagwas)
}

#' Pathway_block_func
#'
#' @param Pagwas NULL
#' @param Pachrom_block_list result for chrome block
#' @param chrom_gwas_list list of gwas block in chrom
#' @param singlecell whether to run singlecell progress
#' @param celltype whether to run celltype progress
#' @param chrom_ld default is block_annotation
#' @param n.cores cores for regression
#' @param backingpath file address for bk files, no "/" in the end.
#' @return Pagwas result list including Pathway block.
#'
Pathway_block_func <- function(Pagwas = NULL,
                               Pachrom_block_list,
                               chrom_gwas_list,
                               singlecell = T,
                               celltype = T,
                               chrom_ld,
                               backingpath,
                               n.cores=1) {
  Pathway_sclm_results <- list()
  Pathway_ld_gwas_data <- list()
  SNP_A <- NULL
  message(paste0(
    "* Start to link gwas and pathway block annotations for ",
    length(Pachrom_block_list), " pathways!"
  ))
  options(bigmemory.allow.dimnames = TRUE)
  pb <- txtProgressBar(style = 3)

  #Rn<-randomStrings(n=length(Pachrom_block_list),len=9,digits=TRUE,
  #                          upperalpha=F,
  #                          loweralpha=TRUE,
  #                          unique=TRUE,
  #                          check=F)
  #Rn<-Rn[,1]
  #names(Rn)<-names(Pachrom_block_list)

  for (pathway in names(Pachrom_block_list)) {
    Pa_chrom_block <- Pachrom_block_list[[pathway]]

    Pa_chrom_data <- lapply(names(Pa_chrom_block), function(chrom) {
      chrom_block <- Pa_chrom_block[[chrom]]
      ld_data <- chrom_ld[[chrom]]

      data.table::setkey(ld_data, SNP_A) # just in case
      if (!(chrom %in% names(chrom_gwas_list))) {
        warning(paste(chrom, " for gwas is missing, could be a problem!",
          sep = ""
        ))
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

    pa_block$ld_data <- bigreadr::rbind_df(lapply(
      seq_len(length(Pa_chrom_data)),
      function(i) Pa_chrom_data[[i]][[3]]
    ))

    rsid_x <- intersect(pa_block$snps$rsid, unique(
      unlist(pa_block$ld_data[, 1:2])
    ))

    pa_block$ld_data <- as.data.frame(
      pa_block$ld_data[pa_block$ld_data$SNP_A
        %in% rsid_x & pa_block$ld_data$SNP_B
        %in% rsid_x, ]
    )
    rm(Pa_chrom_data)

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

    ## compute the single cell result

    Rn<-paste0(pathway,sample(1:1000,1))
    if (singlecell) {

      Pathway_sclm_results[[pathway]] <- get_Pathway_sclm(
        pa_block = pa_block,
        pca_scCell_mat = Pagwas$pca_scCell_mat,
        data_mat = Pagwas$data_mat,
        rawPathway_list = Pagwas$rawPathway_list,
        snp_gene_df = Pagwas$snp_gene_df,
        backingpath=backingpath,
        n.cores=n.cores,
        Rns=Rn
      )
    }

    if (celltype) {
      pa_block <- link_pwpca_block(
        pa_block = pa_block,
        pca_cell_df = data.matrix(Pagwas$pca_cell_df),
        merge_scexpr = Pagwas$merge_scexpr,
        snp_gene_df = Pagwas$snp_gene_df,
        rawPathway_list = Pagwas$rawPathway_list
      )


      Pathway_ld_gwas_data[[pathway]] <- pa_block
    }
    setTxtProgressBar(
      pb,
      which(names(Pachrom_block_list) == pathway) /
        length(names(Pachrom_block_list))
    )
  }

  close(pb)

  if (singlecell) {
    Pathway_sclm_results <- Pathway_sclm_results[!sapply(
      Pathway_sclm_results, is.null
    )]
    Pathway_sclm_results <- data.matrix(as.data.frame(Pathway_sclm_results))
    rownames(Pathway_sclm_results) <- colnames(Pagwas$pca_scCell_mat)
    Pagwas$Pathway_sclm_results <- Pathway_sclm_results#, "dgCMatrix")
  }
  if (celltype) {
    names(Pathway_ld_gwas_data) <- names(Pachrom_block_list)
    Pagwas$Pathway_ld_gwas_data <- Pathway_ld_gwas_data
  }
  return(Pagwas)
}


#' link_pwpca_block
#'
#' @param pa_block snp blocks for pathway.
#' @param pca_cell_df pca matrix for cell types.
#' @param merge_scexpr mean expression for merged cell types of single cell
#' data.
#' @param snp_gene_df Data frame of snp mapping to gene.
#' @param rawPathway_list Raw pathway list before filter
#'
#' @description Link the pca score and expression for each pathway genes
#' for each block
#' Requires rownames that are identitcal to block labels loaded previously.
#' @return pathway block including "x" "y" "snps" "include_in_inference"

link_pwpca_block <- function(pa_block,
                             pca_cell_df,
                             merge_scexpr,
                             snp_gene_df,
                             rawPathway_list) {
  pathway <- unique(pa_block$block_info$pathway)

  x <- pca_cell_df[pathway, ]
  if (length(pathway) == 1) {
    x <- matrix(x, nrow = 1)
    rownames(x) <- pathway
  }

  if (nrow(pa_block$snps) == 0) {
    pa_block$include_in_inference <- F
    pa_block$x <- NULL # to make sure we totally replace previous stuffs
    return(pa_block)
  }
  proper_genes <- rownames(merge_scexpr)
  mg <- intersect(rawPathway_list[[pathway]], proper_genes)
  x2 <- merge_scexpr[mg, ]

  if (ncol(merge_scexpr) == 1) {
    x2 <- data.matrix(x2)
  }
  if (length(mg) > 1) {
    x2 <- apply(x2, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  }

  if (pa_block$n_snps > 1) {
    if (ncol(merge_scexpr) == 1) {
      x2 <- data.matrix(x2[pa_block$snps$label, ])
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- data.matrix(x[rep(1, pa_block$n_snps), ])
      rownames(x) <- pa_block$snps$rsid


      rownames(snp_gene_df) <- snp_gene_df$rsid
      #x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
      x3 <- x2 * x
    } else {
      x2 <- x2[pa_block$snps$label, ]
      pa_block$n_snps <- nrow(pa_block$snps)

      x <- x[rep(1, pa_block$n_snps), ]
      rownames(x) <- pa_block$snps$rsid


      rownames(snp_gene_df) <- snp_gene_df$rsid
      #x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
      x3 <- x2 * x
    }
  } else {
    x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
    rownames(x2) <- pa_block$snps$label
    pa_block$n_snps <- nrow(pa_block$snps)

    rownames(x) <- pa_block$snps$rsid

    rownames(snp_gene_df) <- snp_gene_df$rsid

    #x <- matrix(as.numeric(x) * as.numeric(
    #  snp_gene_df[pa_block$snps$rsid, "slope"]
    #), nrow = 1)
    x3 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
  }

  rm(x)
  rm(x2)

  pa_block$x <- t(pa_block$ld_matrix_squared) %*% x3
  rownames(pa_block$x) <- pa_block$snps$rsid
  colnames(pa_block$x) <- colnames(merge_scexpr)
  rm(x3)

  pa_block$include_in_inference <- T

  return(pa_block[c("x", "y", "snps", "include_in_inference")])
}


#' make_ld_matrix
#' @description construct the final matrix for regression in blocks
#' @param all_snps The snps that were queried
#' @param ld_data A returned LD matrix with SNP
#'
make_ld_matrix <- function(all_snps, ld_data) {
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
