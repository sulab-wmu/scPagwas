
#' link_pwpca_block
#' @description Link the pca score and expression for each pathway genes
#' for each block
#' Requires rownames that are identitcal to block labels loaded previously.
#' @param Pagwas Pagwas format, deault is NULL.
#' @param scPagwasSession "scPagwasSession"
#'
#' @return
#' @export
#' @examples
#' library(scPagwas)
#' # Pagwas should have inhibit data
#' Pagwas <- link_pwpca_block(Pagwas)

link_pwpca_block <- function(Pagwas,scPagwasSession="scPagwasSession") {
  Sys.setenv(R_LOCAL_CACHE=scPagwasSession)
  cell_names <- intersect(colnames(merge_scexpr), colnames(pca_cell_df))

  merge_scexpr <- merge_scexpr[, cell_names]
  pca_cell_df <- data.matrix(pca_cell_df[, cell_names])
  #Pagwas$pca_cell_df <- pca_cell_df

  message("*  merging functional information about blocks")
  pb <- txtProgressBar(style = 3)

  Pathway_ld_gwas_data <- lapply(Pathway_ld_gwas_data, function(pa_block) {

    pathway <- unique(pa_block$block_info$pathway)

    x <- pca_cell_df[pathway, ]
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
    mg <- intersect(Pagwas$rawPathway_list[[pathway]], proper_genes)
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

    pa_block$x <-Matrix::crossprod(t(as(pa_block$ld_matrix_squared,"matrix")), x3)
    rownames(pa_block$x) <- pa_block$snps$rsid
    colnames(pa_block$x) <- colnames(merge_scexpr)
    rm(x3)
    pa_block$include_in_inference <- T

    setTxtProgressBar(pb, which(names(Pathway_ld_gwas_data) == pathway) / length(Pathway_ld_gwas_data))

    return(pa_block)
  })
  close(pb)
  CT_Pathway_ld_gwas_data <- Pathway_ld_gwas_data[!sapply(Pathway_ld_gwas_data, is.null)]

  message("*** Start to store the variables: ")
  message("*1)CT_Pathway_ld_gwas_data")
  SOAR::Store(CT_Pathway_ld_gwas_data)
  message("*2)merge_scexpr")
  SOAR::Store(merge_scexpr)
  message("*3)pca_cell_df")
  SOAR::Store(pca_cell_df)
  gc()
  return(Pagwas)
}

