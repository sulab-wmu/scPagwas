
#' link_pwpca_block
#' @description Link the pca score and expression for each pathway genes for each block
#' Requires rownames that are identitcal to block labels loaded previously.
#' @param Pagwas Pagwas format, deault is NULL.
#'
#' @return
#' @export
#' @examples
#' library(scPagwas)
#' # Pagwas should have inhibit data
#' Pagwas <- link_pwpca_block(Pagwas)
link_pwpca_block <- function(Pagwas) {
  if (is.null(Pagwas$merge_scexpr)) {
    message("* no load merge_scexpr data in single cell data input step!")
    return(Pagwas)
  }
  merge_scexpr <- Pagwas$merge_scexpr

  if (is.null(Pagwas$Pathway_ld_gwas_data)) {
    message("* no loaded Pathway_ld_gwas data")
    return(Pagwas)
  }

  pca_cell_df <- Pagwas$pca_cell_df
  cell_names <- intersect(colnames(merge_scexpr), colnames(pca_cell_df))

  merge_scexpr <- merge_scexpr[, cell_names]
  pca_cell_df <- pca_cell_df[, cell_names]
  Pagwas$pca_cell_df <- pca_cell_df

  # pca_cell_df2 <- pca_cell_df %>% as.data.frame() %>% dplyr::mutate(partition = cut(row_number(), breaks = 10, labels = F))

  message("*  merging functional information about blocks")
  pb <- txtProgressBar(style = 3)

  Pathway_ld_gwas_data <- lapply(Pagwas$Pathway_ld_gwas_data, function(pa_block) {

    pathway <- unique(pa_block$block_info$pathway)

    # print(pathway)
    x <- pca_cell_df[pathway, ]
    # if (nrow(x) != 1) {stop("remove duplicates from pa_block data")}
    # pa_block$partition  <- pca_cell_df2[pathway,"partition"]

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
      x2 <- as.data.frame(apply(x2, 2, function(x) (x - min(x)) / (max(x) - min(x))))
      # x2 <- as.data.frame(apply(x2,2, function(x) x/sum(x)))
    }

    if (pa_block$n_snps > 1) {
      x2 <- Matrix::Matrix(as.matrix(x2[pa_block$snps$label, ]))
      pa_block$n_snps <- nrow(pa_block$snps)
      # pa_block$x2 <- x2
      x <- Matrix::Matrix(as.matrix(x[rep(1, pa_block$n_snps), ], drop = FALSE))
      rownames(x) <- pa_block$snps$rsid

      snp_gene_df <- Pagwas$snp_gene_df
      rownames(snp_gene_df) <- Pagwas$snp_gene_df$rsid
      x <- x * snp_gene_df[pa_block$snps$rsid, "slope"]
      x3 <- Matrix::Matrix(x2 * x)
    } else {
      x2 <- matrix(x2[pa_block$snps$label, ], nrow = 1)
      rownames(x2) <- pa_block$snps$label
      pa_block$n_snps <- nrow(pa_block$snps)
      x <- matrix(x[rep(1, pa_block$n_snps), ], nrow = 1)
      rownames(x) <- pa_block$snps$rsid
      snp_gene_df <- Pagwas$snp_gene_df
      rownames(snp_gene_df) <- Pagwas$snp_gene_df$rsid
      x <- matrix(as.numeric(x) * as.numeric(snp_gene_df[pa_block$snps$rsid, "slope"]), nrow = 1)
      x3 <- matrix(as.numeric(x2) * as.numeric(x), nrow = 1)
    }

    rownames(x3) <- pa_block$snps$rsid
    pa_block$x <- as.matrix(crossprod(t(pa_block$ld_matrix_squared), x3))
    pa_block$include_in_inference <- T
    #message("1")
    setTxtProgressBar(pb, which(names(Pagwas$Pathway_ld_gwas_data) == pathway) / length(Pagwas$Pathway_ld_gwas_data))

    return(pa_block)
  })
  close(pb)
  Pagwas$Pathway_ld_gwas_data <- Pathway_ld_gwas_data[!sapply(Pathway_ld_gwas_data, is.null)]

  return(Pagwas)
}

