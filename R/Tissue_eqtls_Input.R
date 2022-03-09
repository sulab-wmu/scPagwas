
#' Tissue_eqtls_Input
#' @description Get the eqtls data and interact with GWAS summary data.
#'
#' @param Pagwas Pagwas format data list inherit from other functions.
#' @param add_eqtls There is options: "OnlyEqtls","Both"; "OnlyEqtls" means only eqtls snps will choose; "Both" means both eqtls and gene-TSS-window snps will choose;
#' @param eqtls_files eqtls files names,frome GTEx.
#' @param eqtl_p pvalue threshold for eqtls data.
#' @param eqtls_cols c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name","slope")
#' @param marg gene-TSS-window size
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(eqtls_files)
#'
#' Pagwas <- Tissue_eqtls_Input(Pagwas = Pagwas, add_eqtls = "OnlyEqtls", eqtls_files = eqtls_files)
Tissue_eqtls_Input <- function(Pagwas = Pagwas,
                               add_eqtls = "OnlyEqtls",
                               eqtls_files = NULL,
                               eqtl_p = 0.05,
                               eqtls_cols = c("rs_id_dbSNP151_GRCh38p7", "variant_pos", "tss_distance", "gene_chr", "gene_start", "gene_end", "gene_name", "slope"),
                               marg = 10000) {
  if (length(eqtls_files) == 1) {
    eqtls_sig <- read.delim(eqtls_files)
  } else {
    eqtls_sig <- lapply(eqtls_files, read.delim) %>% Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2), .)
  }

  if (!is.character(eqtls_cols)) {
    stop("eqtls_cols is not of character!")
  }

  if (length(eqtls_cols) != 8) {
    stop("eqtls_cols have not 7 characters!")
  }

  if ("pval_true_df" %in% colnames(eqtls_sig)) {
    eqtls_sig <- eqtls_sig[eqtls_sig$pval_true_df < eqtl_p, eqtls_cols]
  }

  colnames(eqtls_sig) <- c("rsid", "pos", "Disstance", "gene_chr", "gene_start", "geng_end", "gene_name", "slope")

  eqtls_sig <- unique(eqtls_sig)
  inter_snps <- unique(intersect(Pagwas$gwas_data$rsid, eqtls_sig$rsid))

  message("There are ", length(inter_snps), " snps for significant eqtls!")

  if (length(inter_snps) == 0) {
    stop("The significant snp for eqtls were none! please set the add_eqtls = FALSE !")
  }
  eqtls_sig <- eqtls_sig[eqtls_sig$rsid %in% inter_snps, ]

  a2 <- Pagwas$gwas_data[!(Pagwas$gwas_data$rsid %in% inter_snps), ]

  # The 1st condition
  snp_gene_df1 <- unique(eqtls_sig[, c("rsid", "gene_name", "pos", "Disstance", "slope")])
  colnames(snp_gene_df1) <- c("rsid", "label", "pos", "Disstance", "slope")
  snp_gene_df1 <- snp_gene_df1[!duplicated(snp_gene_df1$rsid), ]

  if (add_eqtls == "OnlyEqtls") {
    Pagwas$snp_gene_df <- snp_gene_df1
    Pagwas$gwas_data <- Pagwas$gwas_data[Pagwas$gwas_data$rsid %in% inter_snps, ]
  } else if (add_eqtls == "Both") {

    # The 2nd condition
    snp_gene_df2 <- Snp2Gene(snp = a2, refGene = Pagwas$block_annotation, marg = marg)
    snp_gene_df2 <- snp_gene_df2[snp_gene_df2$Disstance == "0", ]
    snp_gene_df2$slope <- median(snp_gene_df1$slope)
    Pagwas$snp_gene_df <- rbind(snp_gene_df1, snp_gene_df2)
  }

  return(Pagwas)
}
