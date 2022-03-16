

#' GWAS_summary_input
#' @description Input and progress the GWAS summary statistics data
#' filter the maf and the sex chrom data
#' @param Pagwas Pagwas format, deault is NULL.
#' @param gwas_data A data set object, need a data frame format, There must be have colomue with 'chrom', 'pos', 'rsid', 'beta', 'se', 'maf'
#' @param maf_filter Filter the maf, default is <0.01
#' @param Sex_filter Filter the SNPs in Sex chromosome,  default is TRUE.
#' @param MHC_filter Filter the SNPs in MHC chromosome,  default is TRUE.
#' @param gwas_z_filter Filter the z-score, calculated by abs(beta/se)
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(GWAS_summ_example)
#' Pagwas <- GWAS_summary_input(Pagwas = NULL, gwas_data = GWAS_summ_example)
GWAS_summary_input <- function(Pagwas=NULL,
                               gwas_data=NULL,
                               maf_filter = 0.01,
                               Sex_filter = TRUE,
                               MHC_filter = TRUE,
                               gwas_z_filter = -1) {
  message("Input gwas summary data frame!")

  if (class(gwas_data)[1] != "data.frame") {
    stop("It is not a data frame format!")
  }

  necessary_cols <- c("chrom", "pos", "rsid", "beta", "se", "maf")
  # try to control for maf

  if (all(!(necessary_cols %in% colnames(gwas_data)))) {
    stop("There are missing the colomn for chrom, pos, beta, or se, maf")
  }

  if (class(gwas_data$pos)[1] != "numeric") {
    gwas_data$pos <- as.numeric(gwas_data$pos)
  }

  if (!is.numeric(gwas_z_filter) | !is.integer(gwas_z_filter)) {
    gwas_z_filter <- as.numeric(gwas_z_filter)
  }

  if (!is.numeric(maf_filter) | !is.integer(maf_filter)) {
    maf_filter <- as.numeric(maf_filter)
  }

  if (!grepl("chr", gwas_data$chrom[1])) {
    message('No "chr" from chrom!, now pasting it!')
    gwas_data$chrom <- paste("chr",gwas_data$chrom,sep = "")

  }

  if ("maf" %in% colnames(gwas_data)) {
    message("Filtering out SNPs with MAF criterion")
    gwas_data <- gwas_data %>%
      dplyr::mutate(maf = ifelse(test = maf > 0.5, yes = 1 - maf, no = maf)) %>%
      dplyr::filter(maf > maf_filter)
  }
  # remove duplicated!
  gwas_data <- gwas_data[!duplicated(gwas_data$rsid), ]

  if (Sex_filter) {
    gwas_data <- gwas_data[!(gwas_data$chrom %in% "chr23"), ]
    gwas_data <- gwas_data[!(gwas_data$chrom %in% "chrX"), ]
    gwas_data <- gwas_data[!(gwas_data$chrom %in% "chrY"), ]
  }

  if (MHC_filter) {
    gwas_data_6 <- gwas_data %>% dplyr::filter(chrom %in% "chr6") %>%
      dplyr::filter(pos > 25000000 & pos < 34000000)

    gwas_data <- gwas_data[!(gwas_data$chrom %in% "chr6"), ]
    gwas_data <- merge(gwas_data, gwas_data_6, all = TRUE)
    rm(gwas_data_6)
  }

  if (gwas_z_filter > 0) {

    message(paste("removing SNPs with |z| > ", gwas_z_filter, sep = ""))
    gwas_data <- gwas_data %>%
      dplyr::mutate(abs_z = abs(beta / se)) %>%
      dplyr::filter(abs_z < gwas_z_filter)
  }

  #Pagwas$gwas_data <- gwas_data
  SOAR::Store(gwas_data)
  gc()
  return(Pagwas)
}
