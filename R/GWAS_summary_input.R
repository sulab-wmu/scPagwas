#' GWAS_summary_input
#' @description Input and progress the GWAS summary statistics data
#' filter the maf and the sex chrom data
#' @param Pagwas Pagwas format, deault is NULL.
#' @param gwas_data A data set object, need a data frame format,
#' There must be have colomue with 'chrom', 'pos', 'rsid', 'beta', 'se',
#' 'maf'
#' @param maf_filter Filter the maf, default is <0.01
#' @param Sex_filter Filter the SNPs in Sex chromosome,  default is TRUE.
#' @param MHC_filter Filter the SNPs in MHC chromosome,  default is TRUE.
#' @param gwas_z_filter Filter the z-score, calculated by abs(beta/se)
#'
#' @return
#' Returns a list:
#' \item{gwas_data}{filtered gwas summary data frame.}
#' @export
#'
#' @examples
#' library(scPagwas)
#' Pagwas <- list()
#' gwas_data <- bigreadr::fread2(system.file("extdata",
#'   "GWAS_summ_example.txt",
#'   package = "scPagwas"
#' ))
#' Pagwas <- GWAS_summary_input(Pagwas = Pagwas, gwas_data = gwas_data)
#' head(Pagwas$gwas_data)
#' # chrom     pos       rsid        se         beta     maf
#' # 1  chr1 1138913  rs3819001 0.0724484 4.774317e-01 0.06761
#' # 2  chr1 1691050 rs35672141 0.0368296 1.108626e-01 0.48200
#' # 3  chr1 1722932  rs3737628 0.0324817 7.881858e-02 0.48820
#' # 4  chr1 2369498  rs4592207 0.0399954 8.822031e-02 0.41790
#' # 5  chr1 3259369  rs2483286 0.0352555 1.541813e-03 0.30470
#' # 6  chr1 3811790 rs12085743 0.0405247 2.465163e-05 0.18500
#' @author Chunyu Deng
#' @aliases GWAS_summary_input
#' @keywords GWAS_summary_input, filter the gwas summary data.

GWAS_summary_input <- function(Pagwas = NULL,
                               gwas_data = NULL,
                               maf_filter = 0.01,
                               Sex_filter = TRUE,
                               MHC_filter = TRUE,
                               gwas_z_filter = -1) {
  message("Input gwas summary data frame!")
  maf <- chrom <- pos <- se <- abs_z <- NULL
  if (class(gwas_data)[1] != "data.frame") {
    stop("It is not a data frame format!")
  }

  necessary_cols <- c("chrom", "pos", "rsid", "beta", "se", "maf")
  # try to control for maf

  if (all(!(necessary_cols %in% colnames(gwas_data)))) {
    stop("There are missing the colomn for chrom, pos, beta, se, maf")
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
    gwas_data$chrom <- paste("chr", gwas_data$chrom, sep = "")
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
    gwas_data_6 <- gwas_data %>%
      dplyr::filter(chrom %in% "chr6") %>%
      dplyr::filter(pos > 25000000 & pos < 34000000)

    gwas_data <- gwas_data[!(gwas_data$chrom %in% "chr6"), ]
    gwas_data <- merge(gwas_data, gwas_data_6, all = TRUE)
    rm(gwas_data_6)
  }

  if (gwas_z_filter > 0 ) {
    message(paste("removing SNPs with |z| > ", gwas_z_filter, sep = ""))
    gwas_data <- gwas_data %>%
      dplyr::mutate(abs_z = abs(beta / se)) %>%
      dplyr::filter(abs_z < gwas_z_filter)
  }
  gwas_data<-gwas_data[gwas_data$beta<1,]
  Pagwas$gwas_data <- gwas_data
  return(Pagwas)
}
