#' Functions for mapping Snps to Genes
#' @description Maps SNPs to their nearest genes within TSS windows.
#'
#' @param gwas_data (data.frame) the result for GWAS_summary_input
#' @param block_annotation (data.frame) Start and end points for block
#'  traits, usually genes.
#' @param marg (integer) region upstream and downstream(default=10000).
#'
#' @return data frame for gwas_data mapping to gene.
#' @export
#'
#' @examples
#' library(scPagwas)
#' Pagwas <- list()
#' gwas_data <- bigreadr::fread(system.file("extdata",
#'   "GWAS_summ_example.txt",
#'   package = "scPagwas"
#' ))
#' Pagwas <- GWAS_summary_input(
#'   Pagwas = Pagwas,
#'   gwas_data = gwas_data,
#'   maf_filter = 0.1
#' )
#' Pagwas$gwas_data_gene_df <- SnpToGene(
#'   gwas_data = Pagwas$gwas_data,
#'   block_annotation = block_annotation,
#'   marg = 10000
#' )
SnpToGene <- function(gwas_data, block_annotation, marg = 10000) {
  gwas_data_geno_range <- GenomicRanges::GRanges(gwas_data[, "chrom"],
    IRanges::IRanges(
      as.numeric(gwas_data[, "pos"]),
      as.numeric(gwas_data[, "pos"])
    ),
    name = gwas_data[, "rsid"]
  )

  gene_geno_range <- GenomicRanges::GRanges(block_annotation[, "chrom"],
    IRanges::IRanges(
      block_annotation[, "start"] + 1 - marg,
      block_annotation[, "start"] + marg
    ),
    name = block_annotation[, "label"]
  )

  start_geno_range <- GenomicRanges::resize(gene_geno_range, fix = "start", width = 1L)
  end_geno_range <- GenomicRanges::resize(gene_geno_range, fix = "end", width = 1L)

  orig_distance <- GenomicRanges::distanceToNearest(gwas_data_geno_range,
    gene_geno_range,
    ignore.strand = TRUE
  )

  distance_n <- orig_distance@elementMetadata$distance
  distance_within <- data.frame(
    stat_dist = orig_distance@from,
    end_dist = orig_distance@to,
    distance = distance_n
  )

  index <- which(distance_n == 0)

  map_result1 <- data.frame(
    rsid=gwas_data_geno_range$name[distance_within$stat_dist[index]],
    label=gene_geno_range$name[distance_within$end_dist[index]],
    pos=as.numeric(gwas_data[distance_within$stat_dist[index], "pos"]),
    Disstance=distance_within$distance[index]
  )
  map_result1 <- map_result1[!is.na(map_result1$rsid),]
  rownames(map_result1) <- map_result1$rsid
  return(map_result1)
}
