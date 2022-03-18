#' Genes_by_pathway_kegg
#'
#' @description
#' A few numbers from \code{
#' library(KEGGREST)
#' pathways.list <- keggList("pathway", "hsa")
#' # Pull all genes for each pathway
#' pathway.codes <- sub("path:", "", names(pathways.list))
#' Genes_by_pathway.kegg <- sapply(pathway.codes,function(pwid){
#' pw <- keggGet(pwid)
#' if (is.null(pw[[1]]$GENE)) return(NA)
#' pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)]
#' # may need to modify this to c(FALSE, TRUE) for other organisms
#' pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
#'  return(pw2)
#' })
#'
#' }
#'
#' @docType data
#' @name Genes_by_pathway_kegg
#' @usage data(Genes_by_pathway_kegg)
#' @format a list.
#' @keywords datasets
#' @examples
#' data(Genes_by_pathway_kegg)
#' str(Genes_by_pathway_kegg)
NULL



#' Genes_by_pathway_hallmark
#'
#' @description A hallmark pathway gene set from \code{
#' x <- readLines("E:/OneDrive/GWAS_Multiomics/tempdata/h.all.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' Genes_by_pathway_hallmark <- lapply(res, "[", -c(1:2))
#' save(Genes_by_pathway_hallmark,file="E:/RPakage/scPagwas/data/Genes_by_pathway.hallmark_RData")
#' }
#' @docType data
#' @name Genes_by_pathway_hallmark
#' @usage data(Genes_by_pathway_hallmark)
#' @format a list.
#' @keywords datasets
#' @examples
#' data(Genes_by_pathway_hallmark)
#' str(Genes_by_pathway_hallmark)
NULL

#' Genes_by_pathway_reactome
#'
#' @description A reactome pathway gene set from \code{
#' x <- readLines("E:/OneDrive/GWAS_Multiomics/tempdata/c2.cp.reactome.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' Genes_by_pathway.reactome <- lapply(res, "[", -c(1:2))
#' save(Genes_by_pathway.reactome,file="E:/RPakage/scPagwas/data/Genes_by_pathway.reactome.RData")
#' }
#'
#' @docType data
#' @name Genes_by_pathway_reactome
#' @usage data(Genes_by_pathway_reactome)
#' @format list
#' @keywords datasets
#' @examples data(Genes_by_pathway_reactome)
#' str(Genes_by_pathway_reactome)
NULL



#' block_annotation
#'
#' @description
#' A gene annotation files from \code{
#' library("rtracklayer")
#' gtf_df<- rtracklayer::import("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/gencode.v34.annotation.gtf.gz")
#' gtf_df <- as.data.frame(gtf)
#' gtf_df <- gtf_df[,c("seqnames","start","end","type","gene_name")]
#' gtf_df <- gtf_df[gtf_df$type=="gene",]
#' block_annotation<-gtf_df[,c(1,2,3,5)]
#' colnames(block_annotation)<-c("chrom", "start","end","label")
#' save(block_annotation,file = "/share/pub/dengcy/GWAS_Multiomics/pagwas/data/block_annotation.RData")
#' }
#'
#' @docType data
#' @name block_annotation
#' @usage data(block_annotation)
#' @format a data.frame,
#' @source Generated from gencode.v34.annotation.gtf.gz
#' @keywords datasets
#' @examples data(block_annotation)
#' str(block_annotation)
NULL


#' chrom_ld
#'
#' @description use vcftools to get ped and map files
#' LD data from \code{
#' /share/apps/vcftools/bin/vcftools --vcf /share/pub/dengcy/Singlecell/COVID19/1000genomes_all_genotypes.vcf --plink-tped --out /share/pub/dengcy/Singlecell/COVID19/1000genomes_all_genotypes
#'
#' /share/pub/dengcy/Singlecell/COVID19/PLINK/plink --tfile /share/pub/dengcy/Singlecell/COVID19/1000genomes_all_genotypes --recode --out /share/pub/dengcy/Singlecell/COVID19/1000genomes_all_genotypes
#' /share/pub/dengcy/Singlecell/COVID19/PLINK/plink --map /share/pub/dengcy/Singlecell/COVID19/1000genomes_all_genotypes.map --ped /share/pub/dengcy/Singlecell/COVID19/1000genomes_all_genotypes.ped --allow-no-sex --autosome --r2 --ld-window-kb 1000 --ld-window-r2 0.2 --out /share/pub/dengcy/Singlecell/COVID19/ld_1000genome
#' covid_ld<-read.delim("/share/pub/dengcy/Singlecell/COVID19/ld_1000genome.ld")
#' covid_ld<-covid_ld[!(covid_ld$%in% 23),]
#' colnames(covid_ld)[7]<-"R"
#' lapply(unique(covid_ld$CHR_A), function(i){
#'   a<-data.table(covid_ld[covid_ld$CHR_A == i,])
#'   file_name <- paste0("/share/pub/dengcy/Singlecell/COVID19/data/LD/",i,".Rds")
#'   saveRDS(a, file = file_name)
#' })
#'
#'
#' chrom_ld<-lapply(as.character(1:22),function(chrom){
#'   chrom_ld_file_path <- paste(ld_folder, '/', chrom, '.Rds', sep = '')
#'  ld_data <- readRDS(chrom_ld_file_path)[(R**2 > r2_threshold), .(SNP_A, SNP_B, R)]
#'   return(ld_data)
#' })
#' save(chrom_ld,file="/share/pub/dengcy/GWAS_Multiomics/pagwas/data/chrom_ld.RData")
#'
#' }
#'
#' @docType data
#' @name chrom_ld
#' @usage data(chrom_ld)
#' @format a list.
#' @source Generated from PLINK 1.90 linux
#' @keywords datasets
#' @examples
#' data(chrom_ld)
#' str(chrom_ld)
NULL
