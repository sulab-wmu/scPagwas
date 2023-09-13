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


#' genes.by.celltype.pathway
#'
#' @description A celltype pathway gene set from \code{
#' x <- readLines("c8.all.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' genes.by.celltype.pathway <- lapply(res, "[", -c(1:2))
#' }
#'
#' @docType data
#' @name genes.by.celltype.pathway
#' @usage data(genes.by.celltype.pathway)
#' @format list
#' @keywords datasets
#' @examples data(genes.by.celltype.pathway)
#' str(genes.by.celltype.pathway)
NULL



#' genes.by.gobp.pathway
#'
#' @description A go term bp  gene set from \code{
#' x <- readLines("c5.go.bp.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' genes.by.gobp.pathway <- lapply(res, "[", -c(1:2))
#' }
#'
#' @docType data
#' @name genes.by.gobp.pathway
#' @usage data(genes.by.gobp.pathway)
#' @format list
#' @keywords datasets
#' @examples data(genes.by.gobp.pathway)
#' str(genes.by.gobp.pathway)
NULL


#' genes.by.hallmark.pathway
#'
#' @description A hallmark pathway gene set from \code{
#' x <- readLines("h.all.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' genes.by.hallmark.pathway <- lapply(res, "[", -c(1:2))
#' }
#'
#' @docType data
#' @name genes.by.hallmark.pathway
#' @usage data(genes.by.hallmark.pathway)
#' @format list
#' @keywords datasets
#' @examples data(genes.by.hallmark.pathway)
#' str(genes.by.hallmark.pathway)
NULL


#' genes.by.reactome.pathway
#'
#' @description A celltype pathway gene set from \code{
#' x <- readLines("c2.cp.reactome.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' genes.by.reactome.pathway <- lapply(res, "[", -c(1:2))
#' save(genes.by.reactome.pathway,file = "genes.by.reactome.pathway.RData")
#' }
#'
#' @docType data
#' @name genes.by.reactome.pathway
#' @usage data(genes.by.reactome.pathway)
#' @format list
#' @keywords datasets
#' @examples data(genes.by.reactome.pathway)
#' str(genes.by.reactome.pathway)
NULL

#' genes.by.regulatory.pathway
#'
#' @description A celltype pathway gene set from \code{
#' x <- readLines("c3.all.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' genes.by.regulatory.pathway <- lapply(res, "[", -c(1:2))
#' }
#'
#' @docType data
#' @name genes.by.regulatory.pathway
#' @usage data(genes.by.regulatory.pathway)
#' @format list
#' @keywords datasets
#' @examples data(genes.by.regulatory.pathway)
#' str(genes.by.regulatory.pathway)
NULL


#' genes.by.tft.pathway
#'
#' @description A celltype pathway gene set from \code{
#' x <- readLines("c3.tft.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' genes.by.tft.pathway <- lapply(res, "[", -c(1:2))
#' }
#'
#' @docType data
#' @name genes.by.tft.pathway
#' @usage data(genes.by.tft.pathway)
#' @format list
#' @keywords datasets
#' @examples data(genes.by.tft.pathway)
#' str(genes.by.tft.pathway)
NULL

#' block_annotation
#'
#' @description GRCh38 gencode.v44.annotation.gff3.gz
#'
#' @docType data
#' @name block_annotation
#' @usage data(block_annotation)
#' @format a data.frame,
#' @source Generated from gencode.v44.annotation.gff3.gz
#' @keywords datasets
#' @examples data(block_annotation)
#' str(block_annotation)
NULL

#' block_annotation_hg37
#'
#' @description GRCh37 gencode.v44lift37.annotation.gff3
#'
#' @docType data
#' @name block_annotation_hg37
#' @usage data(block_annotation_hg37)
#' @format a data.frame,
#' @source Generated from gencode.v44lift37.annotation.gff3
#' @keywords datasets
#' @examples data(block_annotation_hg37)
#' str(block_annotation_hg37)
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
