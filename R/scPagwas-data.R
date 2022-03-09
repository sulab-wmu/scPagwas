
#' KEGG pathway gene list
#'
#' A few numbers from \code{
#' library(KEGGREST)
#' pathways.list <- keggList("pathway", "hsa")
#' # Pull all genes for each pathway
#' pathway.codes <- sub("path:", "", names(pathways.list))
#' Genes_by_pathway.kegg <- sapply(pathway.codes,function(pwid){
#' pw <- keggGet(pwid)
#' if (is.null(pw[[1]]$GENE)) return(NA)
#' pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
#' pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
#'  return(pw2)
#' })
#'
#' }
#'
#' @docType data
#' @name Genes_by_pathway.kegg
#' @format list.
#' @source Generated from KEGGREST
#' @examples data(Genes_by_pathway.kegg)
#' str(Genes_by_pathway.kegg)
NULL



#' Genes_by_pathway.hallmark
#' A hallmark pathway gene set from \code{
#' x <- readLines("E:/OneDrive/GWAS_Multiomics/tempdata/h.all.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' Genes_by_pathway.hallmark <- lapply(res, "[", -c(1:2))
#' save(Genes_by_pathway.hallmark,file="E:/RPakage/scPagwas/data/Genes_by_pathway.hallmark.RData")
#' }
#' @docType data
#' @name Genes_by_pathway.hallmark
#' @format list
#' @source Generated from h.all.v7.5.1.symbols.gmt
#' @examples data(Genes_by_pathway.hallmark)
#' str(Genes_by_pathway.hallmark)
NULL

#' Genes_by_pathway.reactome
#'
#' A reactome pathway gene set from \code{
#' x <- readLines("E:/OneDrive/GWAS_Multiomics/tempdata/c2.cp.reactome.v7.5.1.symbols.gmt")
#' res <- strsplit(x, "\t")
#' names(res) <- vapply(res, function(y) y[1], character(1))
#' Genes_by_pathway.reactome <- lapply(res, "[", -c(1:2))
#' save(Genes_by_pathway.reactome,file="E:/RPakage/scPagwas/data/Genes_by_pathway.reactome.RData")
#' }
#'
#' @docType data
#' @name Genes_by_pathway.reactome
#' @format list
#' @source Generated from c2.cp.reactome.v7.5.1.symbols.gmt
#' @examples data(Genes_by_pathway.reactome)
#' str(Genes_by_pathway.reactome)
NULL




#' GWAS summary data frame
#'
#' A GWAS summary data from \code{
#' library(readr)
#' library(dplyr)
#' GWAS_liver<-read_delim("/share/pub/mayl/Sherlock/XiangbingYu_Primary_biliary_cholangitis/PBC_GWAS_UKBiobank_summary_final",delim="\t")
#' #select the specific coloumn names: chrom，pos，rsid，beta，se，maf
#' GWAS_liver<-GWAS_liver[,c("CHR","POS","SNP","p","SE")]
#' colnames(GWAS_liver)<-c("chrom","pos","rsid","p","se")
#' COVOD19<-read_delim("/share/pub/mayl/01_COVID19_GWAS_round4/COVID19_HGI_B2_ALL_leave_23andme_20201020.b37.txt",delim="\t")
#' COVOD19<-COVOD19[,c("all_meta_AF","rsid")]
#' colnames(COVOD19)<-c("maf","rsid")
#' ##get beta
#' GWAS_liver$beta<-GWAS_liver$se*qchisq(GWAS_liver$p,1,lower.tail=F)
#' GWAS_liver_inner <- inner_join(GWAS_liver,COVOD19,by="rsid")
#' GWAS_liver_inner <- as.data.frame(GWAS_liver_inner)
#' GWAS_summ_example <- as.data.frame(GWAS_liver_inner[sample(nrow(GWAS_liver_inner),10000),])
#' save(GWAS_summ_example,file="/share/pub/dengcy/GWAS_Multiomics/pagwas/data/GWAS_summ_example.RData")
#' }
#'
#' @docType data
#' @name GWAS_summ_example
#' @format data.frame.
#' @source Generated from PBC_GWAS_UKBiobank_summary_final
#' @examples data(GWAS_summ_example)
#' str(GWAS_summ_example)
NULL


#' gtf_df
#'
#' A gene annotation files from \code{
#' library("rtracklayer")
#' gtf_df<- rtracklayer::import("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/gencode.v34.annotation.gtf.gz")
#' gtf_df <- as.data.frame(gtf)
#' gtf_df <- gtf_df[,c("seqnames","start","end","type","gene_name")]
#' gtf_df <- gtf_df[gtf_df$type=="gene",]
#' gtf_df<-gtf_df[,c(1,2,3,5)]
#' colnames(gtf_df)<-c("chrom", "start","end","label")
#' save(gtf_df,file = "/share/pub/dengcy/GWAS_Multiomics/pagwas/data/gtf_df.RData")
#' }
#'
#' @docType data
#' @name gtf_df
#' @format data.frame.
#' @source Generated from gencode.v34.annotation.gtf.gz
#' @examples data(gtf_df)
#' str(gtf_df)
NULL


#' scRNAexample
#' A single cell example data from \code{
#' liver_data<-readRDS("/share/pub/qiuf/COVID/liver/Rpoly.rds")
#' #set Idents
#' liver_data@meta.data$annotation <- liver_data@meta.data$anno
#' Idents(liver_data) <- liver_data@meta.data$anno
#' #normalization
#' liver_data <- NormalizeData(liver_data, normalization.method = "LogNormalize", scale.factor = 10000)
#' #correct the cell names, need no specific signal
#' anno <- liver_data@meta.data$annotation
#'  anno<-str_replace_all(anno,"-",".")
#'  anno<-str_replace_all(anno," ",".")
#' anno[which(anno=="Inflammatory.monocyte/macrophages")]<-"Inflammatory.monocyte.macrophages"
#' Celltype_anno$annotation<-str_replace_all(anno,"\\+",".")
#' Idents(liver_data) <- as.factor(anno)
#' liver_data<-subset(liver_data,idents = c("Stellate.cell","unknow","γδ.T.cell"), invert = TRUE)
#' liver_data2<-liver_data[,sample(ncol(liver_data),1000)]
#' save(scRNAexample,file="/share/pub/dengcy/GWAS_Multiomics/pagwas/data/scRNAexample.RData")
#' }
#'
#' @docType data
#' @name scRNAexample
#' @format Seruat
#' @source Generated from PBC
#' @examples data(scRNAexample)
#' str(scRNAexample)
NULL



#' eqtls_files
#' Generated from GTEx
#' @docType data
#' @name eqtls_files
#' @format vector
#' @source Generated from GTEx
#' @examples data(eqtls_files)
#' str(eqtls_files)
eqtls_files<-"./inst/extdata/Liver.v8.egenes.txt.gz"


#' ld_folder
#' use vcftools to get ped and map files
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
#' @format list
#' @source Generated from PLINK 1.90 linux
#' @examples data(chrom_ld)
#' str(chrom_ld)
NULL

