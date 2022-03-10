#' @import ggsci
#' @import dplyr
#' @import foreach
#' @import data.table
#' @importFrom irlba irlba
#' @import Seurat
#' @import SeuratObject
#' @importFrom Matrix Matrix
#' @import stringr
#' @import future
#' @import future.apply
#' @import glmnet
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @import utils
#' @import ggplot2
#' @import ggthemes
#' @import ggpubr
#' @import ggtext
#' @import ggnewscale
#' @import bigstatsr
#' @importFrom RMTstat qWishartMax
#' @importFrom gridExtra grid.arrange



#' @title Main wrapper functions for celltypes.
#' @name Pagwas_main
#' @description Main Pagwas wrapper function in order.
#' @details The entry point for Pagwas analysis.
#'
#' @param Pagwas (list)Pagwas data list, default is "NULL"
#' @param gwas_data (data.frame)GWAS summary data
#' @param add_eqtls (character)There are three options: "OnlyTSS","OnlyEqtls","Both"; "OnlyTSS" means there is only snps within TSS window;"OnlyEqtls" represents only significant snps in eqtls files;  "Both" represents both significant snps in eqtls and the other snps in GWAS summary files;
#' @param marg (integr)gene-TSS-window size
#' @param eqtls_files (character)eqtls files names,frome GTEx.
#' @param eqtls_cols (character) the needed columns in eqtls_files; In our examples files, the columns including: c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name")
#' @param block_annotation (data.frame) Start and end points for block traits, usually genes.
#' @param Single_data (Seruat)Input the Single data in seruat format, Idents should be the celltypes annotation.
#' @param nfeatures (integr) The parameter for FindVariableFeatures, NULL means select all genes
#' @param Pathway_list (list,character) pathway gene sets list
#' @param chrom_ld (list,numeric)LD data for 22 chromosome.
#' @param maf_filter (numeric)Filter the maf, default is 0.01
#' @param min.reads (integr)Threshold for total reads fo each cells. default is 5
#' @param min.detected (integr)Threshold for total cells fo each gene. default is 1
#' @param min.lib.size (integr)Threshold for single data library. default is 200
#' @param min_clustercells (integr)Threshold for total cells fo each cluster.default is 10
#' @param min.pathway.size (integr)Threshold for min pathway gene size. default is 5
#' @param max.pathway.size (integr)Threshold for max pathway gene size. default is 300
#' @param iters (integr)number of bootstrap iterations to perform, default is 200
#' @param perform_cv (logical)whether to Perform regularization inference.
#' @param n_folds (integr)folds for regularized inference,default is 10.
#' @param n.cores (integr)Parallel cores,default is 1. use detectCores() to check the cores in computer.
#' @param simp_results (logical)Whether to return a simple results.If TURE,there is only finaly result return; if FALSE, all the Intermediate Data will keep in pagwas, this option is used to save time when rerun the function.
#' @param ... other parameters for subfunction
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(Genes_by_pathway.kegg)
#' data(GWAS_summ_example)
#' data(gtf_df)
#' data(scRNAexample)
#' data(eqtls_files)
#' data(chrom_ld)
#' #1. OnlyTSS
#' Pagwas<-Pagwas_main(Pagwas = NULL,
#'                     gwas_data = GWAS_summ_example,
#'                     add_eqtls="OnlyTSS",
#'                     block_annotation = gtf_df,
#'                     Single_data = scRNAexample,
#'                     Pathway_list=Genes_by_pathway.kegg,
#'                     chrom_ld = chrom_ld)
#' #1. Both
#' scPagwas_main <- function(Pagwas = NULL,
#'                           gwas_data = GWAS_summ_example,
#'                           add_eqtls="Both",
#'                           eqtls_files=eqtls_files,
#'                           eqtl_p=0.05,
#'                           block_annotation = gtf_df,
#'                           Single_data = scRNAexample,
#'                           Pathway_list=Genes_by_pathway.kegg,
#'                           chrom_ld = chrom_ld)
#'
Pagwas_main <- function(Pagwas = NULL,
                        gwas_data = NULL,
                        add_eqtls,
                        eqtls_files=NULL,
                        eqtls_cols=c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name","pval_beta"),
                        block_annotation = NULL,
                        Single_data = NULL,
                        nfeatures =NULL,
                        Pathway_list=NULL,
                        chrom_ld=NULL,
                        marg=10000,
                        maf_filter = 0.01,
                        min.reads = 5,
                        min.detected = 1,
                        min.lib.size = 200,
                        min_clustercells=10,
                        min.pathway.size=5,
                        max.pathway.size=300,
                        iters = 200,
                        perform_cv = F,
                        n_folds = 10,
                        n.cores=1,
                        simp_results=F,
                           ...) {
  if (is.null(Pagwas)) {
  Pagwas <- list();
  class(Pagwas) <- 'Pagwas'
  }
  #1.gwas summary data input

  message(paste(utils::timestamp(quiet = T), ' ******* 1st: GWAS_summary_input function start! ********',sep = ''))

  if (!is.null(gwas_data)){

    Pagwas <- GWAS_summary_input(Pagwas=Pagwas,
                                 gwas_data=gwas_data,
                                 maf_filter=maf_filter)  }
  message('done!')
  #2.single data input
  message(paste(utils::timestamp(quiet = T), ' ******* 2nd: Single_data_input function start! ********',sep = ''))

  #suppressMessages(require(SeuratObject))

  if (!is.null(Single_data)){
    Pagwas <- Single_data_input(Pagwas=Pagwas,
                                nfeatures =nfeatures,
                                Single_data=Single_data,
                                min.lib.size = min.lib.size,
                                min.reads =min.reads,
                                min.detected =min.detected,
                                min_clustercells=min_clustercells)

  }
  message(ncol(Pagwas$Single_data)," cells are remain!" )
  message('done!')

  #3.calculated pca score
  message(paste(utils::timestamp(quiet = T), ' ******* 3rd: Pathway_pcascore_run function start!! ********',sep = ''))


  if (!is.null(Pathway_list)){


  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,n.cores=n.cores,
                                 Pathway_list=Pathway_list,
                                 min.pathway.size=min.pathway.size,
                                 max.pathway.size=max.pathway.size
                                 )

   }
   message('done!')

   #4.calculated Snp2Gene
  message(paste(utils::timestamp(quiet = T), ' ******* 4th: Snp2Gene start!! ********',sep = ''))

  if(!is.null(block_annotation)){

    Pagwas$block_annotation<-block_annotation

    if(add_eqtls!="OnlyTSS"){
      if (!is.null(eqtls_files)){

        message("Filter snps for significant eqtls!")
        Pagwas<-  Tissue_eqtls_Input(Pagwas=Pagwas,
                                     add_eqtls=add_eqtls,
                                     eqtls_files=eqtls_files,
                                     eqtl_p=0.05,
                                     eqtls_cols=eqtls_cols,
                                     marg=marg)
        message(nrow(Pagwas$gwas_data)," snps for significant eqtls! Please set 'add_eqtls=F' if the amount of snps were too little!")

       }else{

        stop("Since the add_eqtls is TURE! No parameter 'eqtls_files' Input!  ")
       }
    }else{
      #Pagwas$block_annotation <- block_annotation
       snp_gene_df<-Snp2Gene(snp=as.data.frame(Pagwas$gwas_data),refGene=block_annotation,marg=marg)
       snp_gene_df$slope <- rep(1,nrow(snp_gene_df))
       Pagwas$snp_gene_df <- snp_gene_df[snp_gene_df$Disstance=="0",]

    }
  }


  #3.pathway block data
  message(paste(utils::timestamp(quiet = T), ' ******* 5th: Pathway_annotation_input function start! ********',sep = ''))

  if (!is.null(block_annotation)){

    Pagwas <- Pathway_annotation_input(Pagwas=Pagwas,n.cores=n.cores)
  }

  message('done!')

  #4.ld data folder,which is preprogress

  message(paste(utils::timestamp(quiet = T), ' ******* 6th: Link_pathway_blocks_gwas function start! ********',sep = ''))

  if (!is.null(chrom_ld)){

    Pagwas <- Link_pathway_blocks_gwas(Pagwas=Pagwas,
                                       chrom_ld=chrom_ld,
                                       n.cores=n.cores)

  message('done!')
    # save(merge_scexpr,file="merge_scexpr.RData")
  #Pathwat_sumexpr<-lapply(Pagwas$Pathway_list,function(x) sapply(merge_scexpr[x,],function(y) sum(y,na.rm =T)))

  message(paste(utils::timestamp(quiet = T), ' ******* 7th: link_pwpca_block function start!! ********',sep = ''))

  #timestart<-Sys.time()
  Pagwas <- link_pwpca_block(Pagwas)
  message('done!')
  }


  message(paste(utils::timestamp(quiet = T), ' ******* 8th: Pagwas_perform_regression function start!! ********',sep = ''))

  #timestart<-Sys.time()
  Pagwas <- Pagwas_perform_regression(Pagwas, iters = iters,n.cores=n.cores)

  message('done!')

  #message(paste(utils::timestamp(quiet = T), ' ******* 9th: Pagwas_perform_regression function start!! ********',sep = ''))


  #write.csv(Pagwas$bootstrap_results,file="Pagwas_kegg_rawResult.csv")
  if (perform_cv) {
  Pagwas <- Pagwas_perform_regularized_regression(Pagwas, n_folds = n_folds)
  }

  if(simp_results){

   Pagwas$Pathway_block_info<- NULL
   Pagwas$gwas_data <- NULL
   Pagwas$snp_gene_df <-NULL
   Pagwas$Single_data <- NULL
   Pagwas$pca_cell_df <-NULL
   Pagwas$merge_scexpr <- NULL
   Pagwas$Pathway_ld_gwas_data <- NULL
   Pagwas$block_annotation <- NULL
   Pagwas$chrom_ld<-NULL

  }

  return(Pagwas)
}


#' scPagwas_main
#' @description Main scPagwas wrapper function for single cell data in order.
#'
#' @param Pagwas (list)Pagwas data list, default is "NULL"
#' @param gwas_data (data.frame)GWAS summary data
#' @param add_eqtls (character)There are three options: "OnlyTSS","OnlyEqtls","Both"; "OnlyTSS" means there is only snps within TSS window;"OnlyEqtls" represents only significant snps in eqtls files;  "Both" represents both significant snps in eqtls and the other snps in GWAS summary files;
#' @param eqtls_files (character)eqtls files names,frome GTEx.
#' @param eqtls_cols (character) the needed columns in eqtls_files; In our examples files, the columns including: c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name")
#' @param eqtl_p (numeric)default is 0.05,
#' @param block_annotation (data.frame) Start and end points for block traits, usually genes.
#' @param Single_data (Seruat)Input the Single data in seruat format, Idents should be the celltypes annotation.
#' @param nfeatures (integr) The parameter for FindVariableFeatures, NULL means select all genes
#' @param Pathway_list (list,character) pathway gene sets list
#' @param chrom_ld (list,numeric)LD data for 22 chromosome.
#' @param marg (integr)gene-TSS-window size
#' @param maf_filter (numeric)Filter the maf, default is 0.01
#' @param min.reads (integr)Threshold for total reads fo each cells. default is 5
#' @param min.detected (integr)Threshold for total cells fo each gene. default is 1
#' @param min.lib.size (integr)Threshold for single data library. default is 200
#' @param min_clustercells (integr)Threshold for total cells fo each cluster.default is 10
#' @param min.pathway.size (integr)Threshold for min pathway gene size. default is 5
#' @param max.pathway.size (integr)Threshold for max pathway gene size. default is 300
#' @param perform_cv (logical)whether to Perform regularization inference.
#' @param n_folds (integr)folds for regularized inference,default is 10.
#' @param n.cores (integr)Parallel cores,default is 1. use detectCores() to check the cores in computer.
#' @param regression (logical)default FALSE, whither to run scPagwas_perform_regression, get the pvalue, the step is too slow.
#' @param simp_results (logical)Whether to return a simple results.If TURE,there is only finaly result return; if FALSE, all the Intermediate Data will keep in pagwas, this option is used to save time when rerun the function.
#' @param split_n (integr)default 1e9, When the cell number is too big, there may have memory errors, set split_n=10000 or other number can split the cell data.
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(Genes_by_pathway.kegg)
#' data(GWAS_summ_example)
#' data(gtf_df)
#' data(scRNAexample)
#' data(eqtls_files)
#' data(chrom_ld)
#' #1. OnlyTSS
#' scPagwas_main <- function(Pagwas = NULL,
#'                           gwas_data = GWAS_summ_example,
#'                           add_eqtls="OnlyTSS",
#'                           block_annotation = gtf_df,
#'                           Single_data = scRNAexample,
#'                           Pathway_list=Genes_by_pathway.kegg,
#'                           chrom_ld = chrom_ld)
#' #1. Both
#' scPagwas_main <- function(Pagwas = NULL,
#'                           gwas_data = GWAS_summ_example,
#'                           add_eqtls="Both",
#'                           eqtls_files=eqtls_files,
#'                           eqtl_p=0.05,
#'                           block_annotation = gtf_df,
#'                           Single_data = scRNAexample,
#'                           Pathway_list=Genes_by_pathway.kegg,
#'                           chrom_ld = chrom_ld)
#'
#'
scPagwas_main <- function(Pagwas = NULL,
                        gwas_data = NULL,
                        add_eqtls="OnlyTSS",
                        eqtls_files=NULL,
                        eqtls_cols=c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name","pval_beta"),
                        eqtl_p=0.05,
                        block_annotation = NULL,
                        Single_data = NULL,
                        nfeatures =NULL,
                        Pathway_list=NULL,
                        chrom_ld = NULL,
                        marg=10000,
                        maf_filter = 0.01,
                        min.reads = 5,
                        min.detected = 1,
                        min.lib.size = 200,
                        min_clustercells=10,
                        min.pathway.size=5,
                        max.pathway.size=300,
                        perform_cv = F,
                        n_folds = 10,
                        n.cores=1,
                        regression=FALSE,
                        simp_results=F,
                        split_n=1e9) {
  if (is.null(Pagwas)){
    Pagwas <- list();
    class(Pagwas) <- 'Pagwas'
  }
   #1.gwas summary data input

  message(paste(utils::timestamp(quiet = T), ' ******* 1st: GWAS_summary_input function start! ********',sep = ''))

  if (!is.null(gwas_data)){

    Pagwas <- GWAS_summary_input(Pagwas=Pagwas,
                                 gwas_data=gwas_data,
                                 maf_filter=maf_filter)
  }
  message('done!')
  #2.single data input
  message(paste(utils::timestamp(quiet = T), ' ******* 2nd: Single_data_input function start! ********',sep = ''))

  suppressMessages(require(SeuratObject))
  if (!is.null(Single_data)){
    Pagwas <- Single_data_input(Pagwas=Pagwas,
                                nfeatures =nfeatures,
                                Single_data=Single_data,
                                min.lib.size = min.lib.size,
                                min.reads =min.reads,
                                min.detected =min.detected,
                                min_clustercells=min_clustercells)

  }
  message(ncol(Pagwas$Single_data)," cells are remain!" )
  message('done!')

  #3.calculated pca score
  message(paste(utils::timestamp(quiet = T), ' ******* 3rd: Pathway_pcascore_run function start!! ********',sep = ''))


  if (!is.null(Pathway_list)){
    Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,n.cores=n.cores,
                                   Pathway_list=Pathway_list,
                                   min.pathway.size=min.pathway.size,
                                   max.pathway.size=max.pathway.size
                                   )

   }
   message('done!')

   #4.calculated Snp2Gene
  message(paste(utils::timestamp(quiet = T), ' ******* 4th: Snp2Gene start!! ********',sep = ''))

  if(!is.null(block_annotation)){

    Pagwas$block_annotation<-block_annotation

    if(add_eqtls!="OnlyTSS"){
      if (!is.null(eqtls_files)){

        message("Filter snps for significant eqtls!")
        Pagwas<-  Tissue_eqtls_Input(Pagwas=Pagwas,
                                     add_eqtls=add_eqtls,
                                     eqtls_files=eqtls_files,
                                     eqtl_p=eqtl_p,
                                     eqtls_cols=eqtls_cols,
                                     marg=marg)
        #add the eqtls genes' snp
        message(nrow(Pagwas$gwas_data)," snps for significant eqtls! Please set 'add_eqtls=F' if the amount of snps were too little!")

       }else{

        stop("Since the add_eqtls is TURE! No parameter 'eqtls_files' Input!  ")
       }
    }else{
      #Pagwas$block_annotation <- block_annotation
       snp_gene_df<-Snp2Gene(snp=Pagwas$gwas_data,refGene=block_annotation,marg=marg)
       snp_gene_df$slope <- rep(1,nrow(snp_gene_df))
       Pagwas$snp_gene_df <- snp_gene_df[snp_gene_df$Disstance=="0",]

    }
  }

  #3.pathway block data
  message(paste(utils::timestamp(quiet = T), ' ******* 5th: Pathway_annotation_input function start! ********',sep = ''))

  if (!is.null(block_annotation)){

    Pagwas <- Pathway_annotation_input(Pagwas=Pagwas,n.cores=n.cores)
  }

  message('done!')

  #4.ld data folder,which is preprogress

  message(paste(utils::timestamp(quiet = T), ' ******* 6th: Link_pathway_blocks_gwas function start! ********',sep = ''))

  if (!is.null(ld_folder)){

    Pagwas <- Link_pathway_blocks_gwas(Pagwas=Pagwas,
                                       chrom_ld=chrom_ld,
                                       n.cores=n.cores)

  message('done!')
  }

  message(paste(utils::timestamp(quiet = T), ' ******* 7th: link_scCell_pwpca_block function start!! ********',sep = ''))
  Pagwas <- link_scCell_pwpca_block(Pagwas,n.cores=n.cores)
  message('done!')

  message(paste(utils::timestamp(quiet = T), ' ******* 8th: scPagwas_perform_score function start!! ********',sep = ''))
  Pagwas <- scPagwas_perform_score(Pagwas, n.cores=n.cores,split_n=split_n)
  message('done!')

  if(regression){
  message(paste(utils::timestamp(quiet = T), ' ******* 9th:scPagwas_perform_regression function start!! ********',sep = ''))
  Pagwas <- scPagwas_perform_regression(Pagwas, n.cores=n.cores)
  message('done!')
    #write.csv(Pagwas$bootstrap_results,file="Pagwas_kegg_rawResult.csv")
    if (perform_cv) {
     Pagwas <- scPagwas_perform_regularized_inference(Pagwas, n_folds = n_folds)
    }

  }

  if(simp_results){
  Pagwas$gwas_data <- NULL
  Pagwas$Single_data <- NULL
  #Pagwas$Celltype_anno <-NULL
  Pagwas$snp_gene_df <-NULL
  Pagwas$chrom_ld<-NULL
  #Pagwas$rawPathway_list <-NULL
  Pagwas$pca_cell_df <-NULL
  Pagwas$block_annotation <- NULL
  Pagwas$merge_scexpr <- NULL
  Pagwas$pathway_blocks <- NULL
  Pagwas$Pathway_ld_gwas_data <- NULL

  }

  return(Pagwas)
}
