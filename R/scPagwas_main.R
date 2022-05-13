
#' @importFrom dplyr mutate filter inner_join %>%
#' @importFrom irlba irlba
#' @importFrom Seurat FindVariableFeatures AverageExpression VariableFeatures GetAssayData RunPCA RunTSNE RunUMAP Embeddings
#' @importFrom SeuratObject Idents
#' @importFrom Matrix Matrix colSums rowSums crossprod
#' @importFrom glmnet cv.glmnet
#' @importFrom GenomicRanges GRanges resize resize
#' @importFrom IRanges IRanges
#' @importFrom utils timestamp
#' @import ggplot2
#' @importFrom ggpubr ggscatter
#' @importFrom bigstatsr as_FBM big_apply big_univLinReg covar_from_df big_transpose big_cprodMat
#' @importFrom RMTstat qWishartMax
#' @importFrom gridExtra grid.arrange
#' @importFrom data.table setkey data.table as.data.table
#' @importFrom bigreadr fread2 rbind_df
#' @importFrom reshape2 dcast
#' @importFrom bigmemory as.big.matrix is.big.matrix
#' @importFrom biganalytics apply

#' @title Main wrapper functions for scPagwas
#' @name scPagwas_main
#' @description Main Pagwas wrapper re-progress function.
#' @details The entry point for Pagwas analysis.
#'
#' @param Pagwas (list)Pagwas data list, default is "NULL"
#' @param gwas_data (data.frame)GWAS summary data
#' @param output.prefix output.prefix,prefix used for output tables
#' @param add_eqtls (character)There are three options: "OnlyTSS","OnlyEqtls","Both"; "OnlyTSS" means there is only snps within TSS window;"OnlyEqtls" represents only significant snps in eqtls files;  "Both" represents both significant snps in eqtls and the other snps in GWAS summary files;
#' @param marg (integr)gene-TSS-window size
#' @param eqtls_files (character)eqtls files names,frome GTEx.
#' @param eqtls_cols (character) the needed columns in eqtls_files; In our examples files, the columns including: c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name")
#' @param block_annotation (data.frame) Start and end points for block traits, usually genes.
#' @param Single_data (Seruat)Input the Single data in seruat format, Idents should be the celltypes annotation.
#' @param assay (character)assay data of your single cell data to use,default is "RNA"
#' @param Pathway_list (list,character) pathway gene sets list
#' @param chrom_ld (list,numeric)LD data for 22 chromosome.
#' @param split_n number of times to compute the singlecell result
#' @param maf_filter (numeric)Filter the maf, default is 0.01
#' @param min_clustercells (integr)Only use is when FilterSingleCell is TRUE.Threshold for total cells fo each cluster.default is 10
#' @param min.pathway.size (integr)Threshold for min pathway gene size. default is 5
#' @param max.pathway.size (integr)Threshold for max pathway gene size. default is 300
#' @param iters (integr)number of bootstrap iterations to perform
#' @param remove_outlier (logical)Whether to remove the outlier for scPagwas score.
#' @param param.file (logical)whether save parameters used for scPagwas.
#' @param SimpleResult (logical)whether simplify the scPagwas result.
#' @param log.file (character)log.file
#' @param ncores (integr)Parallel cores,default is 1. use detectCores() to check the cores in computer.
#' @param singlecell
#' @param celltype
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(Genes_by_pathway_kegg)
#' data(GWAS_summ_example)
#' data(gtf_df)
#' data(scRNAexample)
#' #data(eqtls_files)
#' data(chrom_ld)
#' #1. OnlyTSS
#' Pagwas<-scPagwas_main(Pagwas = NULL,
#'                     gwas_data = GWAS_summ_example,
#'                     add_eqtls="OnlyTSS",
#'                     block_annotation = gtf_df,
#'                     Single_data = scRNAexample,
#'                     Pathway_list=Genes_by_pathway_kegg,
#'                     chrom_ld = chrom_ld)
#'
scPagwas_main <- function(Pagwas = NULL,
                        gwas_data = NULL,
                        output.prefix,
                        add_eqtls=c("OnlyTSS","OnlyEqtls","Both"),
                        eqtls_files=NULL,
                        eqtls_cols=c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name","pval_beta"),
                        block_annotation = NULL,
                        Single_data = NULL,
                        assay=c("RNA","SCT"),
                        #Store_CACHE="test",
                        Pathway_list=NULL,
                        chrom_ld=NULL,
                        split_n=1,
                        marg=10000,
                        maf_filter = 0.01,
                        min_clustercells=10,
                        min.pathway.size=5,
                        max.pathway.size=300,
                        iters=200,
                        param.file=T,
                        singlecell=T,
                        celltype=T,
                        remove_outlier=T,
                        SimpleResult=F,
                        log.file='scPagwas.run.log',
                        ncores=1) {
  #debug
  # Pagwas = NULL;
  # gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas");
  # output.prefix="test";
  # add_eqtls="OnlyTSS";
  # eqtls_files=NULL;
  # eqtls_cols=c("rs_id_dbSNP151_GRCh38p7","variant_pos","tss_distance","gene_chr", "gene_start", "gene_end","gene_name","pval_beta");
  # block_annotation = block_annotation;
  # Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas");
  # assay="RNA";
  # Pathway_list=Genes_by_pathway_kegg;
  # chrom_ld=chrom_ld;
  # split_n=3;
  # marg=10000;
  # singlecell=T;
  # celltype=T;
  # maf_filter = 0.01;
  # min_clustercells=10;
  # min.pathway.size=5;
  # max.pathway.size=300;
  # iters=200;
  # param.file=T;
  # remove_outlier=T;
  # SimpleResult=F;
  # log.file='scPagwas.run.log';
  # ncores=2;

  #######
  ## initialize log-file
  cat('##', format(Sys.time()), '\n', file=log.file)

  ## miximal file path lenght;
  ## Windows OS support max. 259 characters
  max.nchar.file.path <- 259

  ## arguments
  #add_eqtls <- match.arg( add_eqtls)
  #assay <- match.arg(assay)

  if(param.file){
    ## save parameters used for ssGSEA
    param.str = c(
      paste('##', Sys.time()),
      paste('input gwas data: ', gwas_data, sep='\t'),
      paste('add_eqtls: ', add_eqtls, sep='\t'),
      paste('eqtls_files: ', eqtls_files, sep='\t'),
      paste('Single_data: ', Single_data, sep='\t'),
      paste('assay: ', assay, sep='\t'),
      #paste('nfeatures: ', nfeatures, sep='\t'),
      paste('Pathway length: ', length(Pathway_list),collapse = " ", sep='\t'),
      paste('split_n: ', split_n, sep='\t'),
      paste('marg: ', marg, sep='\t'),
      paste('maf_filter: ', maf_filter, sep='\t'),
      paste('min_clustercells: ', min_clustercells, sep='\t'),
      paste('min.pathway.size: ', min.pathway.size, sep='\t'),
      paste('max.pathway.size: ', max.pathway.size, sep='\t'),
      paste('remove_outlier: ', remove_outlier, sep='\t'),
      paste('iters: ', iters, sep='\t'),
      paste('ncores: ', ncores, sep='\t')
    )
    writeLines(param.str, con=paste(output.prefix, 'parameters.txt', sep='_'))
  }

  #Sys.setenv(R_LOCAL_CACHE=Store_CACHE)

  tt <- Sys.time()
  if (is.null(Pagwas)) {
  Pagwas <- list();
  class(Pagwas) <- 'Pagwas'
  }
   message(paste(utils::timestamp(quiet = T), ' ******* 1st: Single_data_input function start! ********',sep = ''))

  #suppressMessages(require(SeuratObject))
  if(!is.null(Single_data)){
  if(grepl(".rds",Single_data)){
    message("** Start to read the single cell data!")
    Single_data=readRDS(Single_data)
  }else{
    stop("There is need a data in `.rds` format ")
  }

  if(!assay %in% Assays(Single_data)){
    stop("There is no need assays in your Single_data")
  }
    Pagwas <- Single_data_input(Pagwas=Pagwas,
                                assay=assay,
                                #nfeatures =nfeatures,
                                Single_data=Single_data,
                                Pathway_list=Pathway_list,
                                min_clustercells=min_clustercells)
  rm(Single_data)
  #message(ncol(Pagwas$Single_data)," cells are remain!" )
  message('done!')
  }

  cat('Single_data import: ',  file=log.file, append=T)
  cat(Sys.time()-tt, '\n',  file=log.file, append=T)

  #3.calculated pca score
  message(paste(utils::timestamp(quiet = T), ' ******* 2nd: Pathway_pcascore_run function start!! ********',sep = ''))

  tt <- Sys.time()
  if (!is.null(Pathway_list)){


  Pagwas <- Pathway_pcascore_run(Pagwas=Pagwas,
                                 Pathway_list=Pathway_list,
                                 min.pathway.size=min.pathway.size,
                                 max.pathway.size=max.pathway.size
                                 )

   }
   message('done!')
   cat('Pathway_pcascore_run: ',  file=log.file, append=T)
   cat(Sys.time()-tt, '\n',  file=log.file, append=T)


   message(paste(utils::timestamp(quiet = T), ' ******* 3rd: GWAS_summary_input function start! ********',sep = ''))

   tt <- Sys.time()
   if(class(gwas_data)=="character"){
     message("** Start to read the gwas_data!")

     suppressMessages(gwas_data <- bigreadr::fread2(gwas_data))

   }else{
     stop("There is need a filename and address for gwas_data")
   }

     Pagwas <- GWAS_summary_input(Pagwas=Pagwas,
                                  gwas_data=gwas_data,
                                  maf_filter=maf_filter)
   rm(gwas_data)
   message('done!')
   cat('GWAS_summary_input: ',  file=log.file, append=T)
   cat(Sys.time()-tt, '\n',  file=log.file, append=T)

   #4.calculated Snp2Gene
  message(paste(utils::timestamp(quiet = T), ' ******* 4th: Snp2Gene start!! ********',sep = ''))

  tt <- Sys.time()
  if(!is.null(block_annotation)){
    if(add_eqtls!="OnlyTSS"){
      if (!is.null(eqtls_files)){
        message("Filter snps for significant eqtls!")
        Pagwas<-  Tissue_eqtls_Input(Pagwas=Pagwas,
                                     block_annotation=block_annotation,
                                     add_eqtls=add_eqtls,
                                     eqtls_files=eqtls_files,
                                     eqtl_p=0.05,
                                     eqtls_cols=eqtls_cols,
                                     marg=marg)

       }else{

        stop("Since the add_eqtls is TURE! No parameter 'eqtls_files' Input!  ")
       }
    }else{

       snp_gene_df<-Snp2Gene(snp=Pagwas$gwas_data,refGene=block_annotation,marg=marg)
       snp_gene_df$slope <- rep(1,nrow(snp_gene_df))
       Pagwas$snp_gene_df<-snp_gene_df[snp_gene_df$Disstance=="0",]
    }
  }
  cat('Snp2Gene: ',  file=log.file, append=T)
  cat(Sys.time()-tt, '\n',  file=log.file, append=T)

  #3.pathway block data
  message(paste(utils::timestamp(quiet = T), ' ******* 5th: Pathway_annotation_input function start! ********',sep = ''))

  tt <- Sys.time()
  if (!is.null(block_annotation)){
    Pagwas <- Pathway_annotation_input(Pagwas=Pagwas,
                                       block_annotation=block_annotation)
  }

  message('done!')
  cat('Pathway_annotation_input: ',  file=log.file, append=T)
  cat(Sys.time()-tt, '\n',  file=log.file, append=T)

  #4.ld data folder,which is preprogress

  message(paste(utils::timestamp(quiet = T), ' ******* 6th: Link_pathway_blocks_gwas function start! ********',sep = ''))

  tt <- Sys.time()
  if (!is.null(chrom_ld)){

  Pagwas <- Link_pathway_blocks_gwas(Pagwas=Pagwas,
                                     chrom_ld=chrom_ld,
                                     split_n=split_n,
                                     singlecell=singlecell,
                                     celltype=celltype,
                                     ncores=ncores)

   message('done!')
  }
  cat('Link_pathway_blocks_gwas: ',  file=log.file, append=T)
  cat(Sys.time()-tt, '\n',  file=log.file, append=T)

  if(celltype){

  message(paste(utils::timestamp(quiet = T), ' ******* 7th: Celltype_heritability_contributions function start! ********',sep = ''))

  Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data=Pagwas$Pathway_ld_gwas_data)
  Pagwas <- Boot_evaluate(Pagwas,bootstrap_iters = iters, part = 0.5)

  Pagwas$Pathway_ld_gwas_data<-NULL

  message('done!')
  cat('Celltype_heritability_contributions: ',  file=log.file, append=T)
  cat(Sys.time()-tt, '\n',  file=log.file, append=T)

  }

  if(singlecell){
  message(paste(utils::timestamp(quiet = T), ' ******* 8th: scPagwas_perform_score function start! ********',sep = ''))
  tt <- Sys.time()
  Pagwas<-scPagwas_perform_score(Pagwas=Pagwas,
                                 remove_outlier=TRUE)
  message('done!')
  cat('scPagwas_perform_score: ',  file=log.file, append=T)
  cat(Sys.time()-tt, '\n',  file=log.file, append=T)

  message(paste(utils::timestamp(quiet = T), ' ******* 9th: scGet_gene_heritability_correlation function start! ********',sep = ''))

  message("** Get gene heritability contributions!")
  Pagwas <- scGet_gene_heritability_correlation(Pagwas=Pagwas)
  message("done")
  }
  if(SimpleResult){
    Pagwas[c("VariableFeatures","merge_scexpr",
             "data_mat","rawPathway_list",
             "snp_gene_df")]<-NULL
  }
  gc()
  return(Pagwas)
}
