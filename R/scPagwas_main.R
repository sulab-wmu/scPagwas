#' scPagwas, A single-cell pathway-based principal component (PC)-scoring algorithm named
#' for integrating polygenic signals from GWAS with single cell data
#' to infer disease-relevant cell populations and score the associations of
#' individual cells with complex diseases.
#'
#'
#' @importFrom dplyr mutate filter inner_join %>%
#' @importFrom irlba irlba
#' @importFrom Seurat FindVariableFeatures AverageExpression VariableFeatures GetAssayData RunPCA RunTSNE RunUMAP Embeddings CreateAssayObject DefaultAssay AddModuleScore
#' @importFrom SeuratObject Idents DefaultAssay GetAssayData SetAssayData
#' @importFrom Matrix Matrix colSums rowSums crossprod
#' @importFrom glmnet cv.glmnet
#' @importFrom GenomicRanges GRanges resize resize
#' @importFrom IRanges IRanges
#' @importFrom utils timestamp
#' @import ggplot2
#' @importFrom SOAR Store Objects Remove
#' @importFrom ggpubr ggscatter
#' @importFrom bigstatsr as_FBM big_apply big_univLinReg covar_from_df big_transpose big_cprodMat
#' @importFrom gridExtra grid.arrange
#' @importFrom data.table setkey data.table as.data.table
#' @importFrom bigreadr fread2 rbind_df
#' @importFrom reshape2 dcast
#' @importFrom bigmemory as.big.matrix is.big.matrix
#' @importFrom biganalytics apply

#' @title Main wrapper functions for scPagwas
#' @name scPagwas_main
#' @description Main Pagwas wrapper for workflow function.
#' @details The entry point for Pagwas analysis.
#'
#' @param Pagwas (list)default is "NULL" when you first run the funciton; Pagwas
#' should be list class; Sometimes, It can inherit the result from "scPagwas_main" function
#' last time, when you turn the "seruat_return" to FALSE; It is suitable for circulation
#' running for the same single data.
#' @param gwas_data (data.frame)GWAS summary data; It must have some colmuns such as:
#'  chrom|    pos    |   rsid    |   se  |  beta |  maf
#'     6 | 119968580 | rs1159767 | 0.032 | 0.019 |0.5275
#'    10 | 130566523 |  rs559109 | 0.033 | 0.045 |0.4047
#'     5 | 133328825 | rs6893145 | 0.048 | 0.144 |0.1222
#'     7 | 146652932 | rs13228798| 0.035 | 0.003 | 0.3211
#' @param marg (integr) the distance to TSS site,default is 10000, then gene-TSS-window size is 20000.
#' @param block_annotation (data.frame) Start and end position for block traits, usually genes.
#' @param Single_data (character or Seruat)Input the Single data in seruat format,
#' or the Seruat data address for rds format;the "Idents" should be the celltypes
#' annotation, if the "celltypes" is TRUE.
#' @param assay (character)assay data of your single cell data to use,"RNA" or "SCT",default is "RNA";
#' @param Pathway_list (list,character) pathway gene sets list,the name of list is pahtway names;
#' the element of list including genes;
#' @param chrom_ld (list,numeric)LD data for 22 chromosomes. Processed data is \code{data(chrom_ld)}
#' @param split_n (integr)Number of blocks which are split from singlecell data;
#' sometimes,the data is too big to run; we also advice split the data before run the main function.
#' @param maf_filter (numeric)Filter the "maf" for gwas_data, default is 0.01;
#' @param min_clustercells (integr)Threshold for number of cells in each cluster.default is 10;
#' @param min.pathway.size (integr)Threshold for min pathway gene size. default is 5
#' @param max.pathway.size (integr)Threshold for max pathway gene size. default is 300
#' @param iters (integr)Number of bootstrap iterations to perform, default is 200;
#' @param remove_outlier (logical)Whether to remove the outlier for scPagwas score.
#' @param ncores (integr)Parallel cores,default is 1. use \code{detectCores()} to check the cores in computer.
#' @param seruat_return (logical) Whether return the Seruat format result, if not,will return a list result;
#' @param singlecell (logical)Whether to produce the singlecell result;
#' @param celltype (logical)Whether to produce the celltypes result;
#' @param output.prefix (character)output.prefix,prefix used for output tables;
#' @param output.dirs (character)output directory for result; If the directory is nonexistence,
#' the function will create it;
#' @param n_topgenes (integr)Number of top associated gene selected to calculate the scPagwas score;
#'
#' @return
#' Returns a Seruat data with entries(seruat_return=T):
#' \describe{
#'   \position{assay:}{
#'   \item{scPagwasPaPca:}{An assay for S4 type of data; the svd result for pathways in each cells;}}
#'
#'   \position{meta.data}{
#'   \item{scPagwas.lmtopgenes.Score1:}{ the column for "meta.data";Enrichment socre for inheritance associated top genes.}
#'   \item{sclm_score:}{ the column for "meta.data";Inheritance regression effects for each cells}}
#'
#'   \position{misc: element in result,\code{Pagwas@misc }}{
#'   \item{Pathway_list:}{a list for pathway gene list intersecting with single cell data}
#'   \item{pca_cell_df:}{ a data frame for pathway pca result for each celltype.}
#'   \item{sclm_results:}{ the regression result for each cell.}
#'   \item{allsnp_gene_heritability_correlation:}{
#'   heritability correlation for each gene;}
#'   \item{Pathway_ctlm_results:}{ the regression result for each pathway in each celltype}
#'   \item{lm_results:}{ the regression result for celltypes}
#'   \item{Pathway_ct_results:}{Inheritance conmtribution Matrix for pathways and celltypes, it contribute from "Pathway_ctlm_results"}
#' }
#'
#' }
#' Returns files:
#'
#' \describe{
#'   \item{d:}{ max(nu, nv) approximate singular values}
#'   \item{u:}{ nu approximate left singular vectors (only when right_only=FALSE)}
#'   \item{v:}{ nv approximate right singular vectors}
#'   \item{iter:}{ The number of Lanczos iterations carried out}
#'   \item{mprod:}{ The total number of matrix vector products carried out}
#' }
#'
#' Returns a Seruat data with entries(seruat_return=T):
#' \describe{
#'   \item{scPagwasPaPca:}{Assays for S4 type of data; the svd result for pathways in each cells;}
#'   \item{scPagwas.lmtopgenes.Score1:}{ the column for "meta.data";Enrichment socre for inheritance associated top genes.}
#'   \item{sclm_score:}{ the column for "meta.data";Inheritance regression effects for each cells}
#'   \item{Pathway_list:}{ The number of Lanczos iterations carried out}
#'   \item{pca_cell_df:}{ The total number of matrix vector products carried out}
#'   \item{sclm_results:}{ The total number of matrix vector products carried out}
#'   \item{allsnp_gene_heritability_correlation:}{ The total number of matrix vector products carried out}
#'   \item{Pathway_ctlm_results:}{ The total number of matrix vector products carried out}
#'   \item{lm_results:}{ The total number of matrix vector products carried out}
#'   \item{Pathway_ct_results:}{ The total number of matrix vector products carried out}
#' }
#'
#' @note
#' 1.When you run the package in linux server, you can run
#' \code{export OPENBLAS_NUM_THREADS=1}
#' before enter into R environment;
#'
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(Genes_by_pathway_kegg)
#' data(GWAS_summ_example)
#' data(gtf_df)
#' data(scRNAexample)
#' # data(eqtls_files)
#' data(chrom_ld)
#' # 1. OnlyTSS
#' Pagwas <- scPagwas_main(
#'   Pagwas = NULL,
#'   gwas_data = GWAS_summ_example,
#'   add_eqtls = "OnlyTSS",
#'   block_annotation = gtf_df,
#'   Single_data = scRNAexample,
#'   Pathway_list = Genes_by_pathway_kegg,
#'   chrom_ld = chrom_ld
#' )
scPagwas_main <- function(Pagwas = NULL,
                          gwas_data = NULL,
                          output.prefix = "Test",
                          output.dirs = "scPagwastest_output",
                          block_annotation = NULL,
                          Single_data = NULL,
                          assay = c("RNA", "SCT"),
                          Pathway_list = NULL,
                          chrom_ld = NULL,
                          split_n = 1,
                          marg = 10000,
                          maf_filter = 0.01,
                          min_clustercells = 10,
                          min.pathway.size = 5,
                          max.pathway.size = 300,
                          iters = 200,
                          n_topgenes = 1000,
                          singlecell = T,
                          celltype = T,
                          seruat_return = T,
                          remove_outlier = T,
                          ncores = 1){
  #####################################
  # debug test
  # Pagwas = NULL;
  # gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas");
  # "E:/RPakage/scPagwas/inst/extdata/AD_prune_gwas_data.txt";
  #
  # output.prefix="test";
  # output.dirs="scPagwastest_output";
  # block_annotation = block_annotation;
  # Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas");
  # "E:/RPakage/scPagwas/inst/extdata/GSE138852_ad.rds";
  #
  # assay="RNA";
  # Pathway_list=Genes_by_pathway_kegg;
  # chrom_ld=chrom_ld;
  # split_n=1;
  # marg=10000;
  # singlecell=T;
  # celltype=T;
  # n_topgenes=1000;
  # maf_filter = 0.01;
  # min_clustercells=10;
  # min.pathway.size=5;
  # max.pathway.size=300;
  # iters=200;
  # seruat_return=T;
  # param.file=T;
  # remove_outlier=T;
  # log.file='scPagwas.run.log';
  # ncores=2;

  #######
  ## initialize log-file
  if (!dir.exists(output.dirs)) {
    dir.create(output.dirs)
  }
  cat("## start at:", format(Sys.time()), "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"))
  ## miximal file path lenght;
  ## Windows OS support max. 259 characters
  max.nchar.file.path <- 259

  param.str <- c(
    paste("##", Sys.time()),
    paste("input gwas data: ", gwas_data, sep = "\t"),
    paste0("Single_data: ", if (class(Single_data) == "Seurat") dim(Single_data) else Single_data),
    paste("assay: ", assay, sep = "\t"),
    paste("Pathway length: ", length(Pathway_list), collapse = " ", sep = "\t"),
    paste("split_n: ", split_n, sep = "\t"),
    paste("marg: ", marg, sep = "\t"),
    paste("maf_filter: ", maf_filter, sep = "\t"),
    paste("min_clustercells: ", min_clustercells, sep = "\t"),
    paste("min.pathway.size: ", min.pathway.size, sep = "\t"),
    paste("max.pathway.size: ", max.pathway.size, sep = "\t"),
    paste("singlecell: ", singlecell, sep = "\t"),
    paste("celltype: ", celltype, sep = "\t"),
    paste("n_topgenes: ", n_topgenes, sep = "\t"),
    paste("seruat_return: ", seruat_return, sep = "\t"),
    paste("remove_outlier: ", remove_outlier, sep = "\t"),
    paste("iters: ", iters, sep = "\t"),
    paste("ncores: ", ncores, sep = "\t")
  )
  writeLines(param.str, con = paste0("./", output.dirs, "/", output.prefix, "_parameters.txt"))

  Sys.setenv(R_LOCAL_CACHE = paste0("./", output.dirs, "/scPagwas_cache"))

  tt <- Sys.time()
  if (is.null(Pagwas)) {
    Pagwas <- list()
  } else if (class(Pagwas) != "list") {
    stop("The class of Pagwas for input is wrong! Should be a list data")
  }

  #############################
  ## 1.Single_data_input
  #############################
  message(paste(utils::timestamp(quiet = T), " ******* 1st: Single_data_input function start! ********", sep = ""))

  if (!is.null(Single_data)) {
    if (class(Single_data) == "character") {
      if (grepl(".rds", Single_data)) {
        message("** Start to read the single cell data!")
        Single_data <- readRDS(Single_data)
      } else {
        stop("There is need a data in `.rds` format ")
      }
    } else if (class(Single_data) != "Seurat") {
      stop("There is need a Seurat class for Single_data")
    }

    if (!assay %in% Assays(Single_data)) {
      stop("There is no need assays in your Single_data")
    }

    Pagwas <- Single_data_input(
      Pagwas = Pagwas,
      assay = assay,
      Single_data = Single_data,
      Pathway_list = Pathway_list,
      min_clustercells = min_clustercells
    )

    Single_data <- Single_data[, colnames(Pagwas$data_mat)]
    SOAR::Store(Single_data)
    message("done!")
  } else {
    stop("Single_data must input!")
  }

  cat("Single_data import: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
  cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)

  #############################
  ## 2.Pathway_pcascore_run
  #############################
  message(paste(utils::timestamp(quiet = T), " ******* 2nd: Pathway_pcascore_run function start!! ********", sep = ""))

  tt <- Sys.time()
  if (!is.null(Pathway_list)) {
    Pagwas <- Pathway_pcascore_run(
      Pagwas = Pagwas,
      Pathway_list = Pathway_list,
      min.pathway.size = min.pathway.size,
      max.pathway.size = max.pathway.size
    )
  } else if (!("pca_scCell_mat" %in% names(Pagwas)) ) {
    stop("Pathway_list should input!")
  }
  message("done!")
  cat("Pathway_pcascore_run: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
  cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)

  #############################
  ## 3.GWAS_summary_input
  #############################
  message(paste(utils::timestamp(quiet = T), " ******* 3rd: GWAS_summary_input function start! ********", sep = ""))

  tt <- Sys.time()
  if (class(gwas_data) == "character") {
    message("** Start to read the gwas_data!")
    suppressMessages(gwas_data <- bigreadr::fread2(gwas_data))
  } else {
    stop("There is need a filename and address for gwas_data")
  }

  Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data,
    maf_filter = maf_filter
  )
  rm(gwas_data)
  message("done!")
  cat("GWAS_summary_input: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
  cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)

  #############################
  ## 4.Snp2Gene
  #############################
  message(paste(utils::timestamp(quiet = T), " ******* 4th: Snp2Gene start!! ********", sep = ""))

  tt <- Sys.time()
  if (!is.null(block_annotation)) {
      snp_gene_df <- Snp2Gene(snp = Pagwas$gwas_data, refGene = block_annotation, marg = marg)
      #snp_gene_df$slope <- rep(1, nrow(snp_gene_df))
      Pagwas$snp_gene_df <- snp_gene_df[snp_gene_df$Disstance == "0", ]

  } else if (!("snp_gene_df" %in% names(Pagwas))) {
    stop("block_annotation should input!")
  }
  cat("Snp2Gene: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
  cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)

  #############################
  ## 5.Pathway_annotation_input
  #############################
  message(paste(utils::timestamp(quiet = T), " ******* 5th: Pathway_annotation_input function start! ********", sep = ""))

  tt <- Sys.time()
  if (!is.null(block_annotation)) {
    Pagwas <- Pathway_annotation_input(
      Pagwas = Pagwas,
      block_annotation = block_annotation
    )
  } else if (!("snp_gene_df" %in% names(Pagwas))) {
    stop("block_annotation should input!")
  }

  message("done!")
  cat("Pathway_annotation_input: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
  cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)

  #############################
  ## 6.Link_pathway_blocks_gwas
  #############################
  message(paste(utils::timestamp(quiet = T), " ******* 6th: Link_pathway_blocks_gwas function start! ********", sep = ""))

  tt <- Sys.time()
  if (!is.null(chrom_ld)) {
    Pagwas <- Link_pathway_blocks_gwas(
      Pagwas = Pagwas,
      chrom_ld = chrom_ld,
      split_n = split_n,
      singlecell = singlecell,
      celltype = celltype,
      ncores = ncores
    )

    message("done!")
  } #else if (!("Pathway_sclm_results" %in% names(Pagwas))) {
  #  stop("chrom_ld should input!")
 # }

  cat("Link_pathway_blocks_gwas: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
  cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)

  #############################
  ## 7.Pagwas_perform_regression
  #############################
  if (celltype) {
    message(paste(utils::timestamp(quiet = T), " ******* 7th: Celltype_heritability_contributions function start! ********", sep = ""))

    Pagwas<-celltype_f(Pagwas=Pagwas,iters=iters)

    message("done!")
    cat("Celltype_heritability_contributions: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
    cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)

    if (!singlecell) {
      return(Pagwas)
    }
  }
  #############################
  ## 8.scPagwas_perform_score
  #############################
  if (singlecell) {
    message(paste(utils::timestamp(quiet = T), " ******* 8th: scPagwas_perform_score function start! ********", sep = ""))
    tt <- Sys.time()

    Pagwas<- scPagwas_perform_regression(Pagwas,
      Pathway_ld_gwas_data = scPathway_ld_gwas_data,
      ncores = ncores)

    Pagwas <- scPagwas_perform_score(
      Pagwas = Pagwas,
      remove_outlier = TRUE
    )

    message("done!")
    cat("scPagwas_perform_score: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
    cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
    #############################
    ## 9.scGet_gene_heritability_correlation
    #############################
    message(paste(utils::timestamp(quiet = T), " ******* 9th: scGet_gene_heritability_correlation function start! ********", sep = ""))

    tt <- Sys.time()
    Pagwas <- scGet_gene_heritability_correlation(Pagwas = Pagwas)
    ############# output
    cat("scGet_gene_heritability_correlation: ", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
    cat(Sys.time() - tt, "\n", file = paste0("./", output.dirs, "/scPagwas.run.log"), append = T)
    message("done")
    #############################
    ## 10.Integrate the output write it to files
    #############################
    message(paste(utils::timestamp(quiet = T), " ******* 10th:Integrate the output write it to files! ********", sep = ""))

    Pagwas<-Pagwas_result_integtate(Pagwas=Pagwas,
                                    seruat_return=seruat_return,
                                    output.dirs=output.dirs,
                                    Single_data=Single_data,
                                    output.prefix=output.prefix,
                                    n_topgenes=n_topgenes,
                                    assay=assay)

    SOAR::Remove(SOAR::Objects())
    gc()
    return(Pagwas)
  }
}
