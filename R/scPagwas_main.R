#' scPagwas, A single-cell pathway-based principal component (PC)-scoring
#' algorithm named for integrating polygenic signals from GWAS with
#' single cell data to infer disease-relevant cell populations and score
#' the associations of individual cells with complex diseases.
#'
#'
#' @importFrom dplyr mutate filter inner_join %>%
#' @importFrom irlba irlba
#' @importFrom Seurat Assays FindVariableFeatures AverageExpression VariableFeatures GetAssayData RunPCA RunTSNE RunUMAP Embeddings CreateAssayObject DefaultAssay AddModuleScore VlnPlot
#' @importFrom SeuratObject Idents DefaultAssay GetAssayData SetAssayData
#' @importFrom Matrix Matrix colSums rowSums crossprod
#' @importFrom glmnet cv.glmnet
#' @importFrom GenomicRanges GRanges resize resize
#' @importFrom IRanges IRanges
#' @importFrom utils timestamp
#' @import ggplot2
#' @importFrom SOAR Store Objects Remove
#' @importFrom ggpubr ggscatter
#' @importFrom RMTstat qWishartMax
#' @importFrom gridExtra grid.arrange
#' @importFrom data.table setkey data.table as.data.table
#' @importFrom bigreadr fread2 rbind_df
#' @importFrom reshape2 dcast
#' @importFrom Rcpp sourceCpp
#' @importFrom bigmemory as.big.matrix is.big.matrix
#' @importFrom bigstatsr as_FBM big_apply big_univLinReg covar_from_df big_transpose big_cprodMat
#' @importFrom biganalytics apply
#' @importFrom graphics text
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom utils txtProgressBar setTxtProgressBar write.table write.csv timestamp globalVariables
#' @importFrom stats reorder cor qnorm density smooth.spline predict dnorm quantile offset
#' @importFrom methods as is
#' @title Main wrapper functions for scPagwas
#' @name scPagwas_main
#' @description Main Pagwas wrapper functions.
#' @details The entry point for Pagwas analysis.
#' Including the data input functions and the main progress functions;
#' It can also output the
#' running log and parameter log for scPagwas, and construct the folder
#' for output.
#'
#' @param Pagwas (list)default is "NULL" when you first run the funciton;
#' Pagwas
#' should be list class; Sometimes, It can inherit the result from
#' "scPagwas_main" function
#' last time, when you turn the "seurat_return" to FALSE; It is suitable
#' for circulation
#' running for the same single data.
#' @param gwas_data (data.frame)GWAS summary data; It must have some
#' colmuns such as:
#'  chrom|    pos    |   rsid    |   se  |  beta |  maf
#'     6 | 119968580 | rs1159767 | 0.032 | 0.019 |0.5275
#'    10 | 130566523 |  rs559109 | 0.033 | 0.045 |0.4047
#'     5 | 133328825 | rs6893145 | 0.048 | 0.144 |0.1222
#'     7 | 146652932 | rs13228798| 0.035 | 0.003 | 0.3211
#'
#' @param output.prefix (character)output.prefix,prefix used for output
#' tables; If the run_split if TRUE, the output.dirs should be identical for
#' different sub-data,and the output.prefix should be different.
#' @param marg (integr) the distance to TSS site,default is 10000, then
#' gene-TSS-window size is 20000.
#' @param block_annotation (data.frame) Start and end points for block
#'  traits, usually genes.
#' @param Single_data (character or seurat)Input the Single data in seurat
#' format, or the seurat data address for rds format.Idents should be the celltypes annotation.
#' @param assay (character)assay data of your single cell data to use,
#' default is "RNA"
#' @param Pathway_list (list,character) pathway gene sets list
#' @param chrom_ld (list,numeric)LD data for 22 chromosome.
#' @param run_split (logical) Whether the input single cell data is a split sub-data, if TRUE,
#' one result(gPas score) is return, if FALSE, the whole function is running. default is FALSE.
#' @param n.cores cores for regression
#' @param maf_filter (numeric)Filter the maf, default is 0.01
#' @param min_clustercells (integr)Only use is when FilterSingleCell is
#' TRUE.Threshold for total cells fo each cluster.default is 10
#' @param min.pathway.size (integr)Threshold for min pathway gene size.
#' default is 5
#' @param max.pathway.size (integr)Threshold for max pathway gene size.
#' default is 300
#' @param iters (integr)number of bootstrap iterations to perform
#' @param remove_outlier (logical)Whether to remove the outlier for
#' scPagwas score.
#' @param seurat_return (logical) Whether return the seurat format result,
#' if not,will return a list result;
#' @param singlecell (logical)Whether to produce the singlecell result;
#' @param celltype (logical)Whether to produce the celltypes result;
#' @param output.dirs (character)output directory for result; If the
#' directory is nonexistence,
#' the function will create it;
#' If the run_split if TRUE, the output.dirs should be identical for
#' different sub-data,and the output.prefix should be different.
#' @param n_topgenes (integr)Number of top associated gene selected to
#' calculate the scPagwas score;
#'
#' @return
#' Returns a seurat data with entries(seurat_return=T):
#' \describe{
#'   \item{assay:}{
#'   {scPagwasPaPca:}{An assay for S4 type of data; the svd result
#'   for pathways in each cells;}
#'   {scPagwasPaHeritability:}{An assay for S4 type of data; the
#'   gPas matrix for pathways in each cells;}}
#'   \item{meta.data}{
#'   {scPagwas.TRS.Score1:}{ the column for "meta.data";Enrichment socre
#'   for inheritance associated top genes.}
#'   {scPagwas.gPAS.score:}{ the column for "meta.data";Inheritance
#'   regression effects for each cells}
#'   {CellScalepValue:}{ the column for "meta.data";CellScalepValue for
#'   each cells}
#'   {CellScaleqValue:}{ the column for "meta.data";CellScaleqValue for
#'   each cells}}
#'   \item{misc: element in result,\code{Pagwas@misc }}{
#'   {Pathway_list:}{a list for pathway gene list intersecting with single
#'   cell data}
#'   {pca_cell_df:}{ a data frame for pathway pca result for each celltype.}
#'   {lm_results:}{ the regression result for each cell.}
#'   {gene_heritability_correlation:}{
#'   heritability correlation value for each gene;}
#'   {scPathways_scGene_rankP:}{pvaues for each pathway}
#'   {bootstrap_results:}{The bootstrap data frame results for celltypes
#'   including bootstrap pvalue and confidence interval.}
#' }
#'
#' }
#' Returns files:
#'
#' \describe{
#'   {scPagwas.run.log:}{ the running log file for scPagwas}
#'   {*_parameters.txt:}{parameters log file for scPagwas}
#'   {*_pca_scCell_mat.txt:}{ 1st svd result for pathways in each cell}
#'   {*_pca_celltypes_mat.txt:}{ 1st svd results for pathways in celltypes}
#'   {*_singlecell_scPagwas_score.Result.csv:}{ The final result for
#'   lm and top gene score}
#'   {*_celltypes_bootstrap_results.csv:}{The bootstrap data frame
#'   results for celltypes including bootstrap pvalue and confidence
#'   interval}
#'   {*_gene_heritability_correlation.csv:}{ heritability correlation
#'   value("cor" for pearson) for each gene;}

#' }
#'
#' Returns a list class with entries(seurat_return=F):
#' \describe{
#'   {scPagwasPaPca:}{Assays for S4 type of data; the svd result for
#'   pathways in each cells;}
#'   {scPagwas.topgenes.Score1:}{ the column for "meta.data";
#'   Enrichment socre for inheritance associated top genes.}
#'   {sclm_score:}{ the column for "meta.data";Inheritance regression
#'   effects for each cells}
#'   {Pathway_list:}{ The number of Lanczos iterations carried out}
#'   {pca_cell_df:}{ The total number of matrix vector products carried out}
#'   {sclm_results:}{ The total number of matrix vector products carried out}
#'   {allsnp_gene_heritability_correlation:}{ The total number of
#'   matrix vector products carried out}
#'   {Pathway_ctlm_results:}{ The total number of matrix vector products
#'   carried out}
#'   {lm_results:}{ The total number of matrix vector products carried out}
#'   {Pathway_ct_results:}{ The total number of matrix vector products
#'   carried out}
#' }
#'
#' @note
#' 1.When you run the package in linux server, you can run
#' \code{export OPENBLAS_NUM_THREADS=1}
#' before enter into R environment, It can help you to avoid core errors.
#' 2. When the code is kill for some wrong, you can run
#' \code{SOAR::Objects()} to check the objects, or
#' run \code{SOAR::Remove(SOAR::Objects())} to remove it.
#' 3. When your data is too big to run the function(bus error), you can
#' split your data into several sub-data and run the merge_pagwas
#' to merge them.
#' @references
#' Polygenic regression identifies trait-relevant cell subpopulations
#' through pathway activity transformed single-cell RNA sequencing data.(2022)
#' @export
#'
#' @examples
#' library(scPagwas)
#' Pagwas_data <- scPagwas_main(
#'   Pagwas = NULL,
#'   gwas_data = system.file("extdata", "GWAS_summ_example.txt",
#'     package = "scPagwas"
#'   ),
#'   Single_data = system.file("extdata", "scRNAexample.rds",
#'     package = "scPagwas"
#'   ),
#'   output.prefix = "test",
#'   output.dirs = "scPagwastest_output",
#'   block_annotation = block_annotation,
#'   assay = "RNA",
#'   Pathway_list = Genes_by_pathway_kegg,
#'   chrom_ld = chrom_ld,
#'   singlecell = TRUE,
#'   celltype = TRUE)
#' @author Chunyu Deng
#' @aliases scPagwas_main
#' @keywords scPagwas_main, wrapper of scPagwas functions.
scPagwas_main <- function(Pagwas = NULL,
                          gwas_data = NULL,
                          output.prefix = "Test",
                          output.dirs = "scPagwastest_output",
                          block_annotation = NULL,
                          Single_data = NULL,
                          assay = "RNA",
                          Pathway_list = NULL,
                          chrom_ld = NULL,
                          run_split=FALSE,
                          Correct_BG_p=FALSE,
                          n.cores=1,
                          marg = 10000,
                          maf_filter = 0.01,
                          min_clustercells = 10,
                          min.pathway.size = 5,
                          max.pathway.size = 300,
                          iters_celltype = 200,
                          iters_singlecell = 500,
                          n_topgenes = 1000,
                          singlecell = TRUE,
                          celltype = TRUE,
                          seurat_return = TRUE,
                          remove_outlier = TRUE) {

  # Pagwas =NULL
  # gwas_data = "/share/pub/dengcy/GWAS_Multiomics/bloodtraits/monocytecount_prune_gwas_data.txt"
  #Single_data ="/share/pub/dengcy/GWAS_Multiomics/modelgroudtruth/sim_data_8.16.rds"
  # output.prefix="Test"
  # output.dirs="Test"
  # #load("D:/tempdata/reduce_genes.by.regulatory.pathway.RData")
  # Pathway_list=reduce_genes.by.regulatory.pathway
  # assay="RNA"
  # singlecell=T
  # celltype=F
  # block_annotation = block_annotation
  # chrom_ld = chrom_ld
  #  assay = "RNA"
  # run_split=FALSE
  # n.cores=1
  # marg = 10000
  # pa_method = "SVD"
  # maf_filter = 0.01
  # min_clustercells = 10
  # min.pathway.size = 5
  # max.pathway.size = 300
  # iters_singlecell = 100
  # Correct_BG_p=TRUE
  # iters_celltype = 200
  # n_topgenes = 1000
  # seurat_return = TRUE
  # remove_outlier = TRUE

  #######
  ## initialize log-file
  ########

  if (!dir.exists(output.dirs)) {
    dir.create(output.dirs)
  }
  if (!dir.exists(output.dirs)) {
  dir.create(paste0("./", output.dirs, "/temp"))
  }
  ## miximal file path lenght;
  ## Windows OS support max. 259 characters

  Sys.setenv(R_LOCAL_CACHE = paste0("./", output.dirs, "/scPagwas_cache"))

  tt <- Sys.time()
  if (is.null(Pagwas)) {
    Pagwas <- list()
  } else if (class(Pagwas) == "Seurat" & is.null(Single_data)) {
    Single_data <- Pagwas
    Pagwas <- list()
    Pagwas <- Single_data@misc
    if ("scPagwasPaPca" %in% Seurat::Assays(Single_data)) {
      Pagwas$pca_scCell_mat <- GetAssayData(Single_data,
        assay = "scPagwasPaPca"
      )
    }
    if (assay %in% Seurat::Assays(Single_data)) {
      Pagwas$data_mat <- Seurat::GetAssayData(Single_data, assay = assay)
    } else {
      stop("Error:assay is not in Pagwas!")
    }

    if (is.null(Pathway_list)) {
      stop("Error:Pathway_list should be input!")
    }
    SOAR::Store(Single_data)
    Pagwas$rawPathway_list <- Pathway_list
  } else if (class(Pagwas) == "Seurat" & !is.null(Single_data)) {
    message("Warning:Single_data and Pagwas seurat class are redundant!
              we will keep the new Single_data and rerun the
            Single_data_input and Pathway_pcascore_run function")
    Pagwas <- list()
  } else if (class(Pagwas) == "list" & is.null(Single_data) & singlecell) {
    if("data_mat" %in% names(Pagwas) & run_split){
      message("The single cell data are in the preprocessed pagwas list!")
    }else{
      stop("Error:Single_data should be input!")
    }

  } else if (class(Pagwas) == "list" & !is.null(Single_data)) {
    Pagwas$rawPathway_list <- Pathway_list
  } else if (class(Pagwas) != "list") {
    stop("Error:The class for Pagwas is wrong! Should be NULL,
         list or Seurat class.")
  }

  #############################
  ## 1.Single_data_input
  #############################
  message(paste(utils::timestamp(quiet = T),
    " ******* 1st: Single_data_input function start! ********",
    sep = ""
  ))
  tt <- Sys.time()
  if (class(Single_data) == "character") {
    if (grepl(".rds", Single_data)) {
      message("** Start to read the single cell data!")
      Single_data <- readRDS(Single_data)
    } else {
      stop("Error:There is need a data in `.rds` format ")
    }
    if (!assay %in% Seurat::Assays(Single_data)) {
      stop("Error:There is no need assays in your Single_data")
    }
  } else if (class(Single_data) != "Seurat") {
    if(run_split){
      message("run_split is TRUE!")
    }else{
      stop("Error:When the run_split is FALSE! There is need a Seurat class for Single_data")
    }

  }


  if (is.null(Pathway_list)) {
    stop("Error:Pathway_list should be input!")
  }

  if(!("data_mat" %in% names(Pagwas))){
    message("** Start to filter single cell data!")
    Pagwas <- Single_data_input(
      Pagwas = Pagwas,
      assay = assay,
      # nfeatures =nfeatures,
      Single_data = Single_data,
      Pathway_list = Pathway_list,
      min_clustercells = min_clustercells
    )
  }

  # save the Single_data
  if(!run_split){
   Single_data <- Single_data[, colnames(Pagwas$data_mat)]
   SOAR::Store(Single_data)

  }else{
    if(!is.null(Single_data)){
      #Single_data <- Single_data[, colnames(Pagwas$data_mat)]
      rm(Single_data)
    }
    r_n <- colnames(Pagwas$data_mat)
  }
  message("done!")


  #############################
  ## 2.Pathway_pcascore_run
  #############################
  if (!("pca_cell_df" %in% names(Pagwas))) {
    message(paste(utils::timestamp(quiet = T),
      " ******* 2nd: Pathway_pcascore_run function start!! ********",
      sep = ""
    ))

    tt <- Sys.time()
      Pagwas <- Pathway_pcascore_run(
        Pagwas = Pagwas,
        Pathway_list = Pathway_list,
        min.pathway.size = min.pathway.size,
        max.pathway.size = max.pathway.size
      )

    message("done!")

  }

  #############################
  ## 3.GWAS_summary_input
  #############################
  message(paste(utils::timestamp(quiet = T),
    " ******* 3rd: GWAS_summary_input function start! ********",
    sep = ""
  ))

  if (!is.null(gwas_data)) {
    if (class(gwas_data) == "character") {
      message("** Start to read the gwas_data!")
      suppressMessages(gwas_data <- bigreadr::fread2(gwas_data))
    } else {
      stop("Error:There is need a filename and address for gwas_data")
    }

    if (maf_filter >= 1 & maf_filter < 0) {
      stop("Error:maf_filter should between 0 and 1")
    }

    tt <- Sys.time()
    Pagwas <- GWAS_summary_input(
      Pagwas = Pagwas,
      gwas_data = gwas_data,
      maf_filter = maf_filter
    )
    rm(gwas_data)
    message("done!")


    #############################
    ## 4.SnpToGene
    #############################
    message(paste(utils::timestamp(quiet = T),
      " ******* 4th: SnpToGene start!! ********",
      sep = ""
    ))

    tt <- Sys.time()
    if (!is.null(block_annotation)) {
      Pagwas$snp_gene_df <- SnpToGene(
        gwas_data = Pagwas$gwas_data,
        block_annotation = block_annotation,
        marg = marg
      )
    } else if (!("snp_gene_df" %in% names(Pagwas))) {
      stop("Error: block_annotation should input!")
    }

  } else if (!("gwas_data" %in% names(Pagwas)) | !("snp_gene_df" %in%
    names(Pagwas))) {
    stop("Error: gwas_data should be input!")
  }
  #############################
  ## 5.Pathway_annotation_input
  #############################
  message(paste(utils::timestamp(quiet = T),
    " ******* 5th: Pathway_annotation_input function
                start! ********",
    sep = ""
  ))

  tt <- Sys.time()
  if (!is.null(block_annotation)) {
    Pagwas <- Pathway_annotation_input(
      Pagwas = Pagwas,
      block_annotation = block_annotation
    )
  } else if (!("snp_gene_df" %in% names(Pagwas))) {
    stop("Error: block_annotation should input!")
  }

  message("done!")


  #############################
  ## 6.Link_pathway_blocks_gwas
  #############################
  message(paste(utils::timestamp(quiet = T),
    " ******* 6th: Link_pathway_blocks_gwas
                function start! ********",
    sep = ""
  ))

  tt <- Sys.time()
  if (!is.null(chrom_ld)) {

    Pagwas <- Link_pathway_blocks_gwas(
      Pagwas = Pagwas,
      chrom_ld = chrom_ld,
      singlecell = singlecell,
      celltype = celltype,
      backingpath=paste0("./", output.dirs, "/temp"),
      n.cores=n.cores
      )

    message("done!")

    if (file.exists(paste0("./", output.dirs, "/temp"))) {
      unlink(paste0("./", output.dirs, "/temp"),recursive = TRUE)
      cat("bk files has been deleted")
    }


  } else if (!("Pathway_sclm_results" %in% names(Pagwas))) {
    stop("Error: chrom_ld should input!")
  }


  #############################
  ## 7.Pagwas_perform_regression
  #############################
  if (celltype) {
    message(paste(utils::timestamp(quiet = T),
      " ******* 7th: Celltype_heritability_contributions
                  function start! ********",
      sep = ""
    ))

    Pagwas$lm_results <- Pagwas_perform_regression(
      Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data
    )
    Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = iters_celltype, part = 0.5)

    Pagwas$Pathway_ld_gwas_data <- NULL

    message("done!")

    if (!singlecell) {
      SOAR::Remove(SOAR::Objects())
      return(Pagwas)
    }
    Pagwas$merge_scexpr <- NULL
  }
  #############################
  ## 8.scPagwas_perform_score
  #############################
  if (singlecell) {
    message(paste(utils::timestamp(quiet = T),
      " ******* 8th: scPagwas_perform_score function start! ********",
      sep = ""
    ))
    Pagwas$Pathway_ld_gwas_data <- NULL
    Pagwas <- scPagwas_perform_score(
      Pagwas = Pagwas,
      remove_outlier = TRUE
    )
    message("done!")

  #############################
  ## 9.scGet_gene_heritability_correlation
  #############################
  if(!run_split){
  message(paste(utils::timestamp(quiet = T),
    " ******* 9th: scGet_gene_heritability_correlation function start! ********",
    sep = ""
  ))

  tt <- Sys.time()

  Pagwas$gene_heritability_correlation <- scGet_gene_heritability_correlation(scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score,
                                                                              data_mat=Pagwas$data_mat)

  scPagwas_topgenes <- names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation, decreasing = T), ])[1:n_topgenes]

  message("done")


  #############################
  ## 10.output
  #############################

  utils::write.table(Pagwas$Pathway_sclm_results,
    file = paste0(
      "./", output.dirs, "/",
      output.prefix, "_Pathway_singlecell_lm_results.txt"
    ),
    quote = F, sep = "\t"
  )
  utils::write.csv(Pagwas$bootstrap_results,
    file = paste0(
      "./", output.dirs, "/",
      output.prefix,
      "_celltypes_bootstrap_results.csv"
    ),
    quote = F
  )
  utils::write.csv(Pagwas$scPathways_scGene_rankP,
    file = paste0(
      "./", output.dirs, "/",
      output.prefix,
      "_singlecell_Pathways_scGene_rankP.csv"
    ),
    quote = F
  )
  utils::write.csv(Pagwas$gene_heritability_correlation,
    file = paste0(
      "./", output.dirs, "/",
      output.prefix,
      "_gene_heritability_correlation.csv"
    ),
    quote = F
  )

  message("* Get scaled P for each single cell")
  CellScalepValue <- scGene_scaleP(
    Single_mat = Pagwas$data_mat[scPagwas_topgenes, ]
    )
  Pagwas[c(
    "VariableFeatures", "merge_scexpr", "snp_gene_df",
    "rawPathway_list", "data_mat"
  )] <- NULL


  Single_data <- Seurat::AddModuleScore(Single_data,
    assay = assay,
    list(scPagwas_topgenes),
    name = c("scPagwas.TRS.Score")
  )

  if(Correct_BG_p){
    message("* Get Random Correct background pvalue for each single cell!")
    correct_pdf<-Get_CorrectBg_p(Single_data=Single_data,
                                 scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                                 iters_singlecell=iters_singlecell,
                                 n_topgenes=n_topgenes,
                                 scPagwas_topgenes=scPagwas_topgenes
    )
    Pagwas$Random_Correct_BG_pdf <- correct_pdf
    message("* Get Merged pvalue for each celltype!")
    Pagwas$Merged_celltype_pvalue<-Merge_celltype_p(single_p=correct_pdf$pooled_p,celltype=Pagwas$Celltype_anno$annotation)
  }
  Pagwas$scPagwas.TRS.Score <- Single_data$scPagwas.TRS.Score1
  Pagwas$CellScalepValue <- CellScalepValue
    # CellScalepValue
  if(!Correct_BG_p){
    a <- data.frame(
      scPagwas.TRS.Score = Pagwas$scPagwas.TRS.Score,
      scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
      pValueHighScale = Pagwas$CellScalepValue$pValueHigh,
      qValueHighScale = Pagwas$CellScalepValue$qValueHigh
    )
  }else if(Correct_BG_p){
    a <- data.frame(
      scPagwas.TRS.Score = Pagwas$scPagwas.TRS.Score,
      scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
      pValueHighScale = Pagwas$CellScalepValue$pValueHigh,
      qValueHighScale = Pagwas$CellScalepValue$qValueHigh,
    Random_Correct_BG_p = correct_pdf$pooled_p,
    Random_Correct_BG_adjp = correct_pdf$adj_p,
    Random_Correct_BG_z = correct_pdf$pooled_z)
    utils::write.csv(Pagwas$Merged_celltype_pvalue,
                     file = paste0(
                       "./", output.dirs, "/", output.prefix,
                       "_Merged_celltype_pvalue.csv"
                     ),
                     quote = F
    )
  }


    utils::write.csv(a,
      file = paste0(
        "./", output.dirs, "/",
        output.prefix,
        "_singlecell_scPagwas_score_pvalue.Result.csv"
      ),
      quote = F
    )
    if (!seurat_return) {
    SOAR::Remove(SOAR::Objects())
    return(Pagwas)
  } else {
    Pagwas[c(
      "snp_gene_df",
      "Pathway_sclm_results", "CellScalepValue",
      "scPagwas.TRS.Score"
    )] <- NULL
    # Pagwas[c()] <- NULL
    scPagwas_pathway <- SeuratObject::CreateAssayObject(data = Pagwas$Pathway_single_results)
    scPagwas_pca <- SeuratObject::CreateAssayObject(data = Pagwas$pca_scCell_mat)

    Single_data[["scPagwasPaHeritability"]] <- scPagwas_pathway
    Single_data[["scPagwasPaPca"]] <- scPagwas_pca

    rm(scPagwas_pathway)
    rm(scPagwas_pca)

    Single_data$scPagwas.gPAS.score <- Pagwas$scPagwas.gPAS.score[rownames(Pagwas$Celltype_anno)]
    Single_data$ScalepValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "pValueHigh"]
    Single_data$ScaleqValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "qValueHigh"]
      if(Correct_BG_p){
    Single_data$Random_Correct_BG_p <- correct_pdf$pooled_p
    Single_data$Random_Correct_BG_adjp <- correct_pdf$adj_p
    Single_data$Random_Correct_BG_z <- correct_pdf$pooled_z
      }
    Pagwas[c(
      "scPagwas.gPAS.score", "Celltype_anno",
      "Pathway_single_results", "pca_scCell_mat"
    )] <- NULL
    Single_data@misc <- Pagwas
    rm(Pagwas)
    SOAR::Remove(SOAR::Objects())
    return(Single_data)

  }
  }else{
    if(all(r_n==names(Pagwas$scPagwas.gPAS.score))){

    df<-data.frame(cellnames=r_n,
                scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score)
    utils::write.csv(df,
    file = paste0(
      "./", output.dirs, "/", output.prefix,
      "_singlecell_scPagwas.gPAS.score.Result.csv"
    ),
    quote = F
    )
    }else{
      warning("The single cell' name may be changed for scPagwas.gPAS.score! It will influence
              the result of merge function")
      df<-data.frame(cellnames=names(Pagwas$scPagwas.gPAS.score),
                      scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score)
      utils::write.csv(df,
                        file = paste0(
                          "./", output.dirs, "/", output.prefix,
                          "_singlecell_scPagwas.gPAS.score.Result.csv"
                        ),
                        quote = F
      )

      }
    }
  }
  if (file.exists(paste0("./", output.dirs, "/scPagwas_cache"))) {
    unlink(paste0("./", output.dirs, "/scPagwas_cache"),recursive = TRUE)
  }
}


#' merge_pagwas
#' @description Merging several pagwas result to one result.When your single
#' cell data is too large to run, you can split your data into several
#' split, then use this function to merge them. Split the data randomly or
#' in order of cellytpes.
#' @note 1.The rownames genes and gwas data must uniform.
#' 2.The merge_pagwas function only can use to merge the seurat format
#' result.
#' That is "singlecell = T;celltype =F" can use to all the split data,
#' "singlecell = T;celltype = T" can be used only when the data is splited
#' by celltypes.
#' but "singlecell = F;celltype =T" can not be used.
#'
#' @param Single_data The Whole single cell data set. rds files in seurat format.
#' The filesnames or seurat format data.
#' @param output.dirs The files names for the result of scPagwas_main function.
#' @param n_topgenes (integr)Number of top associated gene selected to
#' calculate the scPagwas score;
#' @param assay The assay used for single cell data.
#' @param random Whether split the Single_data in random.
#' Using id when the Single_data is too big to run.
#' @param seed The seed for random.
#' @param random_times The times for random, default is 10.
#' @param proportion The proportion for spliting the Single_data; the bigger for
#' Single_data the less proportion for spliting, and the more times for random_times.
#'
#' @return seurat format result for scPagwas.
#' @export
#'
#' @examples
#' library(Seurat)
#' scRNAexample <-readRDS(system.file("extdata", "scRNAexample.rds", package = "scPagwas"))
#' #set the output.dirs for any directory.
#' set.dirs="D:/OneDrive/GWAS_Multiomics/Pagwas/test/"
#' #set the number of split time, it depend on the size of your single cell data.
#' #There have two form of spliting the data:
#'
#' #########
#' ##1.By Random
#' #Create the random index number.
#' n_split=2
#' Split_index <- rep(1:n_split, time = ceiling(ncol(scRNAexample)/n_split), length = ncol(scRNAexample))
#'
#' for (i in 1:n_split) {
#'   Example_splice <- scRNAexample[,Split_index==i]
#'   saveRDS(Example_splice,file = paste0(set.dirs,"Example_splice",i,".rds"))
#' }
#'
#' #########
#' ###1.By celltypes
#' #Create the random index number.
#' for (i in levels(Idents(scRNAexample))) {
#'   Example_splice <- subset(scRNAexample,idents = i)
#'   saveRDS(Example_splice,file = paste0(set.dirs,"Example_splice",i,".rds"))
#' }
#' Pagwas<-list()
#' gwas_data <- bigreadr::fread2(system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"))
#' Pagwas <- GWAS_summary_input(
#'   Pagwas = Pagwas,
#'   gwas_data = gwas_data,
#'   maf_filter = 0.1
#' )
#' Pagwas$snp_gene_df <- Snp2Gene(snp = Pagwas$gwas_data, refGene = block_annotation, marg = 10000)
#' names(Pagwas)
#' setwd(set.dirs)
#'
#' for (i in 1:n_split) {
#'   scPagwas_main(Pagwas =Pagwas,
#'                 gwas_data =NULL,
#'                 Single_data = paste0("Example_splice",i,".rds"), #the
#'                 output.prefix=i, #the prefix can be the circulation coefficient i.
#'                 output.dirs="splice_scPagwas", #the output.dirs shoukd be the same for each circulation
#'                 Pathway_list=Genes_by_pathway_kegg,
#'                 run_split=TRUE, #This parameter is set to run single result for each split result.
#'                 min_clustercells=1, #the minimum cluster cell number should be 1.
#'                 assay="RNA",
#'                 block_annotation = block_annotation,
#'                 chrom_ld = chrom_ld)
#' }
#' if(length(SOAR::Objects())>0){
#'   SOAR::Remove(SOAR::Objects())
#' }
#'
#' setwd(set.dirs) #you must set the workfiles the same as before
#' Pagwas_integrate <- merge_pagwas(Single_data = system.file("extdata", "scRNAexample.rds", package = "scPagwas") ,# read the whole single cell data.
#'                                 output.dirs="splice_scPagwas", #The same with scPagwas_main functionn_topgenes = 1000,
#'                                 assay='RNA',
#'                                 random=T,
#'                                 seed=1234,
#'                                 random_times=10,
#'                                 proportion=0.3)
#'
#'
#' @author Chunyu Deng
#' @aliases merge_pagwas
#' @keywords merge_pagwas, integrate the spliced scPagwas result.

merge_pagwas <- function(Single_data = NULL,
                         n_topgenes = 1000,
                         output.dirs="splice_scPagwas",
                         assay='RNA',
                         random=FALSE,
                         seed=1234,
                         Correct_BG_p=FALSE,
                         iters_singlecell=1000,
                         random_times=10,
                         proportion=0.3) {

  if (class(Single_data) == "character") {
    if (grepl(".rds", Single_data)) {
      message("** Start to read the single cell data!")
      Single_data <- readRDS(Single_data)
    } else {
      stop("Error:There is need a data in `.rds` format ")
    }
  } else if (class(Single_data) != "Seurat") {
    stop("Error:There is need a Seurat class for Single_data")
  }

  if (!assay %in% Seurat::Assays(Single_data)) {
    stop("Error:There is no need assays in your Single_data")
  }

  oriDir <- paste0(getwd(),"./",output.dirs)
  files <- list.files(oriDir, pattern="*_singlecell_scPagwas.gPAS.score.Result.csv")

  scPagwas.gPAS.score<-unlist(lapply( files,function(file){
    gs<-read.csv(file=paste0(oriDir,"/",file))
    ga <- gs$scPagwas.gPAS.score
    names(ga) <- gs$cellnames
    return(ga)
  }))
  #colnames(Single_data)
  ins<-intersect(names(scPagwas.gPAS.score),colnames(Single_data))

  if(is.null(ins)){
    stop("The names of scPagwas.gPAS.score are not the same with colnames of Single_data,
         please check the names")
  }
  if(length(ins)<ncol(Single_data)){
    warning("The total length of the result of scPagwas.gPAS.score for your splited single cell data
            is not the same with your Single_data, Please check the code before or the min_clustercells
            parameter(set to 1)")
  }
  Single_data <- Single_data[,names(scPagwas.gPAS.score)]
  Single_data$scPagwas.gPAS.score<-scPagwas.gPAS.score

  if(!random){
  scPagwas.gPAS.score <- data.matrix(Single_data$scPagwas.gPAS.score)
  sparse_cor <- corSparse(
    X = t(as_matrix(GetAssayData(Single_data, slot = "data", assay = "RNA"))),
    Y = scPagwas.gPAS.score
  )
  #rownames(sparse_cor) <- rownames(Single_data)
  colnames(sparse_cor) <- "gene_heritability_correlation"
  sparse_cor[is.nan(sparse_cor)] <- 0
  Single_data@misc$gene_heritability_correlation <- sparse_cor
  scPagwas_topgenes <- names(sparse_cor[order(sparse_cor, decreasing = T), ])[1:n_topgenes]
  Single_data <- Seurat::AddModuleScore(Single_data,
                                         assay = "RNA",
                                         list(scPagwas_topgenes),
                                         name = c("scPagwas.TRS.Score")
  )

  ########### get the CellScalepValue

  CellScalepValue <- scGene_scaleP(
    Single_mat = t(data.matrix(
      GetAssayData(Single_data, assay = "RNA")[scPagwas_topgenes, ]
    ))
  )
  Single_data$CellScalepValue <- CellScalepValue[, "pValueHigh"]
  Single_data$CellScaleqValue <- CellScalepValue[, "qValueHigh"]

  }else{
    message("Start to run the corSparse function in random!")
    set.seed(seed)
    sparse_cor_list<-list()
    for (j in 1:random_times) {
      print(paste0("Randome Times: ",j))
    index<-sample(1:ncol(Single_data),ceiling(ncol(Single_data)*proportion))
    #Single_data[,index]

    scPagwas.gPAS.score <- data.matrix(Single_data$scPagwas.gPAS.score[index])
    sparse_cor <- corSparse(
      X = t(as_matrix(GetAssayData(Single_data[,index], slot = "data", assay = "RNA"))),
      Y = scPagwas.gPAS.score
    )
    #rownames(sparse_cor) <- colnames(Single_data)[index]
    colnames(sparse_cor) <- "gene_heritability_correlation"
    sparse_cor[is.nan(sparse_cor)] <- 0
    sparse_cor_list[[j]]<-unlist(sparse_cor[,1])
    }
    sparse_cor_list<-as.data.frame(sparse_cor_list)
    colnames(sparse_cor_list)<-paste0("Cor_coefficient.ramdom",1:random_times)
    sparse_cor<- apply(sparse_cor_list,1,mean)
    Single_data@misc$sparse_cor_df <- sparse_cor_list
    Single_data@misc$gene_heritability_correlation <- data.matrix(sparse_cor)
    colnames(Single_data@misc$gene_heritability_correlation)<-"gene_heritability_correlation"
    scPagwas_topgenes <- names(sparse_cor[order(sparse_cor, decreasing = T)])[1:n_topgenes]
    Single_data <- Seurat::AddModuleScore(Single_data,
                                          assay = "RNA",
                                          list(scPagwas_topgenes),
                                          name = c("scPagwas.TRS.Score")
    )

    message("Start to run the scGene_rankP function in random!")
    pv<-list()
    qv<-list()
    for (j in 1:random_times) {
      print(paste0("Randome Times: ",j))
      index<-sample(1:ncol(Single_data),ceiling(ncol(Single_data)*proportion))
      CellScalepValue <- scGene_scaleP(
        Single_mat = t(data.matrix(
        GetAssayData(Single_data[,index], assay = "RNA")[scPagwas_topgenes, ]
         ))
       )
      a<-CellScalepValue[, "pValueHighScale"]
      names(a)<-rownames(CellScalepValue)
      pv[[j]]<-a
      a<-CellScalepValue[, "qValueHighScale"]
      names(a)<-rownames(CellScalepValue)
      qv[[j]]<-a
    }

    Single_data$CellScalepValue <- unlist(pv)[colnames(Single_data)]
    Single_data$CellScaleqValue <- unlist(qv)[colnames(Single_data)]

  }

  if(Correct_BG_p){
    message("* Get Random Correct background pvalue for each single cell!")
    correct_pdf<-Get_CorrectBg_p(Single_data=Single_data,
                                 scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                                 iters_singlecell=iters_singlecell,
                                 n_topgenes=n_topgenes,
                                 scPagwas_topgenes=scPagwas_topgenes
    )
    #Pagwas$Random_Correct_BG_pdf <- correct_pdf
    Single_data$Random_Correct_BG_p <- correct_pdf$pooled_p
    Single_data$Random_Correct_BG_adjp <- correct_pdf$adj_p
    Single_data$Random_Correct_BG_z <- correct_pdf$pooled_z
    Single_data@misc$Merged_celltype_pvalue<-Merge_celltype_p(single_p=correct_pdf$adj_p,celltype=Pagwas$Celltype_anno$annotation)

  }
  ########### get the scPagwas.TRS.Score
  return(Single_data)
}

