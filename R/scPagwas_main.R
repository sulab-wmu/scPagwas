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
#' @importFrom bigstatsr as_FBM big_apply big_univLinReg covar_from_df big_transpose big_cprodMat
#' @importFrom RMTstat qWishartMax
#' @importFrom gridExtra grid.arrange
#' @importFrom data.table setkey data.table as.data.table
#' @importFrom bigreadr fread2 rbind_df
#' @importFrom reshape2 dcast
#' @importFrom bigmemory as.big.matrix is.big.matrix
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
#' last time, when you turn the "seruat_return" to FALSE; It is suitable
#' for circulation
#' running for the same single data.
#' @param gwas_data (data.frame)GWAS summary data; It must have some
#' colmuns such as:
#'  chrom|    pos    |   rsid    |   se  |  beta |  maf
#'     6 | 119968580 | rs1159767 | 0.032 | 0.019 |0.5275
#'    10 | 130566523 |  rs559109 | 0.033 | 0.045 |0.4047
#'     5 | 133328825 | rs6893145 | 0.048 | 0.144 |0.1222
#'     7 | 146652932 | rs13228798| 0.035 | 0.003 | 0.3211
#' @param output.prefix (character)output.prefix,prefix used for output
#' tables;
#' @param marg (integr) the distance to TSS site,default is 10000, then
#' gene-TSS-window size is 20000.
#' @param block_annotation (data.frame) Start and end points for block
#'  traits, usually genes.
#' @param Single_data (character or Seruat)Input the Single data in seruat
#' format, or the Seruat data address for rds format.Idents should be the celltypes annotation.
#' @param assay (character)assay data of your single cell data to use,
#' default is "RNA"
#' @param Pathway_list (list,character) pathway gene sets list
#' @param chrom_ld (list,numeric)LD data for 22 chromosome.
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
#' @param ncores (integr)Parallel cores,default is 1. use detectCores()
#' to check the cores in computer.
#' @param seruat_return (logical) Whether return the Seruat format result,
#' if not,will return a list result;
#' @param singlecell (logical)Whether to produce the singlecell result;
#' @param celltype (logical)Whether to produce the celltypes result;
#' @param output.dirs (character)output directory for result; If the
#' directory is nonexistence,
#' the function will create it;
#' @param n_topgenes (integr)Number of top associated gene selected to
#' calculate the scPagwas score;
#'
#' @return
#' Returns a Seruat data with entries(seruat_return=T):
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
#'   {scPathways_rankPvalue:}{pvaues for each pathway}
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
#' Returns a list class with entries(seruat_return=F):
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
#'   celltype = TRUE,
#'   ncores = 1
#' )
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
                          ncores = 1) {

  #######
  ## initialize log-file
  ########
  if (!dir.exists(output.dirs)) {
    dir.create(output.dirs)
  }
  cat("## start at:", format(Sys.time()), "\n",
    file = paste0("./", output.dirs, "/scPagwas.run.log")
  )
  ## miximal file path lenght;
  ## Windows OS support max. 259 characters

  param.str <- c(
    paste("##", Sys.time()),
    paste("input gwas data files: ", gwas_data, sep = "\t"),
    paste0(
      "Single_data files: ",
      if (class(Single_data) == "Seurat") dim(Single_data) else Single_data
    ),
    paste("assay in single data: ", assay, sep = "\t"),
    paste("Pathway length: ", length(Pathway_list),
      collapse = " ",
      sep = "\t"
    ),
    paste("marg distance from TSS: ", marg, sep = "\t"),
    paste("output directory is: ", output.dirs, sep = "\t"),
    paste("Whether to run singlecell step: ", singlecell, sep = "\t"),
    paste("Whether to run celltype step: ", celltype, sep = "\t"),
    paste("The number of top genes for TRS score: ", n_topgenes, sep = "\t"),
    paste("maf_filter: ", maf_filter, sep = "\t"),
    paste("min_clustercells: ", min_clustercells, sep = "\t"),
    paste("min.pathway.size: ", min.pathway.size, sep = "\t"),
    paste("max.pathway.size: ", max.pathway.size, sep = "\t"),
    paste("remove_outlier: ", remove_outlier, sep = "\t"),
    paste("iters: ", iters, sep = "\t"),
    paste("ncores: ", ncores, sep = "\t")
  )
  writeLines(param.str, con = paste0(
    "./", output.dirs, "/", output.prefix,
    "_parameters.txt"
  ))
  # }

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
      Pagwas$data_mat <- GetAssayData(Single_data, assay = assay)
    } else {
      stop("Error:assay is not in Pagwas!")
    }

    if (is.null(Pathway_list)) {
      stop("Error:Pathway_list should be input!")
    }
    SOAR::Store(Single_data)
    Pagwas$rawPathway_list <- Pathway_list
  } else if (class(Pagwas) == "Seurat" & !is.null(Single_data)) {
    message("Warning:Single_data and Pagwas seruat class are redundant!
              we will keep the new Single_data and rerun the
            Single_data_input and Pathway_pcascore_run function")
    Pagwas <- list()
  } else if (class(Pagwas) == "list" & is.null(Single_data) & singlecell) {
    stop("Error:Single_data should be input!, the same as Pagwas")
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
  } else if (class(Single_data) != "Seurat") {
    stop("Error:There is need a Seurat class for Single_data")
  }

  if (!assay %in% Seurat::Assays(Single_data)) {
    stop("Error:There is no need assays in your Single_data")
  }
  if (is.null(Pathway_list)) {
    stop("Error:Pathway_list should be input!")
  }

  Pagwas <- Single_data_input(
    Pagwas = Pagwas,
    assay = assay,
    # nfeatures =nfeatures,
    Single_data = Single_data,
    Pathway_list = Pathway_list,
    min_clustercells = min_clustercells
  )

  Single_data <- Single_data[, colnames(Pagwas$data_mat)]
  # save the Single_data
  SOAR::Store(Single_data)
  message("done!")


  cat("Single_data import: ",
    file = paste0(
      "./", output.dirs,
      "/scPagwas.run.log"
    ),
    append = T
  )
  cat(Sys.time() - tt, "\n", file = paste0(
    "./", output.dirs,
    "/scPagwas.run.log"
  ), append = T)

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
    cat("Pathway_pcascore_run: ",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
    cat(Sys.time() - tt, "\n",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
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
    cat("GWAS_summary_input: ",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
    cat(Sys.time() - tt, "\n",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )

    #############################
    ## 4.Snp2Gene
    #############################
    message(paste(utils::timestamp(quiet = T),
      " ******* 4th: Snp2Gene start!! ********",
      sep = ""
    ))

    tt <- Sys.time()
    if (!is.null(block_annotation)) {
      Pagwas$snp_gene_df <- Snp2Gene(
        snp = Pagwas$gwas_data,
        refGene = block_annotation,
        marg = marg
      )
    } else if (!("snp_gene_df" %in% names(Pagwas))) {
      stop("Error: block_annotation should input!")
    }
    cat("Snp2Gene: ",
      file = paste0("./", output.dirs, "/scPagwas.run.log"),
      append = T
    )
    cat(Sys.time() - tt, "\n",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
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
  cat("Pathway_annotation_input: ",
    file = paste0(
      "./", output.dirs,
      "/scPagwas.run.log"
    ),
    append = T
  )
  cat(Sys.time() - tt, "\n",
    file = paste0(
      "./", output.dirs,
      "/scPagwas.run.log"
    ),
    append = T
  )


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
      ncores = ncores
    )

    message("done!")
  } else if (!("Pathway_sclm_results" %in% names(Pagwas))) {
    stop("Error: chrom_ld should input!")
  }

  cat("Link_pathway_blocks_gwas: ",
    file = paste0(
      "./", output.dirs,
      "/scPagwas.run.log"
    ),
    append = T
  )
  cat(Sys.time() - tt, "\n",
    file = paste0(
      "./", output.dirs,
      "/scPagwas.run.log"
    ),
    append = T
  )

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
    Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = iters, part = 0.5)

    Pagwas$Pathway_ld_gwas_data <- NULL

    message("done!")
    cat("Celltype_heritability_contributions: ",
      file = paste0(
        "./",
        output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
    cat(Sys.time() - tt, "\n",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )

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
    tt <- Sys.time()
    Pagwas$Pathway_ld_gwas_data <- NULL
    Pagwas <- scPagwas_perform_score(
      Pagwas = Pagwas,
      remove_outlier = TRUE
    )
    message("done!")
    cat("scPagwas_perform_score: ",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
    cat(Sys.time() - tt, "\n",
      file = paste0(
        "./", output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
    #############################
    ## 9.scGet_gene_heritability_correlation
    #############################
    message(paste(utils::timestamp(quiet = T),
      " ******* 9th: scGet_gene_heritability_correlation function start! ********",
      sep = ""
    ))

    tt <- Sys.time()
    Pagwas <- scGet_gene_heritability_correlation(Pagwas = Pagwas)
    scPagwas_topgenes <- names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation, decreasing = T), ])[1:n_topgenes]

    message("done")
    cat("9.scGet_gene_heritability_correlation: ",
      file = paste0(
        "./",
        output.dirs,
        "/scPagwas.run.log"
      ),
      append = T
    )
    cat(Sys.time() - tt, "\n",
      file = paste0(
        "./",
        output.dirs, "/scPagwas.run.log"
      ),
      append = T
    )

    #############################
    ## 9.output
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
    utils::write.csv(Pagwas$scPathways_rankPvalue,
      file = paste0(
        "./", output.dirs, "/",
        output.prefix,
        "_singlecell_Pathways_rankPvalue.csv"
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

    Pagwas[c(
      "VariableFeatures", "merge_scexpr", "snp_gene_df",
      "rawPathway_list", "data_mat"
    )] <- NULL


    Single_data <- Seurat::AddModuleScore(Single_data,
      assay = assay,
      list(scPagwas_topgenes),
      name = c("scPagwas.TRS.Score")
    )

    message("* Get rankPvalue for each single cell")
    CellScalepValue <- rankPvalue(
      datS = t(data.matrix(
        GetAssayData(Single_data, assay = assay)[scPagwas_topgenes, ]
      )),
      pValueMethod = "scale"
    )

    if (!seruat_return) {
      Pagwas$scPagwas.TRS.Score <- Single_data$scPagwas.TRS.Score1
      Pagwas$CellScalepValue <- CellScalepValue

      # CellScalepValue
      a <- data.frame(
        scPagwas.TRS.Score = Pagwas$scPagwas.TRS.Score,
        scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
        pValueHighScale = Pagwas$CellScalepValue$pValueHighScale,
        qValueHighScale = Pagwas$CellScalepValue$qValueHighScale
      )
      utils::write.csv(a,
        file = paste0(
          "./", output.dirs, "/",
          output.prefix,
          "_singlecell_scPagwas_score_pvalue.Result.csv"
        ),
        quote = F
      )

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
      Single_data$CellScalepValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "pValueHighScale"]
      Single_data$CellScaleqValue <- CellScalepValue[rownames(Pagwas$Celltype_anno), "qValueHighScale"]
      Pagwas[c(
        "scPagwas.gPAS.score", "Celltype_anno",
        "Pathway_single_results", "pca_scCell_mat"
      )] <- NULL
      Single_data@misc <- Pagwas
      rm(Pagwas)
      utils::write.csv(Single_data@meta.data[, c(
        "scPagwas.gPAS.score",
        "CellScalepValue",
        "scPagwas.TRS.Score1"
      )],
      file = paste0(
        "./", output.dirs, "/", output.prefix,
        "_singlecell_scPagwas_score_pvalue.Result.csv"
      ),
      quote = F
      )
      SOAR::Remove(SOAR::Objects())
      return(Single_data)
    }
  }
}


#' merge_pagwas
#' @description Merging several pagwas result to one result.When your single
#' cell data is too large to run, you can split your data into several
#' split, then use this function to merge them. Split the data randomly or
#' in order of cellytpes.
#' @note 1.The rownames genes and gwas data must uniform.
#' 2.The merge_pagwas function only can use to merge the seruat format
#' result.
#' That is "singlecell = T;celltype =F" can use to all the split data,
#' "singlecell = T;celltype = T" can be used only when the data is splited
#' by celltypes.
#' but "singlecell = F;celltype =T" can not be used.
#' @param Pagwas_list list format, element including seruat format result.
#' @param n_topgenes (integr)Number of top associated gene selected to
#' calculate the scPagwas score;
#'
#' @return seruat format result for scPagwas.
#' @export
#'
#' @examples
#' library(scPagwas)
#' Pagwas <- list()
#' library(Seurat)
#' scRNAexample <- readRDS(system.file("extdata", "scRNAexample.rds",
#'   package = "scPagwas"
#' ))
#' Example_splice1 <- subset(scRNAexample, idents = c(
#'   "Celltype1", "Celltype2", "Celltype3",
#'   "Celltype4", "Celltype5"
#' ))
#' Example_splice2 <- subset(scRNAexample, idents = c(
#'   "Celltype6", "Celltype7",
#'   "Celltype8", "Celltype9", "Celltype10"
#' ))
#'
#' gwas_data <- bigreadr::fread2(system.file("extdata",
#'   "GWAS_summ_example.txt",
#'   package = "scPagwas"
#' ))
#' Pagwas <- GWAS_summary_input(
#'   Pagwas = Pagwas,
#'   gwas_data = gwas_data,
#'   maf_filter = 0.1
#' )
#' Pagwas$snp_gene_df <- Snp2Gene(
#'   snp = Pagwas$gwas_data,
#'   refGene = block_annotation, marg = 10000
#' )
#' # Run the first split data
#' Pagwas1 <- scPagwas_main(
#'   Pagwas = Pagwas,
#'   gwas_data = NULL,
#'   Single_data = Example_splice1,
#'   output.prefix = "splice1",
#'   output.dirs = "splice_scPagwas",
#'   Pathway_list = Genes_by_pathway_kegg,
#'   ncores = 5,
#'   assay = "RNA",
#'   block_annotation = block_annotation,
#'   chrom_ld = chrom_ld
#' )
#' Pagwas2 <- scPagwas_main(
#'   Pagwas = Pagwas,
#'   gwas_data = NULL,
#'   Single_data = Example_splice2,
#'   output.prefix = "splice2 ",
#'   output.dirs = "splice_scPagwas",
#'   Pathway_list = Genes_by_pathway_kegg,
#'   ncores = 5,
#'   assay = "RNA",
#'   block_annotation = block_annotation,
#'   chrom_ld = chrom_ld
#' )
#' Pagwas_integrate <- merge_pagwas(
#'   Pagwas_list = list(Pagwas1, Pagwas2),
#'   n_topgenes = 1000
#' )
#' @author Chunyu Deng
#' @aliases merge_pagwas
#' @keywords merge_pagwas, integrate the spliced scPagwas result.

merge_pagwas <- function(Pagwas_list = NULL,
                         n_topgenes = 1000) {
  if (length(Pagwas_list) < 2) {
    stop("Pagwas_list should have more than two elements!")
  }

  ######### judge the format for Pagwas_list
  Pagwas_list <- lapply(Pagwas_list, function(Pagwas) {
    if (class(Pagwas) != "Seurat") stop("Pagwas_list should be constituted by
    the scPagwas result for single cell, not the celltypes result; The
    celltypes result can integrate by hand!")

    return(Pagwas)
  })
  ########## extract the "bootstrap_results" for cell types

  bootstrap_results <- lapply(Pagwas_list, function(Pagwas) {
    if ("bootstrap_results" %in% names(Pagwas@misc)) {
      return(Pagwas@misc$bootstrap_results[-1, ])
    } else {
      return(NULL)
    }
  })
  bootstrap_results <- Reduce(
    function(dtf1, dtf2) rbind(dtf1, dtf2),
    bootstrap_results
  )
  ######### remove misc
  Pagwas_list <- lapply(Pagwas_list, function(Pagwas) {
    Pagwas@misc <- list()
    return(Pagwas)
  })
  ########### merge the seruat list
  Pagwas_merge <- merge(Pagwas_list[[1]],
    y = Pagwas_list[2:length(Pagwas_list)],
    project = "merged", merge.data = TRUE
  )
  ########### get the gene_heritability_correlation
  scPagwas.gPAS.score <- data.matrix(Pagwas_merge$scPagwas.gPAS.score)
  sparse_cor <- corSparse(
    X = t(as_matrix(GetAssayData(Pagwas_merge, slot = "data", assay = "RNA"))),
    Y = scPagwas.gPAS.score
  )
  rownames(sparse_cor) <- rownames(GetAssayData(Pagwas_merge,
    slot = "data",
    assay = "RNA"
  ))
  colnames(sparse_cor) <- "gene_heritability_correlation"
  sparse_cor[is.nan(sparse_cor)] <- 0
  Pagwas_merge@misc$gene_heritability_correlation <- sparse_cor
  Pagwas_merge@misc$bootstrap_results <- bootstrap_results

  ########### get the scPagwas.TRS.Score

  scPagwas_topgenes <- names(sparse_cor[order(sparse_cor, decreasing = T), ])[1:n_topgenes]
  Pagwas_merge <- Seurat::AddModuleScore(Pagwas_merge,
    assay = "RNA",
    list(scPagwas_topgenes),
    name = c("scPagwas.TRS.Score")
  )

  ########### get the CellScalepValue

  CellScalepValue <- rankPvalue(
    datS = t(data.matrix(
      GetAssayData(Pagwas_merge, assay = "RNA")[scPagwas_topgenes, ]
    )),
    pValueMethod = "scale"
  )
  Pagwas_merge$CellScalepValue <- CellScalepValue[, "pValueHighScale"]
  Pagwas_merge$CellScaleqValue <- CellScalepValue[, "qValueHighScale"]

  return(Pagwas_merge)
}
