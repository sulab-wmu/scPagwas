#'
#' @title Main wrapper functions for scPagwas in version 2.0
#' @description Main Pagwas wrapper functions in version 2.0.
#' @details The entry point for Pagwas analysis. Including the data input...
#' Including the data input functions and the main progress functions;
#' It can also output the
#' running log and parameter log for scPagwas, and construct the folder
#' for output.
#'
#' @param Pagwas = NULL: This parameter is typically not required and does not need any input data. When seurat_return = FALSE, all intermediate data is stored in the "Pagwas" list and returned as the result. This result can be inherited and used as input for subsequent calculations. In certain scenarios, such as when performing two computations with the same single-cell data input but different GWAS data inputs, the result list obtained from the first computation can be used as the "Pagwas" parameter input for the second computation. This allows skipping the single-cell calculations, significantly expediting the process. However, when seurat_return = TRUE, the returned result cannot be manipulated in this manner, as it is the final Seurat result with many intermediate data removed.
#' @param gwas_data (data.frame)GWAS summary data; It must have some
#' colmuns such as:
#'  chrom|    pos    |   rsid    |   se  |  beta |  maf
#'     6 | 119968580 | rs1159767 | 0.032 | 0.019 |0.5275
#'    10 | 130566523 |  rs559109 | 0.033 | 0.045 |0.4047
#'     5 | 133328825 | rs6893145 | 0.048 | 0.144 |0.1222
#'     7 | 146652932 | rs13228798| 0.035 | 0.003 | 0.3211
#'
#' @param output.prefix = "Test": This parameter sets the prefix for the output result files.
#' @param output.dirs = "scPagwastest_output": This parameter specifies the directory for the output result files.
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
#' @param remove_outlier (logical)Whether to remove the outlier for
#' scPagwas score.
#' @param iters_celltype (integr)number of bootstrap iterations for celltype
#' @param iters_singlecell (integr)number of bootstrap iterations for singlecell；
#' The parameter "iters_singlecell" is used to calculate the significance p-value for individual cells. However, we have observed that this step requires a significant amount of computational memory. Therefore, we do not recommend selecting a large value for this parameter initially.If you do not want to waste time calculating the p-value, you can choose to set it as 0.
#' @param seurat_return (logical) Whether return the seurat format result,
#' if not,will return a list result;
#' @param singlecell (logical)Whether to produce the singlecell result;
#' @param celltype (logical)Whether to produce the celltypes result;
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
#'   {Random_Correct_BG_p}{: CellpValue for each cells;}
#'   {Random_Correct_BG_adjp}{: fdr for each cells, adjust p value.}
#'   {Random_Correct_BG_z}{: z score for eahc cells.}
#'   \item{misc: element in result,\code{Pagwas@misc }}{
#'   {Pathway_list:}{a list for pathway gene list intersecting with single
#'   cell data}
#'   {pca_cell_df:}{ a data frame for pathway pca result for each celltype.}
#'   {lm_results:}{ the regression result for each cell.}
#'   {PCC:}{
#'   heritability correlation value for each gene;In the previous version, we referred to it as Pearson correlation coefficients.}
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
#'   {*_singlecell_scPagwas_score.Result.csv:}{ The final result for
#'   lm and top gene score}
#'   {*_celltypes_bootstrap_results.csv:}{The bootstrap data frame
#'   results for celltypes including bootstrap pvalue and confidence
#'   interval}
#'   {*_gene_PCC.csv:}{ heritability correlation
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
#'   {PCC:}{ The total number of
#'   matrix vector products carried out;In the previous version, we referred to it as Pearson correlation coefficients}
#'   {Pathway_ctlm_results:}{ The total number of matrix vector products
#'   carried out}
#'   {lm_results:}{ The total number of matrix vector products carried out}
#'   {Pathway_ct_results:}{ The total number of matrix vector products
#'   carried out}
#' }
#' @export
#'
#' @author Chunyu Deng
#' @aliases scPagwas_main2
#' @keywords scPagwas_main2, wrapper of scPagwas functions in version 2.0.
scPagwas_main2 <- function(Pagwas = NULL,
                          gwas_data = NULL,
                          output.prefix = "Test",
                          output.dirs = "scPagwastest_output",
                          block_annotation = block_annotation,
                          Single_data = NULL,
                          assay = "RNA",
                          Pathway_list = Genes_by_pathway_kegg,
                          chrom_ld = chrom_ld,
                          run_split=FALSE,
                          n.cores=1,
                          marg = 10000,
                          maf_filter = 0.01,
                          min_clustercells = 10,
                          min.pathway.size = 5,
                          max.pathway.size = 300,
                          iters_celltype = 200,
                          iters_singlecell = 100,
                          n_topgenes = 1000,
                          singlecell = TRUE,
                          celltype = TRUE,
                          seurat_return = TRUE,
                          remove_outlier = TRUE) {


  #Pagwas =NULL
  #gwas_data = system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas")
  #Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas")
  # output.prefix=""
  # output.dirs="test"
  # #load("D:/OneDrive/GWAS_Multiomics/Manuscripts/Revise_comments/pathway_add/reduce_genes.by.gobp.pathway.RData")
  # Pathway_list=Genes_by_pathway_kegg
  # assay="RNA"
  # singlecell=T
  # celltype=F
  # block_annotation = block_annotation
  # chrom_ld = chrom_ld
  #  assay = "RNA"
  # run_split=FALSE
  # n.cores=1
  # marg = 10000
  # maf_filter = 0.01
  # min_clustercells = 10
  # min.pathway.size = 5
  # max.pathway.size = 300
  # iters_singlecell = 100
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
  if (dir.exists(paste0("./", output.dirs, "/temp"))) {
    temp_folder_path <- paste0("./", output.dirs, "/temp")
    files_to_remove <- list.files(temp_folder_path)
    invisible(file.remove(file.path(temp_folder_path, files_to_remove)))
    cat("bk files has been deleted")
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
    }else if(singlecell){
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

    } else if(class(gwas_data) == "data.frame"){
      gwas_data$se<-as.numeric(gwas_data$se)
      gwas_data$beta<-as.numeric(gwas_data$beta)
    }else{
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
  ## 8.scPagwas_perform_score
  #############################

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
  ## 9.scGet_PCC
  #############################
  if(!run_split){
    message(paste(utils::timestamp(quiet = T),
                  " ******* 9th: scGet_PCC function start! ********",
                  sep = ""
    ))

    tt <- Sys.time()

    Pagwas <- scGet_PCC2(Pagwas=Pagwas)

    utils::write.csv(Pagwas$PCC,
                     file = paste0(
                       "./", output.dirs, "/",
                       output.prefix,
                       "_gene_PCC.csv"
                     ),
                     quote = F
    )

    scPagwas_topgenes <- rownames(Pagwas$PCC)[order(Pagwas$PCC$weight_pcc, decreasing = T)[1:n_topgenes]]
    message("done!")

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


    Pagwas[c(
      "VariableFeatures", "merge_scexpr", "snp_gene_df",
      "rawPathway_list", "data_mat"
    )] <- NULL


    Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, scPagwas_topgenes, name ="scPagwas.TRS.Score")

    message("* Get Random Correct background pvalue for each single cell!")
    if(iters_singlecell>0){
      correct_pdf<-Get_CorrectBg_p(Single_data=Single_data,
                                   scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                                   iters_singlecell=iters_singlecell,
                                   n_topgenes=n_topgenes,
                                   scPagwas_topgenes=scPagwas_topgenes,
                                   assay=assay
      )
      Pagwas$Random_Correct_BG_pdf <- correct_pdf
      message("* Get Merged pvalue for each celltype!")
      Pagwas$Merged_celltype_pvalue<-Merge_celltype_p(single_p=correct_pdf$pooled_p,celltype=Pagwas$Celltype_anno$annotation)

      Pagwas$scPagwas.TRS.Score <- Single_data$scPagwas.TRS.Score1


      a <- data.frame(
        scPagwas.TRS.Score = Pagwas$scPagwas.TRS.Score,
        scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
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
    }else{
      Pagwas$Random_Correct_BG_pdf <-NA
      Pagwas$Merged_celltype_pvalue <-NA
      Pagwas$scPagwas.TRS.Score <- Single_data$scPagwas.TRS.Score1

      a <- data.frame(
        scPagwas.TRS.Score = Pagwas$scPagwas.TRS.Score,
        scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score)

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
        "Pathway_sclm_results",
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
      if(iters_singlecell>0){
        Single_data$Random_Correct_BG_p <- correct_pdf$pooled_p
        Single_data$Random_Correct_BG_adjp <- correct_pdf$adj_p
        Single_data$Random_Correct_BG_z <- correct_pdf$pooled_z
      }
      Pagwas[c(
        "scPagwas.gPAS.score", "Celltype_anno",
        "Pathway_single_results", "pca_scCell_mat","Random_Correct_BG_pdf"
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

  if (file.exists(paste0("./", output.dirs, "/scPagwas_cache"))) {
    unlink(paste0("./", output.dirs, "/scPagwas_cache"),recursive = TRUE)
  }
}

