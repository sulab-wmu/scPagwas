#' Single_data_input
#' @description Input the Single data in seruat format
#'
#' @param Pagwas Pagwas format
#' @param Single_data Input the Single data in seruat format, Idents should
#' be the celltypes annotation. mainly used the Single_data@assays$RNA@data
#' @param Pathway_list (list,character) pathway gene sets list
#' @param nfeatures The parameter for FindVariableFeatures,
#' NULL means select all genes
#' @param min_clustercells Threshold for total cells fo each cluster.
#' @param assay assay data of your single cell data to use,default is "RNA".
#'
#' @return Pagwas list including:
#' "Celltype_anno"    "data_mat"         "VariableFeatures" "merge_scexpr"
#' @export
#'
#' @examples
#' library(scPagwas)
#' Pagwas <- list()
#' # Start to read the single cell data
#' Single_data <- readRDS(system.file("extdata", "scRNAexample.rds",
#'   package = "scPagwas"
#' ))
#' Pagwas <- Single_data_input(
#'   Pagwas = Pagwas,
#'   assay = "RNA",
#'   Single_data = Single_data,
#'   Pathway_list = Genes_by_pathway_kegg
#' )
#' @author Chunyu Deng
#' @aliases Single_data_input
#' @keywords Single_data_input, extract the single cell matrix.
Single_data_input <- function(Pagwas,
                              Single_data,
                              nfeatures = NULL,
                              Pathway_list = NULL,
                              assay = "RNA",
                              min_clustercells = 5) {
  message("Input single cell data!")

  if (!("Seurat" %in% class(Single_data))) {
    message("Single_data is not of class Seurat!")
    return(Pagwas)
  }
  # 1.Celltype_anno
  Celltype_anno <- data.frame(
    cellnames = rownames(Single_data@meta.data),
    annotation = as.vector(SeuratObject::Idents(Single_data))
  )
  rownames(Celltype_anno) <- Celltype_anno$cellnames

  Afterre_cell_types <- table(Celltype_anno$annotation) > min_clustercells
  Afterre_cell_types <- names(Afterre_cell_types)[Afterre_cell_types]
  message(
    length(Afterre_cell_types),
    " cell types are remain, after filter!"
  )

  Celltype_anno <- Celltype_anno[Celltype_anno$annotation %in%
    Afterre_cell_types, ]
  Single_data <- Single_data[, Celltype_anno$cellnames]

  Pagwas$Celltype_anno <- Celltype_anno

  Pagwas$data_mat <- GetAssayData(Single_data, slot = "data", assay = assay)

  merge_scexpr <- Seurat::AverageExpression(Single_data, assays = assay)[[assay]]

  # 5.VariableFeatures

    Pagwas$VariableFeatures <- rownames(Pagwas$data_mat)
  rm(Single_data)
  pagene <- intersect(
    unique(unlist(Pathway_list)),
    rownames(Pagwas$data_mat)
  )
  if (length(pagene) < 100) {
    stop("There are little match between rownames of Single_data and
         pathway genes!")
  }

  Pagwas$VariableFeatures <- intersect(Pagwas$VariableFeatures, pagene)
  if (ncol(merge_scexpr) == 1) {
    a <- colnames(merge_scexpr)
    Pagwas$merge_scexpr <- data.matrix(merge_scexpr[Pagwas$VariableFeatures, ])
    rownames(Pagwas$merge_scexpr) <- Pagwas$VariableFeatures
    colnames(Pagwas$merge_scexpr) <- a
  } else {
    Pagwas$merge_scexpr <- merge_scexpr[Pagwas$VariableFeatures, ]
  }

  return(Pagwas)
}

#' @title scCount_data_input
#' @description Input the scCount data in dataframe or matrix, the rownames is gene names, the colnames is cell names.
#' @param Pagwas Pagwas format
#' @param count_data (data.frame or matrix) count data of your single cell data to use, the rownames is gene names, the colnames is cell names.
#' @param meta_data (data.frame) meta data of your single cell data to use, the first column is cell names. the cell names should be the same as the first row of count_data. It should include the "cell_type" column.
#' @param Pathway_list (list,character) pathway gene sets list
#' @param min_clustercells Threshold for total cells fo each cluster.
#' 
#' @return Pagwas list including:
#' "Celltype_anno"    "data_mat"  "merge_scexpr"
#' @export
#' @examples
#' library(scPagwas)
#'
#' @author Chunyu Deng
#' @aliases scCount_data_input
#' @keywords scCount_data_input, extract the single cell matrix.

scCount_data_input <- function(Pagwas,count_data,meta_data,Pathway_list,min_clustercells){
  #判断输入的数据是否符合要求
  if(!is.data.frame(count_data) & !is.matrix(count_data)){
    stop("count_data should be data.frame or matrix!")
  }
  if(!is.data.frame(meta_data)){
    stop("meta_data should be data.frame!")
  }
  #判断count_data的列数是否和meta_data的行数一致
  if(ncol(count_data) != nrow(meta_data)){
    stop("The number of columns of count_data should be the same as the number of rows of meta_data!")
  }
  #判断count_data的列名是否和meta_data的行名一致
  if(!all(colnames(count_data) == rownames(meta_data))){
    stop("The column names of count_data should be the same as the row names of meta_data!")
  }
  #判断count_data的行名是否全部不在Pathway_list中，如果全部不在，则报错，说明行名不是基因名
  if(all(!rownames(count_data) %in% unlist(Pathway_list))){
    stop("The row names of count_data should be gene names!")
  }
  #判断meta_data的第一列的列名是否有cell_type，如果没有，则报错
  if(!"cell_type" %in% colnames(meta_data)){
    stop("The first column of meta_data should be cell_type!")
  }
  #判断meta_data的第一列的cell_type是否有NA，如果有，则报错
  if(any(is.na(meta_data$cell_type))){
    stop("The cell_type column of meta_data should not have NA!")
  }
  #判断meta_data的第一列的cell_type是否为character，如果不是，则报错
  if(!is.character(meta_data$cell_type)){
    stop("The cell_type column of meta_data should be character!")
  }
  Celltype_anno <- data.frame( cellnames = colnames(count_data),
                               annotation = meta_data$cell_type)
  rownames(Celltype_anno) <- Celltype_anno$cellnames
  Afterre_cell_types <- table(Celltype_anno$annotation) > min_clustercells
  Afterre_cell_types <- names(Afterre_cell_types)[Afterre_cell_types]
  message(
    length(Afterre_cell_types),
    " cell types are remain, after filter!"
  )

  Celltype_anno <- Celltype_anno[Celltype_anno$annotation %in%
    Afterre_cell_types, ]
  count_data <- count_data[, Celltype_anno$cellnames]

  Pagwas$Celltype_anno <- Celltype_anno
  Pagwas$data_mat <- count_data
  rm(count_data)
  #基于Pagwas$Celltype_anno$annotation计算每组细胞类型的在Pagwas$data_mat相应的列的均值
  #并将其作为Pagwas$merge_scexpr的列名
  merge_scexpr <- lapply(
    unique(Pagwas$Celltype_anno$annotation),
    function(x) {
      rowMeans(Pagwas$data_mat[, Pagwas$Celltype_anno$annotation == x])
    }
  )
  merge_scexpr <- do.call(cbind, merge_scexpr)
  colnames(merge_scexpr) <- unique(Pagwas$Celltype_anno$annotation)
  pagene <- intersect(
    unique(unlist(Pathway_list)),
    rownames(Pagwas$data_mat)
  )
  if (length(pagene) < 100) {
    stop("There are little match between rownames of Single_data and
         pathway genes!")
  }

  Pagwas$VariableFeatures <- intersect(Pagwas$VariableFeatures, pagene)
  if (ncol(merge_scexpr) == 1) {
    a <- colnames(merge_scexpr)
    Pagwas$merge_scexpr <- data.matrix(merge_scexpr[Pagwas$VariableFeatures, ])
    rownames(Pagwas$merge_scexpr) <- Pagwas$VariableFeatures
    colnames(Pagwas$merge_scexpr) <- a
  } else {
    Pagwas$merge_scexpr <- merge_scexpr[Pagwas$VariableFeatures, ]
  }
  return(Pagwas)
}

