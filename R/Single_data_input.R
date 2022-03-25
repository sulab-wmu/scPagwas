

#' Single_data_input
#' @description Input the Single data in seruat format
#'
#' @param Pagwas Pagwas format
#' @param Single_data Input the Single data in seruat format, Idents should be
#' the celltypes annotation.
#' @param FilterSingleCell whther to filter the single cell data,if you
#' filter it before,choose FALSE, otherwise set TRUE.
#' @param Pathway_list (list,character) pathway gene sets list
#' @param nfeatures The parameter for FindVariableFeatures,
#' NULL means select all genes
#' @param min.lib.size Threshold for single data library.
#' @param min.reads Threshold for total reads fo each cells.
#' @param min.detected Threshold for total cells fo each gene.
#' @param min_clustercells Threshold for total cells fo each cluster.
#'
#' @return Pagwas
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(scRNAexample)
#' # Pagwas is better getting from GWAS_summary_input(),or NULL is right,
#' Pagwas <- Single_data_input(Pagwas = Pagwas, Single_data = scRNAexample)
Single_data_input <- function(Pagwas,
                              Single_data,
                              nfeatures = NULL,
                              FilterSingleCell=FALSE,
                              Pathway_list=NULL,
                              min.lib.size = 1000,
                              min.reads = 10,
                              min.detected = 5,
                              min_clustercells = 5) {
  message("Input single cell data!")

  if (!("Seurat" %in% class(Single_data))) {
    message("Single_data is not of class Seurat!")
    return(Pagwas)
  }

  raw_data_mat <- FBM(nrow(Single_data), ncol(Single_data))
  raw_data_mat[] <- as(Seurat::GetAssayData(object = Single_data, slot = "data"),"matrix")
  dim_raw_data_mat<-list(row=rownames(Single_data),col=colnames(Single_data))

  Pagwas$raw_data_mat<-raw_data_mat
  Pagwas$dim_raw_data_mat<-dim_raw_data_mat
  Celltype_anno <- data.frame(cellnames = rownames(Single_data@meta.data),
                              annotation = as.vector(SeuratObject::Idents(Single_data)))

  rownames(Celltype_anno) <- Celltype_anno$cellnames

  if (!is.null(nfeatures)) {
    if (nfeatures < nrow(Single_data)) {
      Single_data <- Seurat::FindVariableFeatures(Single_data,
        selection.method = "vst",
        nfeatures = nfeatures
      )
      Pagwas$VariableFeatures <- Seurat::VariableFeatures(Single_data)
      # Single_data <- Single_data[,]
    } else if (nfeatures == nrow(Single_data)) {
      Pagwas$VariableFeatures <- rownames(Single_data)
    } else {
      stop("Error: nfeatures is too big")
    }
  } else {
    Pagwas$VariableFeatures <- rownames(Single_data)
  }

  Afterre_cell_types <- table(Celltype_anno$annotation) > min_clustercells
  Afterre_cell_types <- names(Afterre_cell_types)[Afterre_cell_types]
  message(length(Afterre_cell_types), " cell types are remain, after filter!")
  Celltype_anno <- Celltype_anno[Celltype_anno$annotation %in% Afterre_cell_types, ]
  Single_data <-Single_data[,Idents(Single_data) %in% Afterre_cell_types]


  pagene<- intersect(unique(unlist(Pathway_list)),rownames(Single_data))
  if(length(pagene)<100){
    stop("There are little match between rownames of Single_data and pathway genes!")
  }
  merge_scexpr <- Seurat::AverageExpression(Single_data)
  merge_scexpr <- merge_scexpr[["RNA"]][pagene,]

  Pagwas$VariableFeatures <- intersect(Pagwas$VariableFeatures,pagene)

  data_mat <- as_FBM(as(Seurat::GetAssayData(object = Single_data[pagene,], slot = "data"),"matrix"))
  dim_data_mat<-list(row=pagene,col=colnames(Single_data))

  #message("*merge_scexpr")
  Pagwas$merge_scexpr<-merge_scexpr
  #dim_raw_data_mat<-list(row=rownames(Single_data),col=colnames(Single_data))

  Pagwas$data_mat<-data_mat
  Pagwas$dim_data_mat<-dim_data_mat
  Pagwas$Single_data<-Single_data

  return(Pagwas)
}
