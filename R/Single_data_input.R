

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
                              min_clustercells = 50) {
  message("Input single cell data!")

  if (!("Seurat" %in% class(Single_data))) {
    message("Single_data is not of class Seurat!")
    return(Pagwas)
  }
  ## get the Variable genes for cells
  #raw_data_mat <-  Seurat::GetAssayData(object = Single_data, slot = "data")
 # if(){
 #   Pagwas$ff.raw_data_mat <- ff::ff(vmode="double", dim=c(dim(Single_data@assays$RNA@data)[1], dim(Single_data@assays$RNA@data)[2]))
 #   Pagwas$ff.raw_data_mat[1:dim(Single_data@assays$RNA@data)[1],1:dim(Single_data@assays$RNA@data)[2]] <- Single_data@assays$RNA@data[1:dim(Single_data@assays$RNA@data)[1],1:dim(Single_data@assays$RNA@data)[2]]
  #  rownames(Pagwas$ff.raw_data_mat)<-rownames(Single_data@assays$RNA@data)
  #  colnames(Pagwas$ff.raw_data_mat)<-colnames(Single_data@assays$RNA@data)

  #}else{
    raw_data_mat <-  Seurat::GetAssayData(object = Single_data, slot = "data")
    #raw_data_mat <- big_copy(raw_data_mat)
    SOAR::Store(raw_data_mat)
  #}


  Celltype_anno <- data.frame(cellnames = rownames(Single_data@meta.data),
                              annotation = as.vector(SeuratObject::Idents(Single_data)))

  rownames(Celltype_anno) <- Celltype_anno$cellnames
  #Pagwas$Celltype_anno <- Celltype_anno

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

  if(FilterSingleCell){

  count <- Seurat::GetAssayData(object = Single_data, slot = "count")
  count <- as_matrix(count)
  # remove cells that don't have enough counts
  count <- count[, colSums(count > 0) > min.lib.size]
  # remove genes that don't have many reads
  count <- count[rowSums(count) > min.reads, ]
  # remove genes that are not seen in a sufficient number of cells
  count <- count[rowSums(count > 0) > min.detected, ]
  # identify for Celltype_anno
  Celltype_anno <- Celltype_anno[colnames(count), ]
  # remove celltypes that don't have enough cell

  Afterre_cell_types <- table(Celltype_anno$annotation) > min_clustercells
  Afterre_cell_types <- names(Afterre_cell_types)[Afterre_cell_types]
  message(length(Afterre_cell_types), "cell types are remain, after filter!")
  Celltype_anno <- Celltype_anno[Celltype_anno$annotation %in% Afterre_cell_types, ]
  count <- count[, Celltype_anno$annotation %in% Afterre_cell_types]

  #Pagwas$Celltype_anno <- Celltype_anno
  Single_data <- Single_data[, colnames(count)]
  rm(count)
  }

  pagene<- intersect(unique(unlist(Pathway_list)),rownames(Single_data))
  if(length(pagene)<100){
    stop("There are little match between rownames of Single_data and pathway genes!")
  }
  merge_scexpr <- Seurat::AverageExpression(Single_data)
  merge_scexpr <- merge_scexpr[["RNA"]][pagene,]
  data_mat <- Seurat::GetAssayData(object = Single_data, slot = "data")
  data_mat <- as_matrix(data_mat)
  data_mat <- data_mat[pagene,]

  #ff.data_mat <- ff::ff(vmode="double", dim=c(dim(data_mat)[1], dim(data_mat)[2]))
  #ff.data_mat[1:dim(data_mat)[1],1:dim(data_mat)[2]] <- data_mat[1:dim(data_mat)[1],1:dim(data_mat)[2]]
  #rownames(ff.data_mat)<-rownames(data_mat)
  #colnames(ff.data_mat)<-colnames(data_mat)

  message("*1)merge_scexpr")
  SOAR::Store(merge_scexpr)
  SOAR::Store(data_mat)

  Pagwas$Single_data<-Single_data

  #Pagwas$ff.data_mat<-ff.data_mat
  return(Pagwas)
}
