

#' Single_data_input
#' @description Input the Single data in seruat format
#' @param Pagwas Pagwas format
#' @param Single_data Input the Single data in seruat format, Idents should be
#' the celltypes annotation.
#' @param FilterSingleCell whther to filter the single cell data,if you
#' filter it before,choose FALSE, otherwise set TRUE.
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
                              min.lib.size = 1000,
                              min.reads = 10,
                              min.detected = 5,
                              min_clustercells = 50) {
  message("Input single cell count data and cell annotation data!")

  if (!("Seurat" %in% class(Single_data))) {
    message("Single_data is not of class Seurat!")
    return(Pagwas)
  }
  ## get the Variable genes for cells

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

  # remove cells that don't have enough counts
  count <- count[, Matrix::colSums(count > 0) > min.lib.size]

  # remove genes that don't have many reads
  count <- count[Matrix::rowSums(count) > min.reads, ]

  # remove genes that are not seen in a sufficient number of cells
  count <- count[Matrix::rowSums(count > 0) > min.detected, ]

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

  merge_scexpr <- Seurat::AverageExpression(Single_data)
  merge_scexpr<-merge_scexpr[["RNA"]]

  #merge_scexpr <- mean_expr(Single_data)
  SOAR::Store(merge_scexpr)
  SOAR::Store(Single_data)
  return(Pagwas)
}
