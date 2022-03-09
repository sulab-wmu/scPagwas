

#' Single_data_input
#' Input the Single data in seruat format
#' @param Pagwas Pagwas format
#' @param Single_data Input the Single data in seruat format, Idents should be the celltypes annotation.
#' @param nfeatures The parameter for FindVariableFeatures, NULL means select all genes
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

  Celltype_anno <- data.frame(cellnames = rownames(Single_data@meta.data), annotation = as.vector(Idents(Single_data)))

  #
  Celltype_anno$annotation <- str_replace_all(Celltype_anno$annotation, "-", ".")
  Celltype_anno$annotation <- str_replace_all(Celltype_anno$annotation, " ", ".")
  Celltype_anno$annotation <- str_replace_all(Celltype_anno$annotation, "\\+", ".")

  rownames(Celltype_anno) <- Celltype_anno$cellnames
  # VSingle_data <- Single_data
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

  count <- Seurat::GetAssayData(object = Single_data, slot = "count")

  # remove cells that don't have enough counts
  count <- count[, colSums(count > 0) > min.lib.size]

  # remove genes that don't have many reads
  count <- count[rowSums(count) > min.reads, ]

  # remove genes that are not seen in a sufficient number of cells
  count <- count[rowSums(count > 0) > min.detected, ]

  # identify for Celltype_anno
  Celltype_anno <- Celltype_anno[colnames(count), ]

  # remove celltypes that don't have enough cell
  celltypes <- unique(Celltype_anno$annotation)
  Afterre_cell_types <- table(Celltype_anno$annotation) > min_clustercells
  Afterre_cell_types <- names(Afterre_cell_types)[Afterre_cell_types]
  message(Afterre_cell_types, " are remain, after filter!")

  Celltype_anno <- Celltype_anno[Celltype_anno$annotation %in% Afterre_cell_types, ]

  count <- count[, Celltype_anno$annotation %in% Afterre_cell_types]
  # Single_data <- GetAssayData(object =Single_data, slot = "data")

  Pagwas$Celltype_anno <- Celltype_anno
  Pagwas$Single_data <- Single_data[, colnames(count)]


  Pagwas$merge_scexpr <- Seurat::AverageExpression(Pagwas$Single_data)[[1]]
  return(Pagwas)
}
