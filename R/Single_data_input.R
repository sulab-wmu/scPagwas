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
  if (!is.null(nfeatures)) {
    if (nfeatures < nrow(Single_data)) {
      Single_data <- Seurat::FindVariableFeatures(Single_data,
        assay = assay,
        selection.method = "vst",
        nfeatures = nfeatures
      )
      Pagwas$VariableFeatures <- Seurat::VariableFeatures(Single_data)
    } else if (nfeatures == nrow(Single_data)) {
      Pagwas$VariableFeatures <- rownames(Pagwas$data_mat)
    } else {
      stop("Error: nfeatures is too big")
    }
  } else {
    Pagwas$VariableFeatures <- rownames(Pagwas$data_mat)
  }
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
