
#' Pagwas_result_integtate
#'
#' @param Pagwas
#' @param seruat_return
#' @param output.dirs
#' @param output.prefix
#' @param n_topgenes
#' @param assay
#' @param Single_data
#'
#' @return
#' @export
#'
#' @examples
Pagwas_result_integtate<-function(Pagwas,
                                  seruat_return,
                                  output.dirs,
                                  Single_data,
                                  output.prefix,
                                  n_topgenes,
                                  assay){
#################remove some unused parameter
Pagwas[c(
  "VariableFeatures", "merge_scexpr",
  "rawPathway_list"
)] <- NULL
if (!seruat_return) {
  return(Pagwas)
}

##########output
write.table(Pagwas$pca_scCell_mat, file = paste0("./", output.dirs, "/", output.prefix, "_pca_scCell_mat.txt"), quote = F)
write.table(Pagwas$pca_cell_df, file = paste0("./", output.dirs, "/", output.prefix, "_pca_celltypes_mat.txt"), quote = F)
write.csv(Pagwas$bootstrap_results, file = paste0("./", output.dirs, "/", output.prefix, "_cellytpes_bootstrap_results.csv"), quote = F)
write.csv(Pagwas$allsnp_gene_heritability_correlation, file = paste0("./", output.dirs, "/", output.prefix, "_gene_heritability_correlation.csv"), quote = F)

#########seruat result
scPagwas_pca <- SeuratObject::CreateAssayObject(data = Pagwas$pca_scCell_mat)
Single_data[["scPagwasPaPca"]] <- scPagwas_pca

rm(scPagwas_pca)
Pagwas[c("pca_scCell_mat","data_mat","snp_gene_df")] <- NULL

#########scPagwaslm_topgenes

scPagwaslm_topgenes <- names(Pagwas$allsnp_gene_heritability_correlation[order(Pagwas$allsnp_gene_heritability_correlation, decreasing = T), ])[1:n_topgenes]

Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(scPagwaslm_topgenes),
                                      name = c("scPagwas.lmtopgenes.Score"))
Single_data$scPagwas.lmtopgenes.Score1 <- scPagwas_score_filter(scPagwas_score = Single_data$scPagwas.lmtopgenes.Score1)

Single_data$sclm_score<-Pagwas$sclm_allsnpscore

write.csv(Single_data@meta.data[, c(
  "sclm_score",
  "scPagwas.lmtopgenes.Score1"
)], file = paste0("./", output.dirs, "/", output.prefix, "_singlecell_scPagwas_score_pvalue.Result.csv"), quote = F)

Pagwas[c("Celltype_anno","sclm_allsnpscore")] <- NULL

Single_data@misc <- Pagwas
return(Single_data)
}
