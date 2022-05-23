
#' reduce_pathway
#'
#' @param pathway_seed
#' @param pathway_list
#' @param remove_proporion
#'
#' @return
#' @export
#'
#' @examples
reduce_pathway <- function(pathway_seed = names(Genes_by_pathway_kegg)[sample(1:length(Genes_by_pathway_kegg), 20)],
                           pathway_list = Genes_by_pathway_kegg,
                           remove_proporion = 0.7) {

  # seed_gene<-unique(unlist(pathway_list[pathway_seed]))

  pre_pathways <- setdiff(names(pathway_list), pathway_seed)
  p_pa <- unlist(lapply(pre_pathways, function(x) {
    p_x <- unlist(lapply(pathway_seed, function(y) {
      p_y <- length(intersect(pathway_list[[x]], pathway_list[[y]])) / length(pathway_list[[y]])
      if (p_y > remove_proporion) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }))
    if (sum(p_x) == 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))
  pre_pathways <- pre_pathways[p_pa]
  return(pathway_list[c(pathway_seed, pre_pathways)])
}

# load("/share/pub/dengcy/GWAS_Multiomics/pagwas/data/genes.by.reactome.pathway.RData")
# reduce_genes.by.reactome.pathway<-reduce_pathway(pathway_seed=names(genes.by.reactome.pathway)[sample(1:length(genes.by.reactome.pathway),500)],
# pathway_list=genes.by.reactome.pathway,
# remove_proporion=0.6)
# save(reduce_genes.by.reactome.pathway,file = "/share/pub/dengcy/GWAS_Multiomics/pagwas/data/reduce_genes.by.reactome.pathway.RData")
