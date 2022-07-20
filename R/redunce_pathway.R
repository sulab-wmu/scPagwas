#' reduce_pathway
#' @description Sometimes, the gene lists for pathways are duplicated,
#' it needs a lot of times to running scPagwas.
#' We can used reduce_pathway to reduce them.
#' @param pathway_seed randome select the seed pathways, sometimes can be the
#' unredundant pathway set.
#' @param pathway_list the list of pathway need to reduce.
#' @param remove_proporion the propotion of pathways after reduce.
#'
#' @return reduced pathway gene list
#' @export
#'
#' @examples
#' reduce_pathway(
#'   pathway_seed = names(Genes_by_pathway_kegg)[sample(
#'     1:length(Genes_by_pathway_kegg), 20
#'   )],
#'   pathway_list = Genes_by_pathway_kegg,
#'   remove_proporion = 0.7
#' )
#' @author Chunyu Deng
#' @aliases reduce_pathway
#' @keywords reduce_pathway, reduce the replicated pathway gene list.

reduce_pathway <- function(pathway_seed,
                           pathway_list,
                           remove_proporion = 0.7) {
  pre_pathways <- setdiff(names(pathway_list), pathway_seed)
  p_pa <- unlist(lapply(pre_pathways, function(x) {
    p_x <- unlist(lapply(pathway_seed, function(y) {
      p_y <- length(intersect(pathway_list[[x]], pathway_list[[y]])) /
        length(pathway_list[[y]])
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
  pathway_list <- pathway_list[c(pathway_seed, pre_pathways)]
  return(pathway_list)
}
