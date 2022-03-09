
#' Pathway_annotation_input
#' @description Deal with the pathway genes with snp blocks
#' @param Pagwas Pagwas format data list inherit from other functions.
#' @param n.cores Parallel cores,default is 1. use detectCores() to check the cores in computer.
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(Genes_by_pathway.kegg)
#' #the Pagwas should be after Single_data_input() and GWAS_summary_input(),Pathway_pcascore_run()
#' Pagwas<-Pathway_annotation_input(Pagwas=Pagwas)
#'

Pathway_annotation_input <- function(Pagwas,n.cores=1) {

  message('adding block annotations')

  block_annotation <- Pagwas$block_annotation

  if (class(Pagwas$Pathway_list) != 'list') { stop('require list of Pathway_list') }

  if (class(block_annotation) != 'data.frame') { stop('require data frame of annotations') }

  if (all(!(c('chrom', 'start', 'end', 'label') %in% colnames(block_annotation)))) {stop('require chrom, start, and end, cols')}

  if (any(duplicated(block_annotation$label))) { block_annotation <- block_annotation[!duplicated(block_annotation$label),]}

  proper.gene.names<- intersect(Pagwas$VariableFeatures,Pagwas$snp_gene_df$label)

  if ( length(intersect(unlist(Pagwas$Pathway_list),proper.gene.names))<1) { stop('no match for Pathway gene and VariableFeatures') }

  # sort and add a random partition label for bootstrap
  Pathway_list <- lapply(Pagwas$Pathway_list,function(Pa) intersect(Pa, proper.gene.names))

  Pa_index<- names(Pathway_list)[unlist(lapply(Pathway_list,function(Pa) length(Pa)))!=0] %>% intersect(.,rownames(Pagwas$pca_cell_df))

  Pathway_list<-Pathway_list[Pa_index]

  pca_cell_df <- Pagwas$pca_cell_df[Pa_index,]
  #zScores<-zScores(Pagwas$gwas_data$pvalue,tails=2, limit=.Machine$double.xmin)
  message("obtain the pathway block information")

  paan_df<-papply(names(Pathway_list),function(pa){
            paan<-block_annotation[block_annotation$label %in% Pathway_list[[pa]],]
            paan$pathway <- rep(pa,nrow(paan))
            return(paan)
          }, n.cores = n.cores)

  Pagwas$pathway_blocks<- lapply(paan_df,function(pa){
                    blocks <- pa %>% dplyr::arrange(chrom, start)
                    # %>%
                    #dplyr::mutate(partition = cut(row_number(), breaks = 10, labels = F))
                    return(blocks)
                    })

  names(Pagwas$pathway_blocks)<-names(Pathway_list)

  Pagwas$Pathway_list <- Pathway_list

  Pagwas$block_annotation <- block_annotation

  return(Pagwas)
}

