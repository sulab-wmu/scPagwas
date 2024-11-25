
#' Pathway_annotation_input
#' @description Deal with the pathway genes with snp blocks
#' Initial the pathway blocks for mapping the pathway to snp.
#'
#' @param Pagwas Pagwas format data list inherit from other functions.
#' @param block_annotation Start and end points for block traits,
#' usually genes.
#'
#' @return list for pathway blocks
#' @export
#'
#' @author Chunyu Deng
#' @aliases Pathway_annotation_input
#' @keywords Pathway_annotation_input, Initial the pathway blocks
#' for mapping the pathway to snp.
Pathway_annotation_input <- function(Pagwas,
                                     block_annotation) {
  if (class(Pagwas$Pathway_list) != "list") {
    stop("require list of Pathway_list")
  }

  if (class(block_annotation) != "data.frame") {
    stop("require data frame of annotations")
  }

  if (all(!(c("chrom", "start", "end", "label") %in%
    colnames(block_annotation)))) {
    stop("require chrom, start, and end, cols")
  }

  if (any(duplicated(block_annotation$label))) {
    block_annotation <- block_annotation[!duplicated(
      block_annotation$label
    ), ]
  }

  proper.gene.names <- intersect(
    Pagwas$VariableFeatures,
    Pagwas$snp_gene_df$label
  )

  if (length(intersect(unlist(Pagwas$Pathway_list), proper.gene.names)) < 1) {
    stop("no match for Pathway gene and VariableFeatures")
  }

  chrom <- start <- NULL
  # sort and add a random partition label for bootstrap
  Pathway_list <- lapply(
    Pagwas$Pathway_list,
    function(Pa) intersect(Pa, proper.gene.names)
  )

  a <- names(Pathway_list)[unlist(lapply(
    Pathway_list,
    function(Pa) length(Pa)
  )) != 0]
  Pa_index <- intersect(a, rownames(Pagwas$pca_cell_df))

  Pathway_list <- Pathway_list[Pa_index]
  message("obtain the pathway block information")

  Pa_index <- unlist(lapply(names(Pathway_list), function(pa) {
    paan <- sum(block_annotation$label %in% Pathway_list[[pa]])
    if(paan<10){
      return(NULL)
    }else{
      return(pa)
    }
  }))
  Pathway_list <- Pathway_list[Pa_index]
  paan_df <- lapply(names(Pathway_list), function(pa) {
    paan <- block_annotation[block_annotation$label %in% Pathway_list[[pa]], ]
    paan$pathway <- rep(pa, nrow(paan))
    return(paan)
  })

  pathway_blocks <- lapply(paan_df, function(pa) {
    blocks <- pa %>% dplyr::arrange(chrom, start)
    return(blocks)
  })

  names(pathway_blocks) <- names(Pathway_list)

  Pagwas$Pathway_list <- Pathway_list

  Pagwas$pathway_blocks <- pathway_blocks
  return(Pagwas)
}

