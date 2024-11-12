
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


import pandas as pd

def pathway_annotation_input(Pagwas, block_annotation):
    # Check the input types
    if not isinstance(Pagwas['Pathway_list'], dict):
        raise ValueError("Pathway_list must be a dictionary")

    if not isinstance(block_annotation, pd.DataFrame):
        raise ValueError("block_annotation must be a pandas DataFrame")

    required_columns = {'chrom', 'start', 'end', 'label'}
    if not required_columns.issubset(block_annotation.columns):
        raise ValueError("block_annotation must contain columns: 'chrom', 'start', 'end', 'label'")

    # Remove duplicated labels in block_annotation
    block_annotation = block_annotation.drop_duplicates(subset='label')

    # Intersect variable features with gene labels
    proper_gene_names = set(Pagwas['VariableFeatures']).intersection(Pagwas['snp_gene_df']['label'])

    if len(set().union(*Pagwas['Pathway_list'].values()).intersection(proper_gene_names)) < 1:
        raise ValueError("No match between Pathway genes and VariableFeatures")

    # Filter the Pathway_list to only include valid genes
    Pathway_list = {
        pa: list(set(genes).intersection(proper_gene_names))
        for pa, genes in Pagwas['Pathway_list'].items()
    }

    # Filter out empty pathways
    Pathway_list = {k: v for k, v in Pathway_list.items() if len(v) > 0}

    # Keep only pathways present in pca_cell_df's rows
    Pa_index = set(Pathway_list.keys()).intersection(Pagwas['pca_cell_df'].index)
    Pathway_list = {k: Pathway_list[k] for k in Pa_index}

    print("Obtaining pathway block information...")

    # Filter pathways based on block_annotation
    valid_pa_index = [
        pa for pa in Pathway_list if block_annotation['label'].isin(Pathway_list[pa]).sum() >= 10
    ]
    Pathway_list = {pa: Pathway_list[pa] for pa in valid_pa_index}

    # Create a DataFrame for each pathway and label it
    paan_df = []
    for pa, genes in Pathway_list.items():
        pa_block = block_annotation[block_annotation['label'].isin(genes)].copy()
        pa_block['pathway'] = pa
        paan_df.append(pa_block)

    # Sort the blocks by chromosome and start position
    pathway_blocks = {
        pa: df.sort_values(['chrom', 'start']) for pa, df in zip(Pathway_list.keys(), paan_df)
    }

    # Update Pagwas object with filtered Pathway_list and pathway_blocks
    Pagwas['Pathway_list'] = Pathway_list
    Pagwas['pathway_blocks'] = pathway_blocks

    return Pagwas
