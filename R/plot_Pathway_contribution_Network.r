
#' Thanks to the SCOPfunctions package
#' https://github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions/blob/main/R/plot.R
#'
NULL


#' plot_pathway_contribution_network
#'
#' @description plot the network based on pathway contribution
#'
#' @param mat_datExpr data matrix for single cell data
#' @param vec_pathwaycontribution vector of pathway contribution score from pagwas
#' @param vec_pathways_highlight vector of pathway highlight from your choose.
#' @param n_max_pathways number of max pathways for showing.
#' @param igraph_algorithm defalt is drl,method for igragh,details see igraph packages.
#' @param fontface_labels fontface for labels
#' @param color_edge edge color.
#' @param fontSize_label_lg fontSize for lg label
#' @param fontSize_legend_lg fontSize for lg legend
#' @param fontSize_legend_xlg fontSize for xlg legend
#' @param edge_thickness thickness of edge
#'
#' @export
#' @examples
#' plot_pathway_contribution_network(
#' mat_datExpr=pca_cell_df,
#' vec_pathwaycontribution=Pagwas$Pathway_block_heritability,
#' vec_pathways_highlight=names(sort(Pagwas$Pathway_block_heritability,decreasing = T)[1:5]),
#' n_max_pathways=20,
#' igraph_algorithm = "drl",
#' fontface_labels="bold.italic",
#' color_edge = "#9D9D9D",
#' fontSize_label_lg=4,
#' fontSize_legend_lg=4,
#' fontSize_legend_xlg=4,
#' edge_thickness = 1
#' )

plot_pathway_contribution_network <- function(
  mat_datExpr,
  vec_pathwaycontribution,
  vec_pathways_highlight=c(),
  n_max_pathways=50,
  igraph_algorithm = "drl",
  fontface_labels="bold.italic",
  color_edge = "#9D9D9D",
  fontSize_label_lg=1,
  fontSize_legend_lg=1,
  fontSize_legend_xlg=1,
  edge_thickness = 1
) {

  # sort
  vec_pathwaycontribution = sort(vec_pathwaycontribution, decreasing = T)
  # take take hub pathways and pathways to highlight (if any)
  vec_logical = names(vec_pathwaycontribution) %in% c(names(vec_pathwaycontribution)[1:min(n_max_pathways, length(vec_pathwaycontribution))],  vec_pathways_highlight)
  vec_pathwaycontribution = vec_pathwaycontribution[vec_logical]

  mat_datExpr <- mat_datExpr[names(vec_pathwaycontribution),]

  # Compute network adjacency
  mat_adj <- WGCNA::adjacency(t(mat_datExpr),
                              power = 1,
                              corFnc = "cor",
                              corOptions = list(use = "p"),
                              type = "signed hybrid")

  df_hub.data <- cbind.data.frame(pathway = names(vec_pathwaycontribution), contribution = vec_pathwaycontribution)
  df_hub.data$color<-rep("none",nrow(df_hub.data))
  df_hub.data$color[which(df_hub.data$pathway %in% vec_pathways_highlight)]<-"Highlight"

  mat_adj %>%
    igraph::graph.adjacency(mode = "undirected", weighted = T, diag = FALSE) %>%
    tidygraph::as_tbl_graph() %>% igraph::upgrade_graph() %>% tidygraph::activate(nodes) %>%
    dplyr::mutate(contribution = df_hub.data$contribution*edge_thickness,color=df_hub.data$color) %>%
    tidygraph::activate(edges) %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(!tidygraph::node_is_isolated())-> hub.plot


  ##plot the netword
    module.plot <- ggraph::ggraph(hub.plot,
                        layout = "igraph",
                        algorithm = igraph_algorithm,
                        #layout = "stress"
  ) +
    ggraph::geom_edge_link(color=color_edge, show.legend=T, aes(alpha=weight)) +
      scale_fill_manual(values=c("#FF4848","#64C9CF"))+
    ggraph::geom_node_point(aes(size = contribution,fill=color), shape=21, alpha=0.8) +
    # scale_size(breaks = c(1,2,3),
    #            limits = c(0,4),
    #            range = c(0, 8),
    #            labels = c(" \U00B1 1", " \U00B1 2", " \U00B1 3")
    #            ) +
    ggraph::geom_node_text(
      aes(label = as.character(name)),
      #fontface="bold.italic",
      size=fontSize_label_lg,
      repel = T,
      colour = "black") +

    #guides(
    #  size = guide_legend(override.aes = list(
     #   #size=c(2,4,6)),
     # ),
     # keywidth = 0.8,
     # keyheight = 0.8, order = 1)) +
    #title = expression(bold(paste("Log"[2],~
    #"fold-change"))))) +
    #ggraph::theme_graph(base_family = 'Arial Narrow') +
    theme(
      legend.title.align=0.5,
      legend.position = "top",
      #legend.margin = margin(0,0,0,0, unit="cm"),
      legend.title = element_text(size=fontSize_legend_xlg, face="bold"),
      legend.text = element_text(size=fontSize_legend_lg, face="bold"),
      #legend.spacing.x = unit(0, "cm"),

      axis.text =element_blank(),
      axis.ticks = element_blank(),
      axis.line=element_blank(),
      axis.title = element_blank(),
      # margin: top, right, bottom, and left
     # plot.margin = unit(c(0, 0, 0, 0), "cm"),

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),

      plot.background = element_blank(),
      panel.grid = element_blank()
    ) #+

    #coord_cartesian(clip="off")

  #if (do_plot) print(module.plot)
  return(module.plot)

}
