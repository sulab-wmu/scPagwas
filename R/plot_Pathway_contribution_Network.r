#'
#'
#'
#' Thanks to the SCOPfunctions package "https://github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions"
#'
#'

plot_network = function(
  mat_datExpr,
  vec_geneImportance,
  vec_genes_highlight=c(),
  n_max_genes=50,
  igraph_algorithm = "drl",
  fontface_labels="bold.italic",
  color_edge = "grey70",
  fontSize_label_lg=1,
  fontSize_legend_lg=1,
  fontSize_legend_xlg=1,
  edge_thickness = 1,
  do_plot=TRUE
) {

  # require("WGCNA")
  # require("dplyr")
  # require("tidygraph")
  # require("ggraph")
  # require("igraph")

  # sort
  vec_geneImportance = sort(vec_geneImportance, decreasing = T)

  # take take hub genes and genes to highlight (if any)
  vec_logical = names(vec_geneImportance) %in% c(names(vec_geneImportance)[1:min(n_max_genes, length(vec_geneImportance))],  vec_genes_highlight)
  vec_geneImportance = vec_geneImportance[vec_logical]

  mat_datExpr <- mat_datExpr[names(vec_geneImportance),]

  # Compute network adjacency
  mat_adj <- WGCNA::adjacency(t(mat_datExpr),
                              power = 1,
                              corFnc = "cor",
                              corOptions = list(use = "p"),
                              type = "signed hybrid")

  df_hub.data <- cbind.data.frame(gene = names(vec_geneImportance), importance = vec_geneImportance)

  mat_adj %>%
    igraph::graph.adjacency(mode = "undirected", weighted = T, diag = FALSE) %>%
    tidygraph::as_tbl_graph() %>% igraph::upgrade_graph() %>% tidygraph::activate(nodes) %>%
    dplyr::mutate(importance = df_hub.data$importance*edge_thickness) %>%
    tidygraph::activate(edges) %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(!tidygraph::node_is_isolated()) -> hub.plot
  module.plot <- ggraph::ggraph(hub.plot,
                        layout = "igraph",
                        algorithm = igraph_algorithm,
                        #layout = "stress"
  ) +
    ggraph::geom_edge_link(color=color_edge, show.legend=F, aes(alpha=weight)) +
    ggraph::geom_node_point(aes(size = weight), fill = "#73BCC9", shape=21, alpha=0.8) +
    # scale_size(breaks = c(1,2,3),
    #            limits = c(0,4),
    #            range = c(0, 8),
    #            labels = c(" \U00B1 1", " \U00B1 2", " \U00B1 3")
    #            ) +
    ggraph::geom_node_text(
      aes(label = as.character(name)),
      fontface="bold.italic",
      size=fontSize_label_lg,
      repel = T,
      # incredibly, there is a mismatch in how the hyphon "-" is encoded in the genes..
      color = c("red", "black") ) +

    guides(
      size = guide_legend(override.aes = list(
        #size=c(2,4,6)),
      ),
      keywidth = 0.8,
      keyheight = 0.8, order = 1)) +
    #title = expression(bold(paste("Log"[2],~
    #"fold-change"))))) +
    ggraph::theme_graph(base_family = 'Helvetica') +
    theme(
      legend.title.align=0.5,
      legend.position = "top",
      legend.margin = margin(0,0,0,0, unit="cm"),
      legend.title = element_text(size=fontSize_legend_xlg, face="bold"),
      legend.text = element_text(size=fontSize_legend_lg, face="bold"),
      legend.spacing.x = unit(0, "cm"),

      axis.text =element_blank(),
      axis.ticks = element_blank(),
      axis.line=element_blank(),
      axis.title = element_blank(),
      # margin: top, right, bottom, and left
      plot.margin = unit(c(0, 0, 0, 0), "cm"),

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),

      plot.background = element_blank(),
      panel.grid = element_blank()
    ) +

    coord_cartesian(clip="off")

  if (do_plot) print(module.plot)
  return(p)

}

