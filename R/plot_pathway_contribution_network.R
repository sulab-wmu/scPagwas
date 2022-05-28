
#' plot_pathway_contribution_network
#'
#' @description plot the network based on pathway contribution
#' oaqc
#'
#' Thanks to the SCOPfunctions package
#' https://github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions/blob/main/R/plot.R
#'
#' @param mat_datExpr data matrix for single cell data
#' @param vec_gene_contribution vector of pathway contribution score from pagwas
#' @param vec_pathways_highlight vector of pathway highlight from your choose.
#' @param n_max_genes number of max pathways for showing.
#' @param igraph_algorithm defalt is drl,method for igragh,details see igraph packages.
#' @param fontface_labels fontface for labels
#' @param color_edge edge color.
#' @param fontSize_label_lg fontSize for lg label
#' @param fontSize_legend_lg fontSize for lg legend
#' @param fontSize_legend_xlg fontSize for xlg legend
#' @param edge_thickness thickness of edge
#' @param do_plot Whether to plot, logical
#' @param figurenames The filename and address of the output plot
#' @param width figure width
#' @param height figure height
#'
#' @export
#' @examples
#' plot_pathway_contribution_network(
#'   mat_datExpr = pca_cell_df,
#'   vec_gene_contribution = Pagwas$Pathway_block_heritability,
#'   vec_pathways_highlight = names(sort(Pagwas$Pathway_block_heritability, decreasing = T)[1:5]),
#'   n_max_genes = 20,
#'   igraph_algorithm = "drl",
#'   fontface_labels = "bold.italic",
#'   color_edge = "#9D9D9D",
#'   fontSize_label_lg = 4,
#'   fontSize_legend_lg = 4,
#'   fontSize_legend_xlg = 4,
#'   edge_thickness = 1
#' )
plot_pathway_contribution_network <- function(mat_datExpr = Pagwas$data_mat,
                                              gene_contribution = Pagwas$gene_heritability_correlation,
                                              Celltype_anno = Pagwas$Celltype_anno,
                                              Pathway_list = Pagwas$Pathway_list,
                                              celltypes,
                                              pathwat_select,
                                              pathways_thre = 0.05,
                                              adjacency.power = 4,
                                              n_max_genes = 200,
                                              text_topn = 10,
                                              igraph_algorithm = "drl",
                                              h_method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                              fontface_labels = "bold.italic",
                                              color_edge = "lightblue",
                                              color_group = pal_d3(alpha = 0.5)(10),
                                              node_alpha = 0.6,
                                              # edge_filter,
                                              fontSize_label_lg = 2,
                                              fontSize_legend_lg = 3,
                                              fontSize_legend_xlg = 3,
                                              limits = c(0, 14),
                                              range = c(1, 5),
                                              edge_thickness = 2,
                                              do_plot = F,
                                              figurenames = NULL,
                                              width = 7,
                                              height = 7) {
  ##### debug


  # pathwat_select=rownames(Pagwas$scPathways_rankPvalue)[order(Pagwas$scPathways_rankPvalue[,celltypes],decreasing = F)[1]];
  # celltypes="OPC";
  # mat_datExpr=Single_data@assays$RNA@data;
  # gene_contribution=Pagwas$gene_heritability_correlation;
  # Celltype_anno=Pagwas$Celltype_anno;
  # Pathway_list=Pagwas$Pathway_list;
  # pathwat_select=pathway;
  # celltypes="OPC";
  # pathways_thre=0.05;
  # adjacency.power=4;
  # n_max_genes=200;
  # n_group=4;
  # text_topn=10;
  # igraph_algorithm = "drl";
  # h_method =c("complete");
  # fontface_labels="bold.italic";
  # color_edge = "lightblue";
  # color_group=pal_d3(alpha =0.5)(10);
  # node_alpha=0.6;
  # edge_filter=0.001;
  # fontSize_label_lg=2;
  # fontSize_legend_lg=3;
  # fontSize_legend_xlg=3;
  # limits=c(0,14);
  # range=c(1,5);
  # edge_thickness = 2;
  # do_plot=F;
  # figurenames=NULL;
  # width=7;
  # height=7;



  tt <- Celltype_anno$annotation == celltypes
  mat_datExpr_cell <- mat_datExpr[, tt]
  # gene_contribution<-gene_contribution[tt]
  # names(gene_contribution)<- rownames(vec_gene_contribution)
  # pathway_select<-
  #

  gene_contribution <- gene_contribution[, 1]
  # names(gene_contribution)<-rownames(vec_gene_contribution)

  gene_contribution <- gene_contribution[unique(unlist(Pathway_list[pathwat_select]))]

  if (length(gene_contribution) > n_max_genes) {
    gene_contribution <- gene_contribution[order(gene_contribution, decreasing = T)[1:n_max_genes]]
  }
  mat_datExpr_cell <- mat_datExpr_cell[names(gene_contribution), ]
  # distance = dist(data.matrix(mat_datExpr_cell))

  pas <- names(sort(gene_contribution, decreasing = T))[text_topn]
  # tf<-rep(F,length(gene_contribution))
  # tf[which(names(gene_contribution)%in% pas )]<-T

  # mydata.hclust = hclust(distance,method = h_method)
  # plot(mydata.hclust,hang=-1)
  # member = as.character(cutree(mydata.hclust,k = n_group))

  # names(member)<-colnames(mat_datExpr_cell)
  mat_datExpr_cell <- t(data.matrix(mat_datExpr_cell))
  mat_adj <- WGCNA::adjacency(mat_datExpr_cell,
    power = adjacency.power,
    corFnc = "cor",
    corOptions = list(use = "p"),
    type = "signed hybrid"
  )

  df_hub.data <- cbind.data.frame(pathway = names(gene_contribution), contribution = gene_contribution)
  # df_hub.data$color<- member
  mat_adj[is.na(mat_adj)] <- 0
  mat_adj[is.nan(mat_adj)] <- 0
  mat_adj[which(mat_adj < quantile(mat_adj)[4])] <- 0
  # mat_adj[ mat_adj< edge_filter]<-0
  mat_adj <- mat_adj[rowSums(mat_adj, na.rm = TRUE) != 0, colSums(mat_adj, na.rm = TRUE) != 0]

  df_hub.data <- df_hub.data[rownames(mat_adj), ]
  a <- df_hub.data$contribution * edge_thickness
  a[a > limits[2]] <- limits[2]

  names(a) <- rownames(mat_adj)

  mat_adj %>%
    igraph::graph.adjacency(mode = "undirected", weighted = T, diag = FALSE) %>%
    tidygraph::as_tbl_graph() %>%
    igraph::upgrade_graph() %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(contribution = a, color = df_hub.data$color) %>%
    tidygraph::activate(edges) %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(!tidygraph::node_is_isolated()) -> hub.plot
  hub.plot <- igraph::simplify(hub.plot)
  V(hub.plot)$grp <- rep(celltypes, length(V(hub.plot)))

  # if(length(E(hub.plot))>10000){
  #  bb <-graphlayouts::layout_as_backbone(hub.plot, keep = 0.1)
  # }else(
  bb <- graphlayouts::layout_as_backbone(hub.plot, keep = 0.4)
  # )

  # E(hub.plot)$col <- F
  # E(hub.plot)$col[bb$backbone] <- T

  module.plot <- ggraph::ggraph(hub.plot, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
    ggraph::geom_edge_link0(aes(edge_width = weight), edge_colour = color_edge) +
    ggraph::geom_node_point(aes(fill = grp, size = contribution),
      alpha = node_alpha,
      shape = 21
    ) +
    ggraph::geom_node_text(aes(
      filter = contribution >= a[pas],
      size = fontSize_label_lg,
      label = name
    ), repel = T) +
    scale_edge_width(range = c(0, 1)) +
    scale_size( # breaks = c(1,2,3,4,5,6,7,8,9,10),
      limits = limits,
      range = range
      # labels = c(" \U00B1 1", " \U00B1 2", " \U00B1 3")
    ) +
    ggforce::geom_mark_hull(
      aes(x, y, group = grp, fill = grp),
      concavity = 4,
      expand = unit(2, "mm"),
      alpha = 0.15
    ) +
    theme_graph() +
    scale_color_manual(values = color_group) +
    scale_fill_manual(values = color_group, guide = "none") +
    ggtitle(paste0(celltypes, "\n", pathwat_select)) +
    # theme(legend.position = "none")
    # theme(legend.position = "bottom")
    theme(
      legend.title.align = 1,
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0, unit = "cm"),
      legend.title = element_text(size = fontSize_legend_xlg, face = "bold"),
      legend.text = element_text(size = fontSize_legend_lg),
      legend.spacing.x = unit(0, "cm"),
      # margin: top, right, bottom, and left
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank()
    ) + theme_graph()


  if (do_plot) print(module.plot)

  if (!is.null(figurenames)) {
    pdf(file = figurenames, width = width, height = height)
    print(module.plot)
    dev.off()
  }
}
