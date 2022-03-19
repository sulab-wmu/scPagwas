#' Thanks to the SCOPfunctions package "https://github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions"
#'
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
#' plot_network(
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

plot_network = function(
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

#' Make grid of violin plots
#'
#' @description produce a n_celltype * n_genes grid of violin plots
#' https://github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions/blob/main/R/plot.R
#'
#' @param seurat_obj Seurat object (Seurat ^3.0)
#' @param assay  seurat_obj assay to use
#' @param var_group  the group variable, character
#' @param slot  seurat_obj slot to use
#' @param vec_features  a vector of features to plot in the violin plot
#' @param vec_group_colors  a vector of colors, named by corresponding group. Length must match number of groups. Character
#' @param f_color  if vec_group_colors is not provided, the user may instead provide a function f_color() that takes as its only argument the number of colors to generate
#' @param flip  if TRUE (default), groups are rows and features are columns, and vice-verso for FALSE
#' @param do_plot  Whether to plot, logical
#' @param pt.size size of jitter in the violin plots. Set to 0 (default) to omit
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' top5genes<-rownames(Pagwas$gene_heritability_correlation)[order(Pagwas$gene_heritability_correlation,decreasing = T)[1:5]]
#' library("RColorBrewer")
#' library("ggplot2")

#' plot_vlnGrid(seurat_obj=scRNAexample,
#'              assay="RNA", slot="data",
#'              var_group="anno",# 细胞cluster注释列
#'              vec_features=top5genes,
#'              vec_group_colors= pal_d3(alpha =0.5)(10)
#' )

plot_vlnGrid = function(seurat_obj,
                        assay,
                        slot,
                        var_group,
                        vec_features,
                        vec_group_colors=NULL,
                        f_color = colorRampPalette(brewer.pal(n=11, name="RdYlBu")),
                        flip = T,
                        do_plot = F,
                        pt.size = 0,
                        feature_fontface = "bold.italic",
                        fontsize_axistext_x=12,
                        fontsize_axistext_y=12,
                        aspect.ratio =NULL
) {

  #=============prepare group and colors==================
  seurat_obj_tmp = seurat_obj
  Idents(seurat_obj_tmp) <- var_group
  levels(x = seurat_obj_tmp) = sort(unique(seurat_obj_tmp@meta.data[[var_group]]), decreasing = if (flip) T else F)

  if (is.null(vec_group_colors)) {
    n_group <- length(levels(x = seurat_obj_tmp))
    vec_group_colors <- f_color(n_group)
    names(vec_group_colors) <- levels(x = seurat_obj_tmp)
  }

  #=============generate plot list==================
  # produces a list of rows of violin plots, one per feature
  list_plot <- VlnPlot(object=seurat_obj_tmp,
                       assay=assay,
                       features = vec_features,
                       pt.size = pt.size,
                       cols=vec_group_colors,
                       sort = F,
                       #group.by = var_group,
                       same.y.lims = F,
                       slot=slot,
                       log = F,
                       combine = F,
                       flip=F)

  names(list_plot) <- vec_features

  if (is.null(aspect.ratio)) {
    aspect.ratio = 1.5*length(vec_group_colors)/length(vec_features)
    message(paste0("Using aspect ratio ", aspect.ratio))
  }

  list_plot_flip <- lapply(1:length(list_plot), function(i) {
    plot_tmp=list_plot[[i]]
    if (flip) {
      plot_tmp = plot_tmp +
        coord_flip()
    }
    plot_tmp = plot_tmp  +
      theme(
        plot.title = element_text(face = feature_fontface, size=fontsize_axistext_y),
        axis.text.y= element_blank(),
        plot.margin = margin(b=1, unit="cm"))
    #margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)

    if (i==1) {
      plot_tmp <- plot_tmp +
        theme(
          axis.text.y=element_text(
            hjust=1,
            vjust=0,
            size=fontsize_axistext_x,
            angle=30
          )
        )
    }


    plot_tmp <- plot_tmp +
      theme(
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(
          size= fontsize_axistext_x,
          angle=30),
        aspect.ratio=aspect.ratio,
        legend.position="none")
  })

  p <- patchwork::wrap_plots(... = list_plot_flip,
                             ncol = length(list_plot_flip)
  )

  if (do_plot) p

  return(p)
}


