
#' plot_pathway_contribution_network
#'
#' @description plot the network based on pathway contribution
#' oaqc
#'
#' Thanks to the SCOPfunctions package
#' https://github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions/blob/main/R/plot.R
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
#' @param do_plot Whether to plot, logical
#' @param figurenames The filename and address of the output plot
#' @param width figure width
#' @param height figure height
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
                                        Celltype_anno,
                                        pathways_thre=0.05,
                                        adjacency.power=4,
                                        adjacency.thre=0.1,
                                        n_max_pathways=10,
                                        igraph_algorithm = "drl",
                                        fontface_labels="bold.italic",
                                        color_point= pal_d3(alpha =0.5)(10),
                                        color_edge = "#9D9D9D",
                                        fontSize_label_lg=3,
                                        fontSize_legend_lg=1,
                                        fontSize_legend_xlg=1,
                                        limits=c(0,10),
                                        edge_thickness = 2,
                                        do_plot=F ,
                                        figurenames=NULL,
                                        width=7,
                                        height=7
                                        ) {

  # sort
  #mat_datExpr
  cl<-unique((Celltype_anno$annotation))

  lapply(cl, function(ss){
    tt<-Celltype_anno$annotation==ss
    mat_datExpr_cell<-mat_datExpr[tt,]
    pathwaycontribution<-vec_pathwaycontribution[,ss]
    #pathwayhighlight = pathwaycontribution[pathwaycontribution<pathways_thre]
    #pathwayhighlight = sort(pathwayhighlight, decreasing = F)
    #if(length(pathwayhighlight)>n_max_pathways){
     # pathwayhighlight<-pathwayhighlight[1:10]
    pathwaycontribution<- -log10(pathwaycontribution)
    #}
    #pathwaycontribution = sort(pathwaycontribution, decreasing = F)
  # take take hub pathways and pathways to highlight (if any)
  #vec_logical = names(vec_pathwaycontribution) %in% c(names(vec_pathwaycontribution)[1:min(n_max_pathways, length(vec_pathwaycontribution))],  vec_pathways_highlight)
  #vec_pathwaycontribution = vec_pathwaycontribution[vec_logical]

    mat_datExpr_cell <- mat_datExpr_cell[,names(pathwaycontribution)]

  # Compute network adjacency
  mat_adj <- WGCNA::adjacency(mat_datExpr_cell,
                              power =adjacency.power,
                              corFnc = "cor",
                              corOptions = list(use = "p"),
                              type = "signed hybrid")

  df_hub.data <- cbind.data.frame(pathway = names(pathwaycontribution), contribution = pathwaycontribution)
  df_hub.data$color<-rep("none",nrow(df_hub.data))
  #df_hub.data$color<-color_point
  #df_hub.data$color[which(df_hub.data$pathway %in% vec_pathways_highlight)]<-"Highlight"
  #df_hub.data$contribution[df_hub.data$contribution>limits[2]]<-10


  mat_adj[which(mat_adj< adjacency.thre)]<-0
  mat_adj<-mat_adj[rowSums(mat_adj,na.rm = TRUE)!=0,colSums(mat_adj,na.rm = TRUE)!=0]


  df_hub.data<-df_hub.data[rownames(mat_adj),]
  a<-df_hub.data$contribution*edge_thickness
  a[a>limits[2]]<-limits[2]

  mat_adj %>%
    igraph::graph.adjacency(mode = "undirected", weighted = T, diag = FALSE) %>%
    tidygraph::as_tbl_graph() %>% igraph::upgrade_graph() %>% tidygraph::activate(nodes) %>%
    dplyr::mutate(contribution = a,color=df_hub.data$color) %>%
    tidygraph::activate(edges) %>%
    tidygraph::activate(nodes) %>%
    dplyr::filter(!tidygraph::node_is_isolated())-> hub.plot
 # library(graphlayouts)
 # library("oaqc")
 # library(ggforce)
 # library("concaveman")

  #set.seed(665)
  #create network with a group structure
  #hub.plot <- sample_islands(9, 40, 0.4, 15)
  hub.plot <- igraph::simplify(hub.plot)
  V(hub.plot)$grp <- as.character(rep(1,length(V(hub.plot))))
  bb <- layout_as_backbone(hub.plot, keep = 0.4)
  E(hub.plot)$col <- F
  E(hub.plot)$col[bb$backbone] <- T

  #bb <- graphlayouts::layout_as_backbone(hub.plot, keep = 0.4)
  #E(hub.plot)$col <- FALSE
  #E(hub.plot)$col[bb$backbone] <- TRUE

  #layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]
#cent = graph.strength(hub.plot)
 # ggraph::ggraph(hub.plot, layout = "stress") +
  cell_p<-  ggraph::ggraph(hub.plot, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2])+
    ggraph::geom_edge_link0(aes(edge_width = weight), edge_colour = "grey66") +
    ggraph::geom_node_point(aes(fill = grp, size = contribution), shape = 21) +
    ggraph::geom_node_text(aes(filter = contribution >=7,
                               size=fontSize_label_lg,
                               label = name)) +
    #scale_fill_manual(values = df_hub.data$color) +
    scale_edge_width(range = c(0, 0.7)) +
    scale_size(#breaks = c(1,2,3,4,5,6,7,8,9,10),
               limits = limits,
               range = limits
               #labels = c(" \U00B1 1", " \U00B1 2", " \U00B1 3")
    )+
    geom_mark_hull(
        aes(x, y, group = grp, fill = grp),
        concavity = 4,
        expand = unit(2, "mm"),
        alpha = 0.25
      )+
    theme_graph() +
      scale_color_brewer(values="") +
      scale_fill_brewer(values="")+
    #theme(legend.position = "none")
  theme(legend.position = "bottom")

  return(cell_p)

  })
    #coord_cartesian(clip="off")

  if (do_plot) print(module.plot)

    if(!is.null(figurenames)){
      pdf(file = figurenames, width = width, height = height)
      print(module.plot)
      dev.off()
    }

}
