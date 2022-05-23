
#' scPagwas_Visualization
#' @description Visualize the scPagwas score in Umap and Tsne.
#'
#' @param scPagwas_score scPagwas_score for scPagwas_main result
#' @param scPagwas_p
#' @param Single_data Single_data in seruat format ,the save with scPagwas_main(), you'd better to run reduction of UMAP AND TSNE
#' @param output.dirs (character)default is NULL.the file folder name for save the figures.NULL means no folder is created, no pdf figure output.
#' @param FigureType (character)"tsne" or "umap
#' @param p_thre (numeric),default is 0.1, Threshold for pecent(<1) of Positive cells for level of scPagwas_score.
#' @param width (numeric)Figure width
#' @param height (numeric)Figure height
#' @param size (numeric)size for scatters
#' @param npcs (integr) the parameter of Seurat::RunPCA.
#' @param title (character)Figure title
#' @param lowColor (character)Color for low scPagwas score
#' @param highColor (character)Color for high scPagwas score
#' @param do_plot  Whether to plot, logical
#'
#' @return
#' @export
#'
#' @examples
#' library(scPagwas)
#' data(scRNAexample)
#' #scPagwas_Visualization(
#' #  scPagwas = Pagwas,
#' #  Single_data = scRNAexample,
#' #  cellpercent = 0.1,
#' #  output.dirs = "scpagwas_pbc_False",
#' #  FigureType = "tsne",
#' #  width = 7,
#' #  height = 7,
#' #  lowColor = "#000957", highColor = "#EBE645",
#' #  size = 0.5,
#' #  title = "scPagwas_score"
#' #)
scPagwas_Visualization <- function(Single_data=NULL,
                                   p_thre = 0.05,
                                   output.dirs="scPagwastest_output",
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", highColor = "#EBE645",
                                   size = 0.5,
                                   do_plot = F) {


  if (is.null(Single_data$scPagwas.lm.score)) {
    stop("ERROR: scPagwas.lm.score is NULL. scPagwas_score can be calculated by scPagwas_perform_score function!")
  }
  if (is.null(Single_data$Cells.lm.rankPvalue)) {
    stop("ERROR: Cells.lm.rankPvalue is NULL.")
  }
  if(!dir.exists(output.dirs)){
    dir.create(output.dirs)
  }

  if (FigureType == "umap") {
    all_fortify_can <- fortify.Seurat.umap(Single_data)
    plot_scPagwas_score <- ggpubr::ggscatter(all_fortify_can,
                                             x = "UMAP_1", y = "UMAP_2",
      color = "scPagwas.lm.score", fill = "scPagwas.lm.score", size = size, title = title,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1)+
      ggtitle("scPagwas.lm.score")

    if(do_plot) print(plot_scPagwas_score)
    if(!is.null(output.dirs)){
    pdf(file = paste0("./", output.dirs, "/scPagwas.lm.score_umap.pdf"), height = height, width = width)
    print(plot_scPagwas_score)
    dev.off()
    }

    plot_scPagwas_score2 <- ggpubr::ggscatter(all_fortify_can,
                                              x = "UMAP_1", y = "UMAP_2",
                                             color = "scPagwas.topgenes.Score1",
                                             fill = "scPagwas.topgenes.Score1",
                                             size = size,
                                             title = title,
                                             repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1)+
      ggtitle("scPagwas.topgenes.Score")

    if(do_plot) print(plot_scPagwas_score2)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.topgenes.Score_umap.pdf"), height = height, width = width)
      print(plot_scPagwas_score2)
      dev.off()
    }

    plots_sigp1 <- ggplot() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Cells.lm.rankPvalue > p_thre, ],
        aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8, color = "gray"
      ) +
      umap_theme() +
      #new_scale_color() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Cells.lm.rankPvalue <= p_thre, ],
        aes(x = UMAP_1, y = UMAP_2), color = "#F90716", size = .2
      ) +
      umap_theme() +
      #new_scale_color() +
      ggtitle(paste0("p<", p_thre, " significant cells"))


    if(do_plot) print(plots_sigp1)
    if(!is.null(output.dirs)){
    pdf(file = paste0("./", output.dirs, "/scPagwas_p", p_thre, "_umap.pdf"), height = height, width = width)
    print(plots_sigp1)
    dev.off()
    }
  }


  if (FigureType == "tsne") {
    all_fortify_can <- fortify.Seurat.tsne(Single_data)

    plot_scPagwas_score <- ggpubr::ggscatter(all_fortify_can,
      x = "TSNE_1", y = "TSNE_2",
      color = "scPagwas.lm.score", fill = "scPagwas.lm.score", size =size,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      ggtitle("scPagwas.lm.score") + theme(aspect.ratio = 1)


    if(do_plot) print(plot_scPagwas_score)
    if(!is.null(output.dirs)){
    pdf(file = paste0("./", output.dirs, "/scPagwas.lm.score_tsne.pdf"), height = height, width = width)
    print(plot_scPagwas_score)
    dev.off()
    }
    plot_scPagwas_score2 <- ggpubr::ggscatter(all_fortify_can,
                                              x = "TSNE_1", y = "TSNE_2",
                                              color = "scPagwas.topgenes.Score1",
                                              fill = "scPagwas.topgenes.Score1",
                                              size = size,
                                              title = title,
                                              repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1)+
      ggtitle("scPagwas.topgenes.Score")

    if(do_plot) print(plot_scPagwas_score2)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.topgenes.Score_tsne.pdf"), height = height, width = width)
      print(plot_scPagwas_score2)
      dev.off()
    }


    plots_sigp1 <- ggplot() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Cells.lm.rankPvalue > p_thre, ],
        aes(x = TSNE_1, y = TSNE_2), size = size, alpha = 0.8, color = "gray"
      ) +
      umap_theme() +
      #new_scale_color() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Cells.lm.rankPvalue <= p_thre, ],
        aes(x = TSNE_1, y = TSNE_2), color = "#F90716", size = size
      ) +
      umap_theme() +
      #new_scale_color() +
      ggtitle(paste0("p<", p_thre, " significant cells"))

    if (do_plot) print(plots_sigp1)

    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas_p", p_thre, "_tsne.pdf"), height = height, width = width)
      print(plots_sigp1)
    dev.off()
    }
  }
  # return(all_fortify_can)
}


#' preprogress_PaSingle_data
#'
#' @param Single_data
#' @param npcs
#' @param nfeatures
#'
#' @return
#' @export
#'
#' @examples
preprogress_PaSingle_data<-function(Single_data,
                       npcs=10,
                       nfeatures=20){

  SeuratObject::DefaultAssay(Single_data) = 'scPagwasPaHeritability'
  b<-SeuratObject::GetAssayData(Single_data, slot="data",assay ="scPagwasPaHeritability")
  # if(b)
  b<-t(apply(b, 1, function(x) (x-min(x))/(max(x)-min(x))))

  Single_data<-SeuratObject::SetAssayData(Single_data,
                                          slot="data",
                                          assay ="scPagwasPaHeritability",
                                          new.data = b)
  Single_data <- suppressWarnings(Seurat::FindVariableFeatures(Single_data,
                                              assay = 'scPagwasPaHeritability',
                                              nfeatures =nfeatures,
                                              selection.method = "mvp"
  ))
  Single_data <- ScaleData(Single_data)
  Single_data <-suppressWarnings(Seurat::RunPCA(object = Single_data,
                                assay = "scPagwasPaHeritability",
                                npcs = npcs))

  Single_data <- suppressWarnings(Seurat::RunTSNE(object = Single_data,
                                 assay = "scPagwasPaHeritability",
                                 reduction = "pca",
                                 dims = 1: npcs))
  Single_data <- suppressWarnings(Seurat::RunUMAP(object = Single_data,
                                                  assay ="scPagwasPaHeritability",
                                                  reduction = "pca",
                                                  dims = 1: npcs))
  return(Single_data)
}
#' scPagwasPaHeritability_Visualization
#'
#' @param Single_data
#' @param p_thre
#' @param output.dirs
#' @param FigureType
#' @param width
#' @param height
#' @param npcs
#' @param nfeatures
#' @param lowColor
#' @param highColor
#' @param size
#' @param title
#' @param do_plot
#'
#' @return
#' @export
#'
#' @examples
scPagwasPaHeritability_Visualization <- function(Single_data=NULL,
                                    p_thre = 0.05,
                                    output.dirs="test",
                                    FigureType = "tsne",
                                    width = 7,
                                    height = 7,
                                    lowColor="#8479E1",
                                    midColor="#F7F5F2",
                                    highColor="#FD5D5D",
                                    p_color= "#E45826",
                                    size = 0.5,
                                    do_plot = F) {
  ##debug
  # Single_data=Pagwas;
  # p_thre = 0.05;
  # output.dirs="adPrunetest";
  # FigureType = "tsne";
  # width = 7;
  # height = 7;
  # lowColor="#8479E1";
  # midColor="#F7F5F2";
  # highColor="#FD5D5D";
  # p_color= "#E45826";
  # size = 0.5;
  # do_plot = T;


  if(!dir.exists(output.dirs)){
    dir.create(output.dirs)
  }

  if (FigureType == "umap") {
    Pagwas_fortify <- fortify.Seurat.umap(Single_data)
    plot1 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = UMAP_1, y = UMAP_2,color =scPagwas.lm.score), size = size, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor)+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
    ggtitle("scPagwas.lm.score in scPagwasPaHeritability")
    if(do_plot) print(plot1)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.lm.score_scPagwasPaHeritability.umap.pdf"), height = height, width = width)
      print(plot1)
      dev.off()
    }

    plot2 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = UMAP_1, y = UMAP_2,color =scPagwas.topgenes.Score1), size = size, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor,
                             midpoint = median(Pagwas_fortify$scPagwas.topgenes.Score1))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("scPagwas.topgenes.Score in scPagwasPaHeritability")
    if(do_plot) print(plot2)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.topgenes.Score_scPagwasPaHeritability.umap.pdf"), height = height, width = width)
      print(plot2)
      dev.off()
    }

    Pagwas_fortify$p_thre<-rep("non",nrow(Pagwas_fortify))
    Pagwas_fortify$p_thre[which(Pagwas_fortify$Cells.lm.rankPvalue<p_thre)]<-"significant"
    plot3 <-  ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = UMAP_1, y = UMAP_2,color =p_thre), size =size, alpha = 1) +
      umap_theme() +
      scale_color_manual(values = c("non" = "#05595B", "significant" = p_color))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      #new_scale_color() +
      ggtitle(paste0("p<", p_thre, " significant cells"))


    if(do_plot) print(plot3)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas_p", p_thre, "_scPagwasPaHeritability.umap.pdf"), height = height, width = width)
      print(plot3)
      dev.off()
    }
  }


  if (FigureType == "tsne") {
    Pagwas_fortify <- fortify.Seurat.tsne(Single_data)

    plot1 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,
                     color =scPagwas.lm.score),
                 size = size) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor,)+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("scPagwas.lm.score in scPagwasPaHeritability")

    if(do_plot) print(plot1)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.lm.score_scPagwasPaHeritability.tsne.pdf"), height = height, width = width)
      print(plot1)
      dev.off()
    }

    plot2 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,color =scPagwas.topgenes.Score1), size = size, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor,
                             midpoint = median(Pagwas_fortify$scPagwas.topgenes.Score1))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("scPagwas.topgenes.Score in scPagwasPaHeritability")

    if(do_plot) print(plot2)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.topgenes.Score_scPagwasPaHeritability.tsne.pdf"), height = height, width = width)
      print(plot2)
      dev.off()
    }

    Pagwas_fortify$p_thre<-rep("non",nrow(Pagwas_fortify))
    Pagwas_fortify$p_thre[which(Pagwas_fortify$Cells.lm.rankPvalue<p_thre)]<-"significant"
    plot3 <-  ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,color =p_thre), size =size, alpha = 1) +
      umap_theme() +
      scale_color_manual(values = c("non" = "#05595B", "significant" = p_color))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle(paste0("p<", p_thre, " significant cells"))


    if(do_plot) print(plot3)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas_p", p_thre, "_scPagwasPaHeritability.tsne.pdf"), height = height, width = width)
      print(plot3)
      dev.off()
    }
  }
  # return(all_fortify_can)
}




scPagwasPaHeritability_Pagene <- function(Single_data=NULL,
                                          gene=NULL,
                                          pathway=NULL,
                                          output.dirs="test",
                                          FigureType = "tsne",
                                          width = 7,
                                          height = 7,
                                          lowColor="#8479E1",
                                          midColor="#F7F5F2",
                                          highColor="#FD5D5D",
                                          p_color= "#E45826",
                                          size = 0.5,
                                          do_plot = F) {

  if(!dir.exists(output.dirs)){
    dir.create(output.dirs)
  }
  data_mat<- Single_data
  if (FigureType == "umap") {
    Pagwas_fortify <- fortify.Seurat.umap(Single_data)

    Pagwas_fortify$gene<-as.vector(data_mat[gene,])
    Pagwas_fortify$Pa<-as.vector(pca_scCell_mat[Pa,])

    plot3 <-  ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = UMAP_1, y = UMAP_2,color =gene), size =0.6, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low="#8479E1",mid="#F7F5F2",high="#FD5D5D")+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      labs(title = gene)
    #pdf("E:/OneDrive/GWAS_Multiomics/ad_test/5.6test/ad_pruneGSE138852/umap_opc_scPagwas_p.pdf")
    print(plot3)


    plot1 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = UMAP_1, y = UMAP_2,color =scPagwas.lm.score), size = size, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor)+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("scPagwas.lm.score in scPagwasPaHeritability")
    if(do_plot) print(plot1)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.lm.score_scPagwasPaHeritability.umap.pdf"), height = height, width = width)
      print(plot1)
      dev.off()
    }

    plot2 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = UMAP_1, y = UMAP_2,color =scPagwas.topgenes.Score1), size = size, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor,
                             midpoint = median(Pagwas_fortify$scPagwas.topgenes.Score1))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("scPagwas.topgenes.Score in scPagwasPaHeritability")
    if(do_plot) print(plot2)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.topgenes.Score_scPagwasPaHeritability.umap.pdf"), height = height, width = width)
      print(plot2)
      dev.off()
    }

    Pagwas_fortify$p_thre<-rep("non",nrow(Pagwas_fortify))
    Pagwas_fortify$p_thre[which(Pagwas_fortify$Cells.lm.rankPvalue<p_thre)]<-"significant"
    plot3 <-  ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = UMAP_1, y = UMAP_2,color =p_thre), size =size, alpha = 1) +
      umap_theme() +
      scale_color_manual(values = c("non" = "#05595B", "significant" = p_color))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      #new_scale_color() +
      ggtitle(paste0("p<", p_thre, " significant cells"))


    if(do_plot) print(plot3)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas_p", p_thre, "_scPagwasPaHeritability.umap.pdf"), height = height, width = width)
      print(plot3)
      dev.off()
    }
  }


  if (FigureType == "tsne") {
    Pagwas_fortify <- fortify.Seurat.tsne(Single_data)

    plot1 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,
                     color =scPagwas.lm.score),
                 size = size) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor,)+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("scPagwas.lm.score in scPagwasPaHeritability")

    if(do_plot) print(plot1)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.lm.score_scPagwasPaHeritability.tsne.pdf"), height = height, width = width)
      print(plot1)
      dev.off()
    }

    plot2 <- ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,color =scPagwas.topgenes.Score1), size = size, alpha = 1) +
      umap_theme() +
      scale_colour_gradient2(low=lowColor,mid=midColor,high=highColor,
                             midpoint = median(Pagwas_fortify$scPagwas.topgenes.Score1))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle("scPagwas.topgenes.Score in scPagwasPaHeritability")

    if(do_plot) print(plot2)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas.topgenes.Score_scPagwasPaHeritability.tsne.pdf"), height = height, width = width)
      print(plot2)
      dev.off()
    }

    Pagwas_fortify$p_thre<-rep("non",nrow(Pagwas_fortify))
    Pagwas_fortify$p_thre[which(Pagwas_fortify$Cells.lm.rankPvalue<p_thre)]<-"significant"
    plot3 <-  ggplot() +
      geom_point(data = Pagwas_fortify,
                 aes(x = TSNE_1, y = TSNE_2,color =p_thre), size =size, alpha = 1) +
      umap_theme() +
      scale_color_manual(values = c("non" = "#05595B", "significant" = p_color))+
      theme(aspect.ratio=1) +
      guides(colour = guide_legend(override.aes = list(size=3)))+
      ggtitle(paste0("p<", p_thre, " significant cells"))


    if(do_plot) print(plot3)
    if(!is.null(output.dirs)){
      pdf(file = paste0("./", output.dirs, "/scPagwas_p", p_thre, "_scPagwasPaHeritability.tsne.pdf"), height = height, width = width)
      print(plot3)
      dev.off()
    }
  }
  # return(all_fortify_can)
}




#' umap_theme
#' @description ggplot2 theme for umap plot
#' @return
#' @export

umap_theme <- function() {
  theme_grey() %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 2),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key = element_blank()
    )
}




#' fortify.Seurat.umap
#' @description set data frame to ggplot
#' @param x Seurat format
#' @export
#' @return

fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap"))
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)

  return(cbind(xy1, as.data.frame(x@meta.data)))
}


#' fortify.Seurat.tsne
#' @description set data frame to ggplot
#' @param x seruat
#' @export
#' @return

fortify.Seurat.tsne <- function(x) {
  xy2 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "tsne"))
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)

  return(cbind(xy2, as.data.frame(x@meta.data)))
}
