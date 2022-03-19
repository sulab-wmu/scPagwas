
#' scPagwas_Visualization
#' @description Visualize the scPagwas score in Umap and Tsne.
#'
#' @param scPagwas_score scPagwas_score for scPagwas_main result
#' @param Single_data Single_data in seruat format ,the save with scPagwas_main(), you'd better to run reduction of UMAP AND TSNE
#' @param Reduction (logical) default is FALSE. Whether to run the Reduction for Single_data.If you are do it before,ignore it.
#' @param assay (character)"RNA" or "SCT", It depens on the Single_data.
#' @param filename (character)default is NULL.the file folder name for save the figures.NULL means no folder is created, print the figure to current environment, check setwd().
#' @param FigureType (character)"tsne" or "umap
#' @param cellpercent (numeric),default is 0.1, Threshold for pecent(<1) of Positive cells for level of scPagwas_score.
#' @param width (numeric)Figure width
#' @param height (numeric)Figure height
#' @param size (numeric)size for scatters
#' @param npcs (integr) the parameter of Seurat::RunPCA.
#' @param title (character)Figure title
#' @param lowColor (character)Color for low scPagwas score
#' @param highColor (character)Color for high scPagwas score
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
#' #  Reduction = FALSE,
#' #  assay = "SCT",
#' #  cellpercent = 0.1,
#' #  filename = "scpagwas_pbc_False",
#' #  FigureType = "tsne",
#' #  width = 7,
#' #  height = 7,
#' #  lowColor = "#000957", highColor = "#EBE645",
#' #  size = 0.5,
#' #  title = "scPagwas_score"
#' #)
scPagwas_Visualization <- function(scPagwas_score = NULL,
                                   Single_data=NULL,
                                   Reduction = FALSE,
                                   assay = "SCT",
                                   cellpercent = 0.1,
                                   filename = NULL,
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   lowColor = "#000957", highColor = "#EBE645",
                                   size = 0.5,
                                   npcs=50,
                                   title = "scPagwas_score") {
  #suppressMessages(require(ggnewscale))
  #suppressMessages(require(ggtext))
  #suppressMessages(require(ggrepel))

  #if(colnames(Single_data@meta.data) %in% "tsne")

  if (is.null(scPagwas_score)) {
    stop("ERROR: scPagwas_score is NULL. scPagwas_score can be calsulated by scPagwas_perform_score function!")
  }

  if (sum(colnames(Single_data) %in% names(scPagwas_score)) != length(names(scPagwas_score))) {
    message("There is blank , '-' or '+' within cell names!")
    colnames(Single_data) <- stringr::str_replace_all(colnames(Single_data), " ", ".")
    colnames(Single_data) <- stringr::str_replace_all(colnames(Single_data), "\\+", ".")
    colnames(Single_data) <- stringr::str_replace_all(colnames(Single_data), "-", ".")
  }

  Single_data <- Single_data[, intersect(colnames(Single_data),names(scPagwas_score))]
  scPagwas_score <- scPagwas_score[intersect(colnames(Single_data),names(scPagwas_score))]
  Single_data$scPagwas_score <- scPagwas_score

  if (Reduction) {
    Single_data <- suppressMessages(Seurat::RunPCA(object = Single_data, assay = assay, npcs = npcs))
    if (FigureType == "tsne") {
      Single_data <- suppressMessages(Seurat::RunTSNE(object = Single_data, assay = assay, reduction = "pca", dims = 1:npcs))
    } else if (FigureType == "umap") {
      Single_data <- suppressMessages(Seurat::RunUMAP(object = Single_data, assay = assay, reduction = "pca", dims = 1:npcs))
    } else {
      stop("ERROR:: The Reduction is TRUE, but no correct FigureType is choose, either tsne or umap")
    }
  }

  if(!is.null(filename)){
  if (!file.exists(paste0("./", filename))) {
    dir.create(filename)
  }
  }
  # num<-length(unique(as.vector(Idents(Single_data))))
  Single_data <- Single_data[, !is.na(Single_data$scPagwas_score)]

  thre <- sort(Single_data$scPagwas_score, decreasing = T)[ncol(Single_data) * cellpercent]

  if (FigureType == "umap") {
    all_fortify_can <- fortify.Seurat.umap(Single_data)
    plot_scPagwas_score <- ggpubr::ggscatter(all_fortify_can,
      x = "UMAP_1", y = "UMAP_2",
      color = "scPagwas_score", fill = "scPagwas_score", size = size, title = title,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1)
    print(plot_scPagwas_score)
    pdf(file = paste0("./", filename, "/scPagwas_score_umap.pdf"), height = height, width = width)
    print(plot_scPagwas_score)
    dev.off()

    plots_sigp1 <- ggplot() +
      geom_point(
        data = all_fortify_can[all_fortify_can$scPagwas_score <= thre, ],
        aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8, color = "gray"
      ) +
      umap_theme() +
      #new_scale_color() +
      geom_point(
        data = all_fortify_can[all_fortify_can$scPagwas_score > thre, ],
        aes(x = UMAP_1, y = UMAP_2), color = "#F90716", size = .2
      ) +
      umap_theme() +
      #new_scale_color() +
      ggtitle(paste0("Top ", cellpercent * 100, "% cells"))

    print(plots_sigp1)
    pdf(file = paste0("./", filename, "/scPagwas_TOP", cellpercent, "_umap.pdf"), height = height, width = width)
    print(plots_sigp1)
    dev.off()
  }


  if (FigureType == "tsne") {
    all_fortify_can <- fortify.Seurat.tsne(Single_data)

    plot_scPagwas_score <- ggpubr::ggscatter(all_fortify_can,
      x = "TSNE_1", y = "TSNE_2",
      color = "scPagwas_score", fill = "scPagwas_score", size =size,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      ggtitle("scPagwas_score") + theme(aspect.ratio = 1)

    print(plot_scPagwas_score)
    pdf(file = paste0("./", filename, "/scPagwas_score_tsne.pdf"), height = height, width = width)
    print(plot_scPagwas_score)
    dev.off()

    plots_sigp1 <- ggplot() +
      geom_point(
        data = all_fortify_can[all_fortify_can$scPagwas_score <= thre, ],
        aes(x = TSNE_1, y = TSNE_2), size = size, alpha = 0.8, color = "gray"
      ) +
      umap_theme() +
      #new_scale_color() +
      geom_point(
        data = all_fortify_can[all_fortify_can$scPagwas_score > thre, ],
        aes(x = TSNE_1, y = TSNE_2), color = "#F90716", size = size
      ) +
      umap_theme() +
      #new_scale_color() +
      ggtitle(paste0("Top ", cellpercent * 100, "% cells"))

    print(plots_sigp1)
    pdf(file = paste0("./", filename, "/scPagwas_TOP", cellpercent, "_tsne.pdf"), height = height, width = width)
    print(plots_sigp1)
    dev.off()
  }
  # return(all_fortify_can)
}


#' umap_theme
#' @description ggplot2 theme for umap plot
#' @return

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
#'
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
#'
#' @return

fortify.Seurat.tsne <- function(x) {
  xy2 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "tsne"))
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)

  return(cbind(xy2, as.data.frame(x@meta.data)))
}
