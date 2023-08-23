
#' scPagwas_Visualization
#' @description Visualize the scPagwas score in Umap and Tsne.
#'
#' @param Single_data Single_data in seruat format ,the save with
#' scPagwas_main(), you'd better to run reduction of UMAP AND TSNE
#' @param output.dirs (character)default is NULL.the file folder name
#' for save the figures.NULL means no folder is created, no pdf figure
#' output.
#' @param FigureType (character)"tsne" or "umap
#' @param p_thre (numeric),default is 0.1, Threshold for pecent(<1) of
#' Positive cells for level of scPagwas_score.
#' @param width (numeric)Figure width
#' @param height (numeric)Figure height
#' @param size (numeric)size for scatters
#' @param title (character)Figure title
#' @param lowColor (character)Color for low scPagwas score
#' @param highColor (character)Color for high scPagwas score
#' @param do_plot  Whether to plot, logical
#'
#' @return figures for TRS score, gPas score,and p value reduction plot
#' @export
#'
#' @examples
#' load(system.file("extdata", "Pagwas_data.RData", package = "scPagwas"))
#' scPagwas_Visualization(
#'   Single_data = Pagwas_data,
#'   p_thre = 0.05,
#'   FigureType = "umap",
#'   width = 7,
#'   height = 7,
#'   lowColor = "white",
#'   highColor = "red",
#'   output.dirs = "figure",
#'   size = 0.5,
#'   do_plot = TRUE
#' )
#' @author Chunyu Deng
#' @aliases scPagwas_Visualization
#' @keywords scPagwas_Visualization, plot the reduction plot for
#' scPagwas result.

scPagwas_Visualization <- function(Single_data = NULL,
                                   p_thre = 0.05,
                                   output.dirs = "scPagwastest_output",
                                   FigureType = "tsne",
                                   width = 7,
                                   height = 7,
                                   title = "",
                                   lowColor = "#000957",
                                   highColor = "#EBE645",
                                   size = 0.5,
                                   do_plot = TRUE) {
  if (is.null(Single_data$scPagwas.gPAS.score)) {
    stop("ERROR: scPagwas.gPAS.score is NULL. scPagwas_score can be
         calculated by scPagwas_perform_score function!")
  }

  if (!dir.exists(output.dirs)) {
    dir.create(output.dirs)
  }

  UMAP_1 <- UMAP_2 <- TSNE_1 <- TSNE_2 <- NULL

  if (FigureType == "umap") {
    all_fortify_can <- fortify.Seurat.umap(Single_data)

    plot_scPagwas_score <- ggpubr::ggscatter(all_fortify_can,
      x = "UMAP_1", y = "UMAP_2",
      color = "scPagwas.gPAS.score", fill = "scPagwas.gPAS.score",
      size = size, title = title,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1) +
      ggtitle("scPagwas.gPAS.score")

    if (do_plot) print(plot_scPagwas_score)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.gPAS.score_umap.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score)
      grDevices::dev.off()
    }

    plot_scPagwas_score2 <- ggpubr::ggscatter(all_fortify_can,
      x = "UMAP_1", y = "UMAP_2",
      color = "scPagwas.TRS.Score1",
      fill = "scPagwas.TRS.Score1",
      size = size,
      title = title,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1) +
      ggtitle("scPagwas.TRS.Score")

    if (do_plot) print(plot_scPagwas_score2)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.TRS.Score_umap.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score2)
      grDevices::dev.off()
    }

    plot_scPagwas_score3 <- ggpubr::ggscatter(all_fortify_can,
                                              x = "UMAP_1", y = "UMAP_2",
                                              color = "scPagwas.downTRS.Score1",
                                              fill = "scPagwas.downTRS.Score1",
                                              size = size,
                                              title = title,
                                              repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = "#146C94") +
      scale_color_gradient(low = lowColor, high = "#146C94") +
      theme(aspect.ratio = 1) +
      ggtitle("scPagwas.downTRS.Score")

    if (do_plot) print(plot_scPagwas_score3)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.downTRS.Score_umap.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score3)
      grDevices::dev.off()
    }

    plot_scPagwas_score4 <- ggpubr::ggscatter(all_fortify_can,
                                              x = "UMAP_1", y = "UMAP_2",
                                              color = "scPagwas.upTRS.Score1",
                                              fill = "scPagwas.upTRS.Score1",
                                              size = size,
                                              title = title,
                                              repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1) +
      ggtitle("scPagwas.upTRS.Score")

    if (do_plot) print(plot_scPagwas_score4)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.upTRS.Score_umap.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score4)
      grDevices::dev.off()
    }
    plots_sigp1 <- ggplot() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Random_Correct_BG_adjp > p_thre, ],
        aes(x = UMAP_1, y = UMAP_2), size = size, alpha = 0.8,
        color = "#E4DCCF"
      ) +
      umap_theme() +
      # new_scale_color() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Random_Correct_BG_adjp <= p_thre, ],
        aes(x = UMAP_1, y = UMAP_2), color = "#EA5455", size = .2
      ) +
      umap_theme() +
      # new_scale_color() +
      ggtitle(paste0("Random_Correct_BG for TRS score, adj p<", p_thre, " significant cells"))


    if (do_plot) print(plots_sigp1)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas_Correct_BG_adjp",
          p_thre, "_umap.pdf"
        ),
        height = height, width = width-1
      )
      print(plots_sigp1)
      grDevices::dev.off()
    }

  }


  if (FigureType == "tsne") {
    all_fortify_can <- fortify.Seurat.tsne(Single_data)
    # globalVariables(names(all_fortify_can))
    plot_scPagwas_score <- ggpubr::ggscatter(all_fortify_can,
      x = "TSNE_1", y = "TSNE_2",
      color = "scPagwas.gPAS.score",
      fill = "scPagwas.gPAS.score",
      size = size,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      ggtitle("scPagwas.gPAS.score") + theme(aspect.ratio = 1)


    if (do_plot) print(plot_scPagwas_score)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.gPAS.score_tsne.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score)
      grDevices::dev.off()
    }
    plot_scPagwas_score2 <- ggpubr::ggscatter(all_fortify_can,
      x = "TSNE_1", y = "TSNE_2",
      color = "scPagwas.TRS.Score1",
      fill = "scPagwas.TRS.Score1",
      size = size,
      title = title,
      repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1) +
      ggtitle("scPagwas.TRS.Score")

    if (do_plot) print(plot_scPagwas_score2)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.TRS.Score_tsne.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score2)
      grDevices::dev.off()
    }

    plot_scPagwas_score3 <- ggpubr::ggscatter(all_fortify_can,
                                              x = "TSNE_1", y = "TSNE_2",
                                              color = "scPagwas.downTRS.Score1",
                                              fill = "scPagwas.downTRS.Score1",
                                              size = size,
                                              title = title,
                                              repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = "#146C94") +
      scale_color_gradient(low = lowColor, high = "#146C94") +
      theme(aspect.ratio = 1) +
      ggtitle("scPagwas.downTRS.Score")

    if (do_plot) print(plot_scPagwas_score3)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.downTRS.Score_tsne.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score3)
      grDevices::dev.off()
    }

    plot_scPagwas_score4 <- ggpubr::ggscatter(all_fortify_can,
                                              x = "TSNE_1", y = "TSNE_2",
                                              color = "scPagwas.upTRS.Score1",
                                              fill = "scPagwas.upTRS.Score1",
                                              size = size,
                                              title = title,
                                              repel = TRUE
    ) + umap_theme() +
      scale_fill_gradient(low = lowColor, high = highColor) +
      scale_color_gradient(low = lowColor, high = highColor) +
      theme(aspect.ratio = 1) +
      ggtitle("scPagwas.upTRS.Score")

    if (do_plot) print(plot_scPagwas_score4)
    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas.upTRS.Score_tsne.pdf"
        ),
        height = height, width = width
      )
      print(plot_scPagwas_score3)
      grDevices::dev.off()
    }
    plots_sigp1 <- ggplot() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Random_Correct_BG_adjp > p_thre, ],
        aes(x = TSNE_1, y = TSNE_2),
        size = size, alpha = 0.8,
        color = "#E4DCCF"
      ) +
      umap_theme() +
      # new_scale_color() +
      geom_point(
        data = all_fortify_can[all_fortify_can$Random_Correct_BG_adjp <= p_thre, ],
        aes(x = TSNE_1, y = TSNE_2), color = "#EA5455", size = size
      ) +
      umap_theme() +
      # new_scale_color() +
      ggtitle(paste0("Random_Correct_BG_adjp<", p_thre, " significant cells"))

    if (do_plot) print(plots_sigp1)

    if (!is.null(output.dirs)) {
      grDevices::pdf(
        file = paste0(
          "./", output.dirs,
          "/scPagwas_Random_Correct_BG_adjp",
          p_thre, "_tsne.pdf"
        ),
        height = height, width = width-1
      )
      print(plots_sigp1)
      grDevices::dev.off()
    }
  }
  # return(all_fortify_can)
}

#' umap_theme
#' @description ggplot2 theme for umap plot
#' @export
#'

umap_theme <- function() {
  theme_grey() %+replace%
    theme(
      panel.background = element_rect(
        fill = "white",
        colour = "black",
        size = 1
      ),
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
#'
#' @param x Seurat format
#' @export

fortify.Seurat.umap <- function(x) {
  xy1 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "umap")
  )
  colnames(xy1) <- c("UMAP_1", "UMAP_2")
  xy1$UMAP_1 <- as.numeric(xy1$UMAP_1)
  xy1$UMAP_2 <- as.numeric(xy1$UMAP_2)

  return(cbind(xy1, as.data.frame(x@meta.data)))
}


#' fortify.Seurat.tsne
#' @description set data frame to ggplot
#'
#' @param x seruat
#' @export

fortify.Seurat.tsne <- function(x) {
  xy2 <- as.data.frame(
    Seurat::Embeddings(x, reduction = "tsne")
  )
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)

  return(cbind(xy2, as.data.frame(x@meta.data)))
}
