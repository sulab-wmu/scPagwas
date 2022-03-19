
#' Thanks to the SCOPfunctions package
#' https://github.com/CBMR-Single-Cell-Omics-Platform/SCOPfunctions/blob/main/R/plot.R
NULL

#' plot_bar_positie_nagtive
#' @description Generate barplot of identity group composition
#'
#' Generate a percentage barplot that shows the composition of each identity
#' (e.g. sample)in terms of groups (e.g. positive and negative cells for scPagwas)
#'
#' @param seurat_obj Seurat object (Seurat ^3.0)
#' @param var_ident the identify variable, character
#' @param var_group the group variable, character
#' @param vec_group_colors a vector of colors, named by corresponding group. Length must match number of groups. Character
#' @param f_color if vec_group_colors is not provided, the user may instead provide a function f_color() that takes as its only argument the number of colors to generate
#' @param do_plot Whether to plot, logical
#' @param title NULL to leave out
#' @param fontsize_title NULL to leave out
#' @param fontsize_axistitle_x NULL to leave out
#' @param fontsize_axistitle_y NULL to leave out
#' @param fontsize_axistext_x NULL to leave out
#' @param fontsize_axistext_y NULL to leave out
#' @param fontsize_legendtitle NULL to leave out
#' @param fontsize_legendtext NULL to leave out
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' p <- plot_barIdentGroup(seurat_obj=seu, var_ident="sample",var_group="cluster")
plot_bar_positie_nagtive <- function(seurat_obj,
                              var_ident,
                              var_group,
                              vec_group_colors=NULL,
                              f_color=colorRampPalette(brewer.pal(n=11, name="RdYlBu")),
                              do_plot = F,
                              title = NULL,
                              fontsize_title = 24,
                              fontsize_axistitle_x = 18,
                              fontsize_axistitle_y = 18,
                              fontsize_axistext_x = 12,
                              fontsize_axistext_y = 12,
                              fontsize_legendtitle = 12,
                              fontsize_legendtext = 10,
                              aspect.ratio=1.2) {

  #===============data.table with sums==================
  dt = data.table("ident" = as.character(seurat_obj@meta.data[[var_ident]]),
                  "group" = as.character(seurat_obj@meta.data[[var_group]]))
  dt[,n_ident := paste0(ident," (n=",.N, ")"), by=ident]
  vec_factorLevels <- dt$n_ident[gsub("\\ .*","",dt$n_ident) %>% as.numeric %>% order] %>% unique
  dt[,n_ident := factor(n_ident, levels = vec_factorLevels, ordered=T),]
  dt_sum <- dt[,.N, by=.(n_ident,group)]

  #===============ggplot==================
  # colors
  if (is.null(vec_group_colors)) {
    n_group <- length(unique(dt$group))
    vec_group_colors <- f_color(n_group)
    names(vec_group_colors) <- unique(dt$group)
  }

  p <- ggplot(dt_sum,
              aes(x = n_ident, y=N, fill = factor(group))) +

    geom_bar(
      position="fill",
      stat="identity",
      width=0.6,
      show.legend = if (!is.null(fontsize_legendtext)) TRUE else FALSE
      #position=position_dodge()
    ) +

    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values=vec_group_colors) +

    theme(
      axis.title.x = if (is.null(fontsize_axistitle_x)) element_blank() else element_text(size=fontsize_axistitle_x, vjust=0),
      axis.text.x = if (is.null(fontsize_axistext_x)) element_blank() else element_text(angle = 90, size=fontsize_axistext_x,vjust=0.5),
      axis.title.y = if (is.null(fontsize_axistitle_y)) element_blank() else element_text(size=fontsize_axistitle_y),
      axis.text.y = if (is.null(fontsize_axistext_y)) element_blank() else  element_text(size=fontsize_axistext_y),
      legend.title = if (is.null(fontsize_legendtext)) element_blank() else element_text(size=fontsize_legendtitle),
      legend.text = if (is.null(fontsize_legendtext)) element_blank() else element_text(size = fontsize_legendtext),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      plot.background=element_blank(),
      aspect.ratio = aspect.ratio) +

    labs(x=var_ident, y="proportion", fill = var_group)

  if (do_plot) p

  return(p)
}


