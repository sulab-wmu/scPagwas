
#' heritability_cor_scatterplot
#'
#' @param gene_heri_cor
#' @param topn_genes_label
#' @param color_low low color
#' @param color_high high color
#' @param color_mid mid color
#' @param text_size text size
#' @param do_plot Whether to plot, logical
#' @param figurenames The filename and address of the output plot
#' @param width figure width
#' @param height figure height
#'
#' @return
#' @export
#'
#' @examples
heritability_cor_scatterplot <- function(gene_heri_cor = Pagwas$gene_heritability_correlation,
                                         topn_genes_label = 10,
                                         color_low = "#035397",
                                         color_high = "#F32424",
                                         color_mid = "white",
                                         text_size = 3,
                                         do_plot = T,
                                         figurenames = NULL,
                                         width = 7,
                                         height = 7) {
  cor_df <- data.frame(genes = rownames(gene_heri_cor), cor = gene_heri_cor[, 1], text = rep(NA, nrow(gene_heri_cor)))

  cor_df <- cor_df[order(cor_df$cor, decreasing = T), ]
  cor_df$text[1:topn_genes_label] <- cor_df$genes[1:topn_genes_label]
  # save(stage_df1,stage_df2,file="stage_df.RData")
  cor_df$order <- 1:nrow(cor_df)

  p <- ggplot(data = cor_df) +
    geom_point(mapping = aes(x = order, y = cor, color = cor)) +
    scale_colour_gradient2(
      low = color_low,
      mid = color_mid,
      high = color_high,
      midpoint = 0
    ) +
    theme_classic() +
    ggrepel::geom_text_repel(aes(x = order, y = cor, label = text), na.rm = T, size = text_size)

  if (do_plot) print(p)
  if (!is.null(figurenames)) {
    pdf(file = figurenames, width = width, height = height)
    print(p)
    dev.off()
  }
}
