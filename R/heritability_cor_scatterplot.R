
#' heritability_cor_scatterplot
#' @description Plot the scatterplot for the Genetic Expression Index(GEI) for each genes.
#' pearson methods for each gene expression and gPAS score in each cells.
#' @param gene_heri_cor Genetic Expression Index result for scPagwas.
#' @param topn_genes_label 10, the number of top genes for correlation.
#' @param color_low low color
#' @param color_high high color
#' @param color_mid mid color
#' @param text_size text size
#' @param max.overlaps mac over
#' @param do_plot parametr for "geom_label_repel", Exclude text labels that overlap too many things. Defaults to 10.
#' @param figurenames The filename and address of the output plot
#' @param width figure width
#' @param height figure height
#'
#' @return scatterplot
#' @export
#'
#' @examples
#' load(system.file("extdata", "Pagwas_data.RData", package = "scPagwas"))
#' heritability_cor_scatterplot(
#'   gene_heri_cor = Pagwas_data@misc$GeneticExpressionIndex,
#'   topn_genes_label = 10,
#'   color_low = "#035397",
#'   color_high = "#F32424",
#'   color_mid = "white",
#'   text_size = 2,
#'   do_plot = TRUE,
#'   max.overlaps = 20,
#'   width = 7,
#'   height = 7
#' )
#' @author Chunyu Deng
#' @aliases heritability_cor_scatterplot
#' @keywords heritability_cor_scatterplot, plot the scatterplot for the correlation coefficient for each genes.
heritability_cor_scatterplot <- function(gene_heri_cor = NULL,
                                         topn_genes_label = 10,
                                         color_low = "#035397",
                                         color_high = "#F32424",
                                         color_mid = "white",
                                         text_size = 3,
                                         do_plot = TRUE,
                                         max.overlaps = 10,
                                         figurenames = NULL,
                                         width = 7,
                                         height = 7) {
  text <- NULL
  cor_df <- data.frame(genes = rownames(gene_heri_cor), GEI = gene_heri_cor[, 1], text = rep(NA, nrow(gene_heri_cor)))

  cor_df <- cor_df[order(cor_df$GEI, decreasing = T), ]
  cor_df$text[1:topn_genes_label] <- cor_df$genes[1:topn_genes_label]

  cor_df$order <- seq_len(nrow(cor_df))

  p <- ggplot(data = cor_df) +
    geom_point(mapping = aes(x = order, y = GEI, color = GEI)) +
    scale_colour_gradient2(
      low = color_low,
      mid = color_mid,
      high = color_high,
      midpoint = 0
    ) +
    theme_classic() +
    ggrepel::geom_text_repel(aes(
      x = order, y = GEI,
      label = text
    ),
    max.overlaps = max.overlaps,
    na.rm = T, size = text_size
    )

  if (do_plot) print(p)
  if (!is.null(figurenames)) {
    grDevices::pdf(file = figurenames, width = width, height = height)
    print(p)
    grDevices::dev.off()
  }
}
