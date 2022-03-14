
#' Bootstrap_P_Barplot
#' @description This bar plot shows the -log2(p value) for bootstrap result,
#' using the ggplot packages
#'
#' @param Pagwas Pagwas format of result in Pagwas_main()
#' @param title The title names of the plot
#' @param figurenames The filename and address of the output plot,
#' default is "test_barplot.pdf".IF figurenames= NULL, only plot the figure
#' and have not pdf figure.
#' @param width figure width, default is 5
#' @param height figure height,default is 7
#'
#' @return A figure of barplot in pdf format, red color is significant.
#' @export
#'
#' @examples
#' library(scPagwas)
#' # Pagwas is the result of Pagwas_main()
#' Bootstrap_P_Barplot(Pagwas=Pagwas,
#'                     title = "Test scPagwas",
#'                     figurenames = "test_barplot.pdf",
#'                     width = 5, height = 7)
Bootstrap_P_Barplot <- function(Pagwas,
                                title = "Test scPagwas",
                                figurenames = NULL,
                                width = 5, height = 7) {
  cell_severe <- Pagwas$bootstrap_results[-1, ]
  cell_severe$logp <- -log2(cell_severe$bp_value)
  cell_severe$sig <- rep("b", nrow(cell_severe))
  cell_severe$sig[which(cell_severe$bp_value < 0.05)] <- "a"

  if (sum(cell_severe$bp_value < 0.05) > 0) {
    p1 <- ggplot2::ggplot(cell_severe, aes(x = reorder(annotation, logp),
                                           y = logp,
                                           fill = sig)) +
      geom_bar(position = "dodge", stat = "identity") +
      theme_classic() +
      labs(x = "", y = "-log2(p)", title = title) +
      coord_flip() +
      # scale_fill_discrete()+
      geom_hline(aes(yintercept = 4.321),
                 colour = "#990000",
                 linetype = "dashed") +
      scale_fill_manual(values = c("#BB6464", "#C3DBD9")) +
      theme(legend.position = "none")
  } else {
    p1 <- ggplot2::ggplot(cell_severe,
                          aes(x = reorder(annotation, logp),
                              y = logp)) +
      geom_bar(position = "dodge", stat = "identity", color = "#C3DBD9") +
      theme_classic() +
      labs(x = "", y = "-log2(p)", title = title) +
      coord_flip() +
      theme(legend.position = "none")
    # scale_fill_discrete()+
    # geom_hline(aes(yintercept=4.321), colour="#990000", linetype="dashed")+
  }
  print(p1)
  if(!is.null(figurenames)){
    pdf(figurenames, width = width, height = height)
    print(p1)
    dev.off()
  }

 }


#' Bootstrap_estimate_Plot
#' @description This forest plot shows the correct estimate values and
#' 95% CI for different celltyppes, using the ggplot packages
#' @param Pagwas Pagwas format of result from Pagwas_main()
#' @param figurenames The filename and address of the output plot,
#' default is "test_forest.pdf".IF figurenames= NULL,
#' only plot the figure and have not pdf figure.
#' @param width figure width
#' @param height figure height
#'
#' @return A forest plot with the table of p values
#' @export
#'
#' @examples
#' library(scPagwas)
#' # Pagwas is the result of Pagwas_main()
#' Bootstrap_estimate_Plot(Pagwas=Pagwas,
#'                         figurenames = "test_forest.pdf",
#'                         width = 9, height = 7)
Bootstrap_estimate_Plot <- function(Pagwas,
                                    figurenames = NULL,
                                    width = 9,
                                    height = 7) {
  bootstrap_results <- Pagwas$bootstrap_results[-1, c("bp_value",
                                                      "bias_corrected_estimate",
                                                      "CI_lo",
                                                      "CI_hi")]
  bootstrap_results <- bootstrap_results[order(bootstrap_results$bias_corrected_estimate, decreasing = T), ]


  dat <- data.frame(
    Index = seq_len(nrow(bootstrap_results)), ## This provides an order to the data
    label = rownames(bootstrap_results),
    estimate = bootstrap_results$bias_corrected_estimate,
    lower = bootstrap_results$CI_lo,
    upper = bootstrap_results$CI_hi,
    Pvalue = bootstrap_results$bp_value
  )

  plot1 <- ggplot2::ggplot(dat, aes(y = Index,
                                    x = estimate)) +
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper),
                   color = "#6D8299",
                   height = 0.25) +
    geom_point(shape = 18,
               size = 5,
               color = "#D57E7E") +
    geom_vline(xintercept = 0,
               color = "#444444",
               linetype = "dashed",
               cex = 1,
               alpha = 0.5) +
    scale_y_continuous(name = "", breaks = seq_len(nrow(dat)),
                       labels = dat$label,
                       trans = "reverse") +
    xlab("bias_corrected_estimate (95% CI)") +
    ylab(" ") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x.bottom = element_text(size = 12, colour = "black"),
      axis.title.x = element_text(size = 12, colour = "black")
    )
  # plot1

  table_base <- ggplot2::ggplot(dat, aes(y = label)) +
    ylab(NULL) +
    xlab("  ") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(color = "white",
                                 hjust = -3,
                                 size = 25), ## This is used to help with alignment
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )

  ## OR point estimate table
  tab1 <- table_base +
    labs(title = "space") +
    geom_text(aes(y = rev(Index), x = 1,
                  label = sprintf("%.3f", round(Pvalue, digits = 3))),
              size = 4) + ## decimal places
    ggtitle("Pvalue")

  lay <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3), nrow = 1)
  plot2 <- gridExtra::grid.arrange(plot1, tab1, layout_matrix = lay)
  print(plot2)
  ## save the pdf figure
  if(!is.null(figurenames)){
  pdf(file = figurenames, width = width, height = height)
  print(plot2)
  dev.off()
  }
}
