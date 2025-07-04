#' Generate a Volcano Plot for Differential Expression Results
#'
#' This function generates a customizable volcano plot from differential expression analysis results,
#' highlighting significant genes based on fold-change and statistical criteria, and optionally labeling
#' top or user-specified genes.
#'
#' @param data A data frame containing at least `log2FoldChange`, `pvalue`, `stat`, `Description`, and `SYMBOL` columns.
#' @param FC Numeric; log2 fold-change threshold to define up- and down-regulated genes.
#' @param statSig Character; which statistic to use: `"t-test"` or `"p-value"`.
#' @param PV Numeric; p-value threshold to define significance (used if `statSig = "p-value"`).
#' @param STAT Numeric; statistic threshold to define significance (used if `statSig = "t-test"`).
#' @param Description Logical (`TRUE` or `FALSE`); whether to filter rows based on the `Description` column.
#' @param Filtering Character; value in the `Description` column to filter by if `Description = TRUE`.
#' @param maxGenes Integer; maximum number of top up- and down-regulated genes to label.
#' @param TopGenes Logical (`TRUE` or `FALSE`); whether to highlight the top `maxGenes` genes based on significance.
#' @param List Logical (`TRUE` or `FALSE`); whether to highlight a specific list of genes.
#' @param GeneList Character vector; list of gene symbols to highlight if `List = TRUE`.
#' @param colors1 Character vector of length 3; colors for `"Up Regulated"`, `"Down Regulated"`, and `"Not Significant"`.
#' @param colors2 Character; color used for axes, lines, and top gene labels.
#' @param colors3 Character; color used for the user-specified gene list.
#' @param Ymax,Ymin Numeric; y-axis limits.
#' @param Xmax,Xmin Numeric; x-axis limits.
#'
#' @return A `ggplot2` object representing the volcano plot.
#'
#' @details
#' The function classifies genes into `"Up Regulated"`, `"Down Regulated"`, or `"Not Significant"` based on `FC`, `PV`, and/or `STAT`.
#' Top significant genes and/or user-specified genes (`GeneList`) can be highlighted and labeled.
#' Two modes of statistical significance are supported: `"t-test"` (using `STAT`) and `"p-value"` (using `PV`).
#'
#' @examples
#' \dontrun{
#' plot <- volcanoPlot(
#'   data = my_data,
#'   FC = 1,
#'   statSig = "p-value",
#'   PV = 0.05,
#'   STAT = 2,
#'   Description = TRUE,
#'   Filtering = "some category",
#'   maxGenes = 10,
#'   TopGenes = TRUE,
#'   List = TRUE,
#'   GeneList = c("TP53", "BRCA1", "EGFR"),
#'   colors1 = c("red", "blue", "gray"),
#'   colors2 = "black",
#'   colors3 = "green",
#'   Ymax = 10,
#'   Ymin = 0,
#'   Xmax = 5,
#'   Xmin = -5
#' )
#' print(plot)
#' }
#'
#' @export
volcanoPlot <- function(data,
                        FC,
                        statSig,
                        PV,
                        STAT,
                        Description = TRUE,
                        Filtering,
                        maxGenes,
                        TopGenes = TRUE,
                        List = TRUE,
                        GeneList,
                        colors1,
                        colors2,
                        colors3,
                        Ymax,
                        Xmax,
                        Ymin,
                        Xmin) {

  custom_colors1 <- c("Up Regulated" = colors1[1], "Down Regulated" = colors1[2], "Not Significant" = colors1[3])
  custom_colors2 <- colors2
  custom_colors3 <- colors3

  if (statSig == "t-test") {
    data$ABSTAT <- abs(data$stat)
    data$diffexpressed <- "Not Significant"
    data$diffexpressed[data$log2FoldChange > FC & data$ABSTAT > STAT] <- "Up Regulated"
    data$diffexpressed[data$log2FoldChange < -FC & data$ABSTAT > STAT] <- "Down Regulated"

    NoNAs <- if (Description) data[data$Description == Filtering, ] else data

    if (TopGenes) {
      downregulated_genes <- NoNAs[NoNAs$log2FoldChange < -FC & NoNAs$ABSTAT > STAT, ]
      upregulated_genes <- NoNAs[NoNAs$log2FoldChange > FC & NoNAs$ABSTAT > STAT, ]
      top10_down <- utils::head(downregulated_genes[order(downregulated_genes$ABSTAT, decreasing = TRUE), ], maxGenes)
      top10_up <- utils::head(upregulated_genes[order(upregulated_genes$ABSTAT, decreasing = TRUE), ], maxGenes)
    }

    if (List) {
      NoNAs2 <- NoNAs[NoNAs$SYMBOL %in% GeneList, ]
      genes.down <- NoNAs2[NoNAs2$log2FoldChange < -FC & NoNAs2$ABSTAT > STAT, ]
      genes.up <- NoNAs2[NoNAs2$log2FoldChange > FC & NoNAs2$ABSTAT > STAT, ]
    }

    volcanoPlot <- ggplot2::ggplot(NoNAs, ggplot2::aes(x = log2FoldChange, y = ABSTAT, color = diffexpressed)) +
      ggplot2::geom_hline(yintercept = -log10(0.05), col = custom_colors2, linetype = 'dashed') +
      ggplot2::geom_point(size = 4, alpha = 0.5) +
      ggplot2::scale_y_continuous(limits = c(Ymin, Ymax)) +
      ggplot2::scale_color_manual(values = custom_colors1) +
      ggplot2::coord_cartesian(xlim = c(-Xmax, Xmax)) +
      ggplot2::labs(color = 'Different Expression',
                    x = expression("log"[2]*" Fold Change"),
                    y = expression("Wald Test")) +
      ggplot2::scale_x_continuous(breaks = seq(-Xmax, Xmax, 2)) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white"),
        panel.grid.major = ggplot2::element_line(colour = "grey90"),
        panel.grid.minor = ggplot2::element_line(colour = "grey90"),
        axis.line = ggplot2::element_line(colour = custom_colors2),
        axis.text = ggplot2::element_text(colour = custom_colors2, size = 20, face = "bold"),
        axis.title = ggplot2::element_text(colour = custom_colors2, size = 20, face = "bold"),
        legend.text = ggplot2::element_text(size = 16, face = "bold"),
        legend.title = ggplot2::element_text(size = 16, face = "bold")
      ) +
      ggplot2::geom_vline(xintercept = c(-1, 1), col = custom_colors2, linetype = 'dashed', size = 1)

    if (TopGenes) {
      volcanoPlot <- volcanoPlot +
        ggplot2::geom_point(data = top10_up, ggplot2::aes(x = log2FoldChange, y = ABSTAT), color = custom_colors2, size = 6, shape = 21) +
        ggrepel::geom_label_repel(data = top10_up, ggplot2::aes(label = SYMBOL), color = custom_colors2,
                                  label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
                                  size = 6, xlim = c((Xmin + 1), Xmax), force = 10, direction = "both") +
        ggplot2::geom_point(data = top10_down, ggplot2::aes(x = log2FoldChange, y = ABSTAT), color = custom_colors2, size = 6, shape = 21) +
        ggrepel::geom_label_repel(data = top10_down, ggplot2::aes(label = SYMBOL), color = custom_colors2,
                                  label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
                                  size = 6, xlim = c(-(Xmin + 1), -Xmax), force = 10, direction = "both")
    }

    if (List) {
      volcanoPlot <- volcanoPlot +
        ggplot2::geom_point(data = genes.up, ggplot2::aes(x = log2FoldChange, y = ABSTAT), color = custom_colors3, size = 6, shape = 21) +
        ggrepel::geom_label_repel(data = genes.up, ggplot2::aes(label = SYMBOL), color = custom_colors3,
                                  label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
                                  size = 6, xlim = c((Xmin + 1), Xmax), force = 10, direction = "both") +
        ggplot2::geom_point(data = genes.down, ggplot2::aes(x = log2FoldChange, y = ABSTAT), color = custom_colors3, size = 6, shape = 21) +
        ggrepel::geom_label_repel(data = genes.down, ggplot2::aes(label = SYMBOL), color = custom_colors3,
                                  label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
                                  size = 6, xlim = c(-(Xmin + 1), -Xmax), force = 10, direction = "both")
    }

  } else if (statSig == "p-value") {

    data$diffexpressed <- "Not Significant"
    data$diffexpressed[data$log2FoldChange > FC & data$pvalue < PV]  <- "Up Regulated"
    data$diffexpressed[data$log2FoldChange < -FC & data$pvalue < PV] <- "Down Regulated"

    NoNAs <- if (Description) {
      data[data$Description == Filtering, ]
    } else {
      data
    }

    if (TopGenes) {
      NoNAs$neg_log10_pvalue <- -log10(NoNAs$pvalue)
      downregulated_genes <- NoNAs[NoNAs$log2FoldChange < -FC & NoNAs$pvalue < PV, ]
      upregulated_genes   <- NoNAs[NoNAs$log2FoldChange >  FC & NoNAs$pvalue < PV, ]
      top10_down <- utils::head(
        downregulated_genes[order(downregulated_genes$neg_log10_pvalue, decreasing = TRUE), ],
        maxGenes
      )
      top10_up <- utils::head(
        upregulated_genes[order(upregulated_genes$neg_log10_pvalue, decreasing = TRUE), ],
        maxGenes
      )
    }

    if (List) {
      NoNAs2 <- NoNAs[NoNAs$SYMBOL %in% GeneList, ]
      genes.down <- NoNAs2[NoNAs2$log2FoldChange < -FC & NoNAs2$pvalue < PV, ]
      genes.up   <- NoNAs2[NoNAs2$log2FoldChange >  FC & NoNAs2$pvalue < PV, ]
    }

    volcanoPlot <- ggplot2::ggplot(
      data = NoNAs,
      ggplot2::aes(x = log2FoldChange, y = -log10(pvalue), color = diffexpressed)
    ) +
      ggplot2::geom_hline(yintercept = -log10(0.05), col = custom_colors2, linetype = 'dashed') +
      ggplot2::geom_point(size = 4, alpha = 0.5) +
      ggplot2::scale_y_continuous(limits = c(Ymin, Ymax)) +
      ggplot2::scale_color_manual(values = custom_colors1) +
      ggplot2::coord_cartesian(xlim = c(-Xmax, Xmax)) +
      ggplot2::labs(
        color = 'Different Expression',
        x = expression("log"[2]*" Fold Change"),
        y = expression("-log"[10]*" p-value")
      ) +
      ggplot2::scale_x_continuous(breaks = seq(-Xmax, Xmax, 2)) +
      ggplot2::theme(
        panel.background  = ggplot2::element_rect(fill = "white"),
        panel.grid.major  = ggplot2::element_line(colour = "grey90"),
        panel.grid.minor  = ggplot2::element_line(colour = "grey90"),
        axis.line          = ggplot2::element_line(colour = custom_colors2),
        axis.text          = ggplot2::element_text(colour = custom_colors2, size = 20, face = "bold"),
        axis.title         = ggplot2::element_text(colour = custom_colors2, size = 20, face = "bold"),
        legend.text        = ggplot2::element_text(size = 16, face = "bold"),
        legend.title       = ggplot2::element_text(size = 16, face = "bold")
      ) +
      ggplot2::geom_vline(xintercept = c(-1, 1), col = custom_colors2, linetype = 'dashed', size = 1)

    if (TopGenes) {
      volcanoPlot <- volcanoPlot +
        ggplot2::geom_point(data = top10_up, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue)),
                            color = custom_colors2, size = 6, shape = 21) +
        ggrepel::geom_label_repel(
          data = top10_up, ggplot2::aes(label = SYMBOL), color = custom_colors2,
          label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
          size = 6, xlim = c((Xmin + 1), Xmax), force = 10, direction = "both"
        ) +
        ggplot2::geom_point(data = top10_down, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue)),
                            color = custom_colors2, size = 6, shape = 21) +
        ggrepel::geom_label_repel(
          data = top10_down, ggplot2::aes(label = SYMBOL), color = custom_colors2,
          label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
          size = 6, xlim = c(-(Xmin + 1), -Xmax), force = 10, direction = "both"
        )
    }

    if (List) {
      volcanoPlot <- volcanoPlot +
        ggplot2::geom_point(data = genes.up, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue)),
                            color = custom_colors3, size = 6, shape = 21) +
        ggrepel::geom_label_repel(
          data = genes.up, ggplot2::aes(label = SYMBOL), color = custom_colors3,
          label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
          size = 6, xlim = c((Xmin + 1), Xmax), force = 10, direction = "both"
        ) +
        ggplot2::geom_point(data = genes.down, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue)),
                            color = custom_colors3, size = 6, shape = 21) +
        ggrepel::geom_label_repel(
          data = genes.down, ggplot2::aes(label = SYMBOL), color = custom_colors3,
          label.padding = grid::unit(0.2, "lines"), box.padding = grid::unit(0.1, "lines"),
          size = 6, xlim = c(-(Xmin + 1), -Xmax), force = 10, direction = "both"
        )
    }
  }

  return(volcanoPlot)
}
