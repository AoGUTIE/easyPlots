#' Create and Arrange Ontology Enrichment Barplots
#'
#' This function generates barplots from a list of data frames containing ontology enrichment results.
#' Each data frame is plotted with customized font size, margins, and a color gradient. The plots are combined
#' vertically into a single figure.
#'
#' @param data A list of data frames. Each data frame must include the columns: `Position`, `Count`, `pvalue`, and `title`.
#' @param fontSize Numeric; font size for the y-axis text labels.
#' @param margin_left Numeric vector; left margins for each plot.
#' @param margin_right Numeric vector; right margins for each plot.
#' @param margin_top Numeric vector; top margins for each plot.
#' @param margin_bottom Numeric vector; bottom margins for each plot.
#' @param color Character vector of length 2; colors defining the gradient for the `pvalue` fill (e.g., `c("blue", "red")`).
#'
#' @return A single ggpubr object combining the individual barplots arranged in one column.
#'
#' @details
#' The function iterates through each data frame in `data`, generating a horizontal barplot for each one using `ggplot2`,
#' and adjusting the margins individually according to the corresponding values in `margin_left`, `margin_right`,
#' `margin_top`, and `margin_bottom`. The plots are then combined vertically using `ggpubr::ggarrange()`.
#'
#' @examples
#' \dontrun{
#' final_plot <- ontologyPlot(
#'   data = list(df1, df2, df3),
#'   fontSize = 12,
#'   margin_left = c(5, 5, 5),
#'   margin_right = c(5, 5, 5),
#'   margin_top = c(5, 5, 5),
#'   margin_bottom = c(5, 5, 5),
#'   color = c("blue", "red")
#' )
#'
#' print(final_plot)
#' }
#'
#' @export
ontologyPlot <- function(data,
                         fontSize,
                         margin_left,
                         margin_right,
                         margin_top,
                         margin_bottom,
                         color) {

  plots <- list()

  for (i in seq_along(data)) {

    # márgenes para este plot
    mL <- margin_left[i]
    mR <- margin_right[i]
    mT <- margin_top[i]
    mB <- margin_bottom[i]

    ggplotbasic <- ggplot2::ggplot(data[[i]], ggplot2::aes(x = Position, y = Count, fill = pvalue)) +
      ggplot2::geom_bar(stat = "identity", width = 0.8) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = fontSize),
        strip.background = ggplot2::element_rect(colour = "black"),
        strip.text.y = ggplot2::element_text(angle = 0),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.key.size = grid::unit(0.5, "line"),
        legend.margin = ggplot2::margin(1, 1, 1, 1),
        legend.box.background = ggplot2::element_rect(color = "black", size = 1),
        plot.margin = ggplot2::margin(mT, mR, mB, mL)
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::facet_grid(title ~ .) +
      ggplot2::scale_fill_gradientn(colors = c(color[1], color[2]), oob = scales::squish)

    plots[[i]] <- ggplotbasic
  }

  # Combina los gráficos
  finalPlot <- ggpubr::ggarrange(plotlist = plots, ncol = 1, nrow = length(plots))

  return(finalPlot)
}
