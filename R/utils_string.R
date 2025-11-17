#' @keywords internal
#' Produces a histogram of Hamming distances to the nearest whitelist barcode.
barcode_hamming_plot <- function(
  data,
  title
) {
  plot <- ggplot(data, aes(x = distance, fill = type)) +
    geom_histogram(position = "dodge", binwidth = 0.5, color = "black", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = title,
      subtitle = "Distribution of distances by barcode match type",
      x = "Hamming Distance",
      y = "Number of Cells",
      fill = "Barcode Match"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      panel.grid.minor.y = element_line(color = "grey95", size = 0.2)
    )

  return(plot)
}
