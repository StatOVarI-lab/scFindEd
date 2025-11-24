#' @keywords internal
#' Converts a list of barcode distances into a data frame, better suited for
#' building the downstream QC ggplot histogram
prepare_hamming_plot_dataframe <- function(
  BC_list,
  BC_name,
  verbose
) {
  if (verbose) cli::cli_alert_info("Preparing data frame for {BC_name}...")

  # Distance from nearest whitelist barcode
  closest <- as.numeric(unlist(lapply(BC_list, function(x) x[1])))

  # Distance from second nearest whitelist barcode
  second_closest <- as.numeric(unlist(lapply(BC_list, function(x) x[2])))

  # Combine into a data frame
  BC_df <- data.frame(
    distance = c(closest, second_closest),
    type = factor(c(rep("Closest", length(closest)),
                    rep("Second Closest", length(second_closest))))
  )

  if (verbose) {
    cli::cli_alert_success("Successfully created data frame for {BC_name}.")
  }

  return(BC_df)
}


#' @keywords internal
#' Produces a histogram of Hamming distances to the nearest whitelist barcode.
barcode_hamming_plot <- function(
  data,
  title
) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = distance, fill = type)) +
    ggplot2::geom_histogram(position = "dodge", binwidth = 0.5, color = "black", alpha = 0.8) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = title,
      subtitle = "Distribution of distances by barcode match type",
      x = "Hamming Distance",
      y = "Number of Cells",
      fill = "Barcode Match"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", size = 0.3),
      panel.grid.minor.y = ggplot2::element_line(color = "grey95", size = 0.2)
    )

  return(plot)
}
