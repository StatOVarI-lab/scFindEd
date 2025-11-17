#' Barcoding QC Plots
#'
#' This function produces illustrative QC plots displaying the distribution of
#' Hamming distances for the closest and second closest barcode matches for both
#' BC1 and BC2/3.
#'
#' @param list_BC1 A list of closest barcode matches for BC1, with Hamming
#' distances.
#' @param list_BC2_3 A list of closest barcode matches for BC2/3, with Hamming
#' distances.
#' @param output_location A string detailing the path to the directory the plots
#' should be exported to.

#' @return Two sanity-check histograms that compare the Hamming distance of the
#' closest barcode matches with the Hamming distance of the second closest
#' barcode matches.

#' @examples
#' barcodes_QC_plots(
#' list_BC1 = list_BC1,
#' list_BC2_3 = list_BC2_3,
#' output_dir = "plots"
#' )

#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_classic ggsave
#' @importFrom glue glue

#' @export
barcodes_QC_plots <- function(
    list_BC1,
    list_BC2_3,
    output_dir = NULL
) {
  # -- BC1 QC Plot --

  message("Constructing QC Plots for BC1...")

  # Distance from nearest whitelist barcode
  BC1_first <- as.numeric(unlist(lapply(list_BC1, function(x) x[1])))

  # Distance from second nearest whitelist barcode
  BC1_second <- as.numeric(unlist(lapply(list_BC1, function(x) x[2])))

  # Combine into a data frame
  BC1_df <- data.frame(
    distance = c(BC1_first, BC1_second),
    type = factor(c(rep("Closest", length(BC1_first)),
                    rep("Second Closest", length(BC1_second))))
  )

  # Plot histogram of BC1 Hamming distances onto a single axis
  BC1_hamming_plot <-
    barcode_hamming_plot(BC1_df, title = "BC1: Hamming Distance Comparison")

  if (!is.null(output_dir)) {
    # Create directory if provided directory name does not exist
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Export BC1_hamming_plot
    ggsave(file.path(output_dir, "BC1_hamming_plot.png"),
           plot = BC1_hamming_plot,
           width = 8, height = 5, dpi = 300)

    message(glue::glue("QC plot for BC1 exported to {normalizePath(output_dir)}"))
  }

  # -- BC2/BC3 QC Plot --

  message("Constructing QC Plots for BC2/3...")

  # Distance from nearest whitelist barcode
  BC2_3_first <- as.numeric(unlist(lapply(list_BC2_3,function(x) x[1])))

  # Distance from second nearest whitelist barcode
  BC2_3_second <- as.numeric(unlist(lapply(list_BC2_3,function(x) x[2])))

  # Combine into a data frame
  BC2_3_df <- data.frame(
    distance = c(BC2_3_first, BC2_3_second),
    type = factor(c(rep("Closest", length(BC2_3_first)),
                    rep("Second Closest", length(BC2_3_second))))
  )

  # Plot side-by-side bars
  BC2_3_hamming_plot <-
    barcode_hamming_plot(BC2_3_df, title = "BC2/3: Hamming Distance Comparison")

  if (!is.null(output_dir)) {
    # Create directory if provided directory name does not exist
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Export BC2_3_hamming_plot
    ggsave(file.path(output_dir, "BC2_3_hamming_plot.png"),
           plot = BC2_3_hamming_plot,
           width = 8, height = 5, dpi = 300)

    message(glue::glue("QC plot for BC2/3 exported to {normalizePath(output_dir)}"))
  }
}
