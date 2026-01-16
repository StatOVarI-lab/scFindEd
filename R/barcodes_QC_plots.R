#' Barcoding QC Plots
#'
#' This function produces illustrative QC plots displaying the distribution of
#' Hamming distances for the closest and second closest barcode matches for both
#' BC1 and BC2/3.
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param list_BC1_file_name String containing the file name structure for the
#' list of closest barcode matches for BC1, with Hamming distances.
#' @param list_BC2_3_file_name String containing the file name structure for the
#' list of closest barcode matches for BC2/3, with Hamming distances.
#' @param output_dir String containing the relative path to the directory that
#' the plots should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return Two sanity-check histograms that compare the Hamming distance of the
#' closest barcode matches with the Hamming distance of the second closest
#' barcode matches.

#' @examples
#' barcodes_QC_plots(
#' source_dir = "data",
#' list_BC1_file_name = "list_BC1.rds",
#' list_BC2_3_file_name = "list_BC2_3.rds",
#' output_dir = "plots",
#' verbose = TRUE
#' )

#' @importFrom dplyr %>%

#' @export
barcodes_QC_plots <- function(
  source_dir,
  list_BC1_file_name,
  list_BC2_3_file_name,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import Datasets --

  imported_data <- import_barcode_distance_lists(
    source_dir,
    list_BC1_file_name,
    list_BC2_3_file_name,
    verbose
  )

  invisible(list2env(imported_data, envir = environment()))

  cli::cli_alert_success("Successfully imported datasets.")

  # -- BC1 QC Plot --

  if (verbose) cli::cli_inform("Constructing QC plots for BC1...")

  # Prepare dataframe for BC1
  BC1_df <- prepare_hamming_plot_dataframe(list_BC1, "BC1", verbose)

  # Plot histogram of BC1 Hamming distances onto a single axis
  BC1_hamming_plot <-
    barcode_hamming_plot(BC1_df, title = "BC1: Hamming Distance Comparison")

  if (!is.null(output_dir)) {
    # Create directory if provided directory name does not exist
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Export BC1_hamming_plot
    ggplot2::ggsave(file.path(output_dir, "BC1_hamming_plot.png"),
           plot = BC1_hamming_plot,
           width = 8, height = 5, dpi = 300)

    cli::cli_alert_success(
      "QC plot for BC1 exported to {normalizePath(output_dir)}")
  }

  # -- BC2/BC3 QC Plot --

  if (verbose) cli::cli_inform("Constructing QC plots for BC2/3...")

  # Prepare dataframe for BC1
  BC2_3_df <- prepare_hamming_plot_dataframe(list_BC2_3, "BC2/3", verbose)

  # Plot histogram of BC2/3 Hamming distances onto a single axis
  BC2_3_hamming_plot <-
    barcode_hamming_plot(BC2_3_df, title = "BC2/3: Hamming Distance Comparison")

  if (!is.null(output_dir)) {
    # Create directory if provided directory name does not exist
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Export BC2_3_hamming_plot
    ggplot2::ggsave(file.path(output_dir, "BC2_3_hamming_plot.png"),
           plot = BC2_3_hamming_plot,
           width = 8, height = 5, dpi = 300)

    cli::cli_alert_success(
      "QC plot for BC2/3 exported to {normalizePath(output_dir)}")
  }
}



