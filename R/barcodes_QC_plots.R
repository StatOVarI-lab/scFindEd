#' Barcoding Quality Control Plots
#'
#' This function produces illustrative quality control plots displaying the
#' distribution of Hamming distances for the closest and second closest barcode
#' matches for both barcode 1 and barcode 2-3.
#'
#' @inheritParams global_params source_dir output_dir verbose
#'
#' @param list_BC1_file_name String containing the exact file name (including the
#' file type) for the list of closest matches for barcode 1, with Hamming distances.
#' @param list_BC2_3_file_name String containing the exact file name (including
#' the file type) for the list of closest matches for barcode 2-3, with Hamming
#' distances.
#'
#' @return Two sanity-check histograms that compare the Hamming distance of the
#' closest barcode matches with the Hamming distance of the second closest
#' barcode matches.
#'
#' @export
barcodes_QC_plots <- function(
  source_dir,
  list_BC1_file_name,
  list_BC2_3_file_name,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import Datasets --

  invisible(list2env(import_barcode_distance_lists(
    source_dir,
    list_BC1_file_name,
    list_BC2_3_file_name,
    verbose
  ), envir = environment()))

  # -- BC1 QC Plot --

  if (verbose) cli::cli_process_start(
    msg = "Constructing QC plots for BC1...",
    msg_done = "QC plot for BC1 exported to {normalizePath(output_dir)}.")

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

    if (verbose) cli::cli_process_done()
  }

  # -- BC2/BC3 QC Plot --

  if (verbose) cli::cli_process_start(
    msg = "Constructing QC plots for BC2_3 pairs...",
    msg_done = "QC plot for BC2_3 pairs exported to {normalizePath(output_dir)}.")

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

    if (verbose) cli::cli_process_done()
  }
}



