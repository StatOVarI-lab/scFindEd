#' Assign UTR Barcodes
#'
#' This function begins by correcting UTR barcode combinations by applying
#' Hamming distance-based threshold correction. It then assigns cell-UTR barcode
#' combinations as follows: a UTR barcode combination is assigned to a cell if
#' and only if there are at least `min_reads` reads of that barcode combination
#' in the cell, and they constitute more than `min_reads_proportion` of the
#' barcode reads in that cell.
#'
#' @inheritParams global_params source_dir output_dir verbose
#'
#' @param BC_thresholds Named list containing the following barcode thresholds:
#' `BC1_threshold`, `BC2_3_threshold`.
#' @param reads_thresholds Named list containing the following barcode
#' correction thresholds: `min_reads`, `min_reads_proportion`.
#'
#' @return A filtered data frame of barcode assignments that meet the read count
#' and proportion thresholds.
#'
#' @export
assign_UTR_barcodes <- function(
  source_dir,
  BC_thresholds,
  reads_thresholds,
  output_dir,
  verbose
) {
  # -- Import Datasets --

  cell_to_extracted_BCs <- import_cell_to_extracted_BCs(
    source_dir,
    "cell_to_extracted_BCs.rds",
    verbose
  )

  invisible(list2env(import_barcode_distance_lists(
    source_dir,
    "list_BC1.rds",
    "list_BC2_3.rds",
    verbose
  ), envir = environment()))

  # -- Barcode Correction --

  cell_to_BCs_corrected <- correct_UTR_barcodes(
    cell_to_extracted_BCs,
    list_BC1,
    list_BC2_3,
    BC_thresholds,
    verbose
  )

  # -- UTR Barcode Assignment --

  cell_to_BCs_assigned <- assign_UTR_barcodes_to_cells(cell_to_BCs_corrected,
                                                       verbose)

  # -- Filtration --

  cell_to_BCs_filtered <- filter_assigned_barcodes(cell_to_BCs_assigned,
                                                   reads_thresholds,
                                                   verbose)

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    file_save_name <- file.path(output_dir, "cell_to_UTR_barcode_system.csv")

    readr::write_csv(cell_to_BCs_filtered, file_save_name)

    cli::cli_alert_success(
      "Cell to UTR barcode system saved to: {normalizePath(output_dir)}."
    )
  }
}


#' This function corrects barcode assignments by applying Hamming distance-based
#' threshold correction.
#' @keywords internal
correct_UTR_barcodes <- function(
  cell_to_extracted_BCs,
  list_BC1,
  list_BC2_3,
  BC_thresholds,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Correcting UTR barcodes...",
                                      msg_done = "UTR barcodes corrected.")

  # Extract correction thresholds
  BC1_threshold <- BC_thresholds$BC1_threshold
  BC2_3_threshold <- BC_thresholds$BC2_3_threshold

  # Correct BC1 barcodes using threshold BC1_threshold
  df_correction_BC1 <- data.frame(
    BC1_extracted = names(list_BC1),
    BC1_corrected_ID = unlist(lapply(list_BC1, correct_BC,
                                     thresh = BC1_threshold))
  )

  # Correct BC2_3 barcodes using threshold BC2_3_threshold
  df_correction_BC2_3 <- data.frame(
    BC2_3_extracted = names(list_BC2_3),
    BC2_3_corrected_ID = unlist(lapply(list_BC2_3, correct_BC,
                                       thresh = BC2_3_threshold))
  )

  # Map corrected barcodes back to each cell by matching IDs
  cell_to_BCs <- data.frame(
    cell_barcode = cell_to_extracted_BCs$cell_barcode,
    BC1 = df_correction_BC1$BC1_corrected_ID[match(
      cell_to_extracted_BCs$BC1,df_correction_BC1$BC1_extracted)],
    BC2 = df_correction_BC2_3$BC2_3_corrected_ID[match(
      cell_to_extracted_BCs$BC2,df_correction_BC2_3$BC2_3_extracted)],
    BC3 = df_correction_BC2_3$BC2_3_corrected_ID[match(
      cell_to_extracted_BCs$BC3,df_correction_BC2_3$BC2_3_extracted)],
    n = cell_to_extracted_BCs$n
  )

  # Filter out rows with missing barcode corrections
  cell_to_BCs_all_BCs <- cell_to_BCs |>
    dplyr::filter(BC1 != "", BC2 != "", BC3 != "")

  # Group by full barcode combination and sum the counts
  cell_to_BCs_corrected <- cell_to_BCs_all_BCs |>
    dplyr::group_by(cell_barcode, BC1, BC2, BC3) |>
    dplyr::summarise(n_corrected = sum(n), .groups = "drop")

  if (verbose) cli::cli_process_done()

  return(cell_to_BCs_corrected)
}


#' @keywords internal
assign_UTR_barcodes_to_cells <- function(
  cell_to_BCs_corrected,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Assigning corrected UTR barcodes...",
                                      msg_done = "Corrected UTR barcodes assigned.")

  # For each cell barcode, calculate the proportion of reads for each UTR barcode
  # assignment within the cell
  cell_to_BCs_corrected <- cell_to_BCs_corrected |>
    dplyr::group_by(cell_barcode) |>
    dplyr::mutate(proportion_reads_in_cell = n_corrected / sum(n_corrected, na.rm = TRUE)) |>
    dplyr::ungroup()

  if (verbose) cli::cli_process_done()

  return (cell_to_BCs_corrected)
}


#' @keywords internal
filter_assigned_barcodes <- function(
  cell_to_BCs_assigned,
  reads_thresholds,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Filtering assigned barcodes...",
                                      msg_done = "Assigned barcodes filtered.")

  # Filters barcode assignments by read count, read count proportion thresholds
  cell_to_BCs_filtered <- cell_to_BCs_assigned |>
    dplyr::filter(
      n_corrected >= reads_thresholds$min_reads &
        proportion_reads_in_cell >= reads_thresholds$min_reads_proportion)

  # Creates a combined barcode string for each assignment
  cell_to_BCs_filtered <- cell_to_BCs_filtered |>
    dplyr::mutate(all_BCs = paste(BC1, BC2, BC3, sep = "_"))

  if (verbose) cli::cli_process_done()

  return(cell_to_BCs_filtered)
}
