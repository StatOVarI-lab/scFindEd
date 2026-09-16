#' Assign Identity Barcodes
#'
#' This function begins by correcting identity barcode combinations by applying
#' Hamming distance-based threshold correction. It then assigns cell-identity
#' barcode combinations as follows: an identity barcode combination is assigned
#' to a cell if and only if there are at least `min_reads` reads of that barcode
#' combination in the cell, and they constitute more than `min_reads_proportion`
#' of the barcode reads in that cell.
#'
#' @inheritParams global_params source_dir sample_name output_dir verbose
#'
#' @param BC_thresholds A named list pairing each identity barcode ID with the
#' barcode-specific distance threshold below which an extracted barcode may be
#' corrected to a whitelist barcode.
#' @param reads_thresholds Named list containing the following barcode
#' correction thresholds: `min_reads`, `min_reads_proportion`.
#'
#' @return A filtered data frame of barcode assignments that meet the read count
#' and proportion thresholds.
#'
#' @export
assign_identity_barcodes <- function(
  source_dir,
  sample_name,
  BC_thresholds,
  reads_thresholds,
  output_dir,
  verbose
) {
  # -- Import cell to extracted BCs --

  cell_to_extracted_BCs_file_name <- paste0(sample_name, "_cell_to_extracted_BCs.csv")
  cell_to_extracted_BCs <- import_cell_to_extracted_BCs(source_dir,
                                                        cell_to_extracted_BCs_file_name,
                                                        verbose)

  # Ensure each barcode in cell_to_extracted_BCs has an associated BC_threshold
  if (!is.list(BC_thresholds)) {
    cli::cli_abort("BC_thresholds must be a named list of barcode Hamming distance thresholds.")
  }
  extracted_BC_IDs = setdiff(colnames(cell_to_extracted_BCs), c("cell_barcode", "n"))
  threshold_BC_IDs = names(BC_thresholds)

  if (!identical(sort(extracted_BC_IDs), sort(threshold_BC_IDs))) {
    cli::cli_abort('The key names for BC_thresholds must exactly match the barcode IDs in "cell_to_extracted_BCs.csv".')
  }

  # -- Import barcode distance lists --

  barcode_distance_lists <- import_barcode_distance_lists(source_dir,
                                                          sample_name,
                                                          extracted_BC_IDs,
                                                          verbose)

  # -- Barcode Correction --

  cell_to_corrected_BCs <- correct_identity_barcodes(
    cell_to_extracted_BCs,
    barcode_distance_lists,
    BC_thresholds,
    verbose
  )

  # -- Identity Barcode Assignment --

  cell_to_assigned_BCs <- assign_identity_barcodes_to_cells(cell_to_corrected_BCs,
                                                            verbose)

  # -- Filtration --

  barcode_system <- filter_assigned_barcodes(cell_to_assigned_BCs,
                                             reads_thresholds,
                                             verbose)

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    barcode_system_file_name <- paste0(sample_name, "_barcode_system.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(barcode_system,
                     file = file.path(output_dir, barcode_system_file_name))

    cli::cli_alert_success(
      "Cell to identity barcode system saved to: {normalizePath(output_dir)}."
    )
  }
}


#' This function corrects barcode assignments by applying Hamming distance-based
#' threshold correction.
#' @keywords internal
correct_identity_barcodes <- function(
  cell_to_extracted_BCs,
  barcode_distance_lists,
  BC_thresholds,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Correcting Identity barcodes...",
                                      msg_done = "Identity barcodes corrected.")

  # Initialise corrected barcodes data frame
  cell_to_corrected_BCs <- data.frame(cell_barcode = cell_to_extracted_BCs$cell_barcode)

  for (BC_ID in names(BC_thresholds)) {
    # Extract correction threshold and distance list from dynamic lists
    correction_threshold <- BC_thresholds[[BC_ID]]
    distance_list <- barcode_distance_lists[[BC_ID]]

    # Dynamically construct data frame with corrected barcode IDs
    df_correction <- setNames(
      data.frame(
        names(distance_list),
        unlist(lapply(distance_list, correct_BC, threshold = correction_threshold))
      ), c(
        paste0(BC_ID, "_extracted"),
        paste0(BC_ID, "_corrected_ID")
      ))

    # Map corrected barcodes back to each cell by matching IDs
    cell_to_corrected_BCs[[BC_ID]] <- df_correction[[paste0(BC_ID, "_corrected_ID")]][
      match(cell_to_extracted_BCs[[BC_ID]],
            df_correction[[paste0(BC_ID, "_extracted")]])]
  }

  # Re-append "n" column from extracted barcodes data frame
  cell_to_corrected_BCs[["n"]] <- cell_to_extracted_BCs$n

  # Filter out rows with missing barcode corrections
  for (BC_ID in names(BC_thresholds)) {
    cell_to_corrected_BCs <- subset(cell_to_corrected_BCs, get(BC_ID) != "")
  }

  # Group by full barcode combination and sum the counts
  group_by_cols <- setdiff(colnames(cell_to_corrected_BCs), "n")
  cell_to_corrected_BCs <- cell_to_corrected_BCs |>
    dplyr::group_by(dplyr::pick(dplyr::all_of(group_by_cols))) |>
    dplyr::summarise(n_corrected = sum(n), .groups = "drop")

  if (verbose) {
    cli::cli_process_done()

    # Calculate pre-correction cells and reads total
    pre_correction_cells <- length(unique(cell_to_extracted_BCs[["cell_barcode"]]))
    pre_correction_reads <- sum(cell_to_extracted_BCs$n)
    # Calculate post-correction cells and reads total
    post_correction_cells <- length(unique(cell_to_corrected_BCs[["cell_barcode"]]))
    post_correction_reads <- sum(cell_to_corrected_BCs$n_corrected)

    cli::cli_alert_info("Cells retained: {post_correction_cells} / {pre_correction_cells}.")
    cli::cli_alert_info("Total read pairs retained: {post_correction_reads} / {pre_correction_reads}.")
  }

  return(cell_to_corrected_BCs)
}


#' @keywords internal
assign_identity_barcodes_to_cells <- function(
  cell_to_corrected_BCs,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Assigning corrected identity barcodes...",
                                      msg_done = "Corrected identity barcodes assigned to cells.")

  # For each cell barcode, calculate the proportion of reads for each identity barcode
  # assignment within the cell
  cell_to_corrected_BCs <- cell_to_corrected_BCs |>
    dplyr::group_by(cell_barcode) |>
    dplyr::mutate(proportion_reads_in_cell = n_corrected / sum(n_corrected, na.rm = TRUE)) |>
    dplyr::ungroup()

  if (verbose) cli::cli_process_done()

  return (cell_to_corrected_BCs)
}


#' @keywords internal
filter_assigned_barcodes <- function(
  cell_to_assigned_BCs,
  reads_thresholds,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Filtering assigned barcodes...",
                                      msg_done = "Assigned barcodes filtered.")

  # Filters barcode assignments by read count, read count proportion thresholds
  cell_to_filtered_BCs <- cell_to_assigned_BCs |>
    dplyr::filter(
      n_corrected >= reads_thresholds$min_reads &
        proportion_reads_in_cell >= reads_thresholds$min_reads_proportion)

  # Creates a combined barcode string for each assignment
  cell_to_filtered_BCs <- cell_to_filtered_BCs |>
    tidyr::unite(col = "all_BCs", tidyr::starts_with("BC"), sep = "_", remove = FALSE)

  if (verbose) {
    cli::cli_process_done()

    # Calculate pre-filtration cells and reads total
    pre_filtration_cells <- length(unique(cell_to_assigned_BCs[["cell_barcode"]]))
    pre_filtration_reads <- sum(cell_to_assigned_BCs$n_corrected)
    # Calculate post-filtration cells and reads total
    post_filtration_cells <- length(unique(cell_to_filtered_BCs[["cell_barcode"]]))
    post_filtration_reads <- sum(cell_to_filtered_BCs$n_corrected)

    cli::cli_alert_info("Cells retained: {post_filtration_cells} / {pre_filtration_cells}.")
    cli::cli_alert_info("Total read pairs retained: {post_filtration_reads} / {pre_filtration_reads}.")
  }

  return(cell_to_filtered_BCs)
}
