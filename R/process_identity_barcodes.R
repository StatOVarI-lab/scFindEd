#' Process Identity Barcodes
#'
#' This function takes barcode reads data and first performs some basic QC
#' filtering. It then extracts the identity barcodes using the optimal flank
#' info generated in `evaluate_optimal_barcode_flanks()`, before grouping cells
#' by cell_barcode and counting reads.
#'
#' @inheritParams global_params source_dir sample_name identity_AmpID min_reads output_dir verbose
#'
#' @return A table linking cells to extracted barcodes.
#'
#' @export
process_identity_barcodes <- function(
  source_dir,
  sample_name,
  identity_AmpID,
  min_reads,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import Identity Barcode Info --

  identity_barcode_info_file_name <- paste0(sample_name, "_identity_barcode_info.csv")
  identity_barcode_info <- import_identity_barcode_info(source_dir,
                                                        identity_barcode_info_file_name,
                                                        verbose)

  if (!all(identity_AmpID %in% identity_barcode_info[["AmpID"]])) {
    cli::cli_abort('All AmpIDs in identity_AmpID must be present in "_identity_barcode_info.csv".')
  }

  # Initialise list for extracted_barcodes tables
  extracted_barcodes_tables <- list()

  for (AmpID in identity_AmpID) {
    if (length(identity_AmpID) > 1 & verbose) {
      cli::cli_h2(AmpID)
    }

    # -- Import Identity Amplicon Reads --

    identity_amplicon_reads_file_name <- paste0(sample_name, "_", AmpID, "_reads.tsv.gz")
    identity_amplicon_reads <- import_identity_amplicon_reads(source_dir,
                                                              identity_amplicon_reads_file_name,
                                                              verbose)

    # If cell_barcode is mis-named as barcode, convert the column name
    if ("barcode" %in% colnames(identity_amplicon_reads)) {
      colnames(identity_amplicon_reads)[colnames(identity_amplicon_reads) == "barcode"] <- "cell_barcode"
    }

    # Filter identity_barcode_info to only the rows containing barcodes in AmpID
    AmpID_identity_barcode_info <- identity_barcode_info[identity_barcode_info$AmpID == AmpID, ]

    # -- Filter reads according to min_reads threshold --

    invisible(list2env(filter_identity_amplicon_reads(
      identity_amplicon_reads,
      min_reads,
      verbose
    ), envir = environment()))

    # -- Extract identity barcodes using known flanking sequences --

    extracted_barcodes_table <- extract_identity_barcodes(
      identity_amplicon_reads_forward,
      identity_amplicon_reads_reverse,
      AmpID_identity_barcode_info,
      verbose
    )

    extracted_barcodes_tables[[AmpID]] <- extracted_barcodes_table
  }

  if (length(identity_AmpID) > 1 & verbose) {
    cli::cli_h2("Collating Across Identity Amplicons")
  }

  # -- Group identity barcodes --

  cell_to_extracted_BCs <- group_identity_barcodes(extracted_barcodes_tables,
                                                   identity_barcode_info,
                                                   verbose)

  # -- Optional file saving (only if output_dir is provided) --
  if (!is.null(output_dir)) {
    cell_to_extracted_BCS_file_name <- paste0(sample_name, "_cell_to_extracted_BCs.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(cell_to_extracted_BCs,
                     file = file.path(output_dir, cell_to_extracted_BCS_file_name))

    cli::cli_alert_success(
      "Identity barcodes grouped and exported to: {normalizePath(output_dir)}")
  }
}


#' @keywords internal
filter_identity_amplicon_reads <- function(
  identity_amplicon_reads,
  min_reads,
  verbose
) {
  if (verbose) {
    cli::cli_process_start(msg = "Filtering identity amplicon reads...",
                           msg_done = "Identity amplicon reads filtered.")

    # Calculate pre-filter cells and reads total
    pre_filter_cells <- length(unique(identity_amplicon_reads$cell_barcode))
    pre_filter_reads <- length(unique(identity_amplicon_reads$read_name))
  }

  # Mark and filter for barcodes that appear more than min_reads times
  identity_amplicon_reads <- identity_amplicon_reads |>
    dplyr::group_by(cell_barcode) |>
    dplyr::filter(dplyr::n() > min_reads) |>
    dplyr::ungroup()

  # Filter out read names where the forward and reverse read positions are the same
  duplicate_pos <- identity_amplicon_reads$read_name[
    duplicated(identity_amplicon_reads[, c("read_name", "pos")])]
  identity_amplicon_reads <- identity_amplicon_reads[
    !identity_amplicon_reads$read_name %in% duplicate_pos, ]

  # Sort the data frame by read_name, then by pos
  # Such that the first occurrence of each read_name is the forward read
  sorted_reads <- identity_amplicon_reads[
    order(identity_amplicon_reads$read_name, identity_amplicon_reads$pos), ]
  is_second_read <- duplicated(sorted_reads$read_name)

  # First occurrence = forward, second occurrence = reverse
  identity_amplicon_reads_forward <- sorted_reads[!is_second_read, ]
  identity_amplicon_reads_reverse <- sorted_reads[is_second_read, ]

  # Calculate acceptable read pos bounds
  forward_bounds <- quantile(identity_amplicon_reads_forward$pos,
                             probs = c(0.025, 0.975))
  reverse_bounds <- quantile(identity_amplicon_reads_reverse$pos,
                             probs = c(0.025, 0.975))

  # Subset reads with start positions between acceptable bounds
  identity_amplicon_reads_forward <- identity_amplicon_reads_forward |>
    dplyr::filter(dplyr::between(pos, forward_bounds[1], forward_bounds[2]))
  identity_amplicon_reads_reverse <- identity_amplicon_reads_reverse |>
    dplyr::filter(dplyr::between(pos, reverse_bounds[1], reverse_bounds[2]))

  # Remove any now unpaired reads
  identity_amplicon_reads_forward <- dplyr::semi_join(identity_amplicon_reads_forward,
                                                      identity_amplicon_reads_reverse,
                                                      by = "read_name")
  identity_amplicon_reads_reverse <- dplyr::semi_join(identity_amplicon_reads_reverse,
                                                      identity_amplicon_reads_forward,
                                                      by = "read_name")

  if (verbose) {
    cli::cli_process_done()

    # Calculate post-filter cells and reads total
    post_filter_reads <- nrow(identity_amplicon_reads_forward)
    post_filter_cells <- length(unique(identity_amplicon_reads_forward$cell_barcode))

    cli::cli_alert_info("Total reads retained: {post_filter_reads} / {pre_filter_reads}.")
    cli::cli_alert_info("Cells retained: {post_filter_cells} / {pre_filter_cells}.")

    cli::cli_alert_info("Retained forward read start positions: {forward_bounds[1]} - {forward_bounds[2]}.")
    cli::cli_alert_info("Retained reverse read start positions: {reverse_bounds[1]} - {reverse_bounds[2]}.")
  }

  return(list(
    identity_amplicon_reads_forward = identity_amplicon_reads_forward,
    identity_amplicon_reads_reverse = identity_amplicon_reads_reverse
  ))
}


#' @keywords internal
extract_identity_barcodes <- function(
  identity_amplicon_reads_forward,
  identity_amplicon_reads_reverse,
  identity_barcode_info,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Extracting identity barcodes...",
                                      msg_done = "Identity barcodes extracted.")

  reads_by_orientation <- list(forward = identity_amplicon_reads_forward,
                               reverse = identity_amplicon_reads_reverse)

  # Iterate over each barcode in identity_barcode_info
  extracted_barcodes <- purrr::pmap_dfr(
    identity_barcode_info,
    function(barcode_ID, amplicon_seq, start, end, length, optimal_flank_size, orientation, direction, ...) {

      # Extract actual flank sequence, ensuring direction is a valid value
      if (direction == "upstream") {
        flank_seq <- stringr::str_sub(amplicon_seq, start - optimal_flank_size, start - 1)
      } else if (direction == "downstream") {
        flank_seq <- stringr::str_sub(amplicon_seq, end + 1, end + optimal_flank_size)
      } else {
        cli::cli_abort('Each direction in identity_barcode_info must be either "upstream" or "downstream".')
      }

      # Count and locate flank matches in the correctly oriented reads
      reads <- reads_by_orientation[[orientation]]
      match_count <- stringr::str_count(reads$seq, stringr::fixed(flank_seq))
      match_pos   <- stringr::str_locate(reads$seq, stringr::fixed(flank_seq))

      # Only extract from reads where the flank matches exactly once, not 0 or >1
      unambiguous_match <- match_count == 1

      if (direction == "upstream") {
        # Flank is upstream of the barcode, extract coordinates *after* flank
        BC_start <- match_pos[, "end"] + 1
        BC_end   <- match_pos[, "end"] + length
      } else if (direction == "downstream") {
        # Flank is downstream of the barcode, extract coordinates *before* flank
        BC_start <- match_pos[, "start"] - length
        BC_end   <- match_pos[, "start"] - 1
      }

      # Extract barcode sequences of correct length from only single hit reads
      barcode_seq <- rep(NA_character_, nrow(reads))
      barcode_seq[unambiguous_match] <- stringr::str_sub(reads$seq[unambiguous_match],
                                                         BC_start[unambiguous_match],
                                                         BC_end[unambiguous_match])

      output_df <- tibble::tibble(barcode_ID = barcode_ID,
                                  barcode_seq = barcode_seq)
      output_df[["read_name"]] <- reads[["read_name"]]
      output_df
    })

  # Convert to wide table - one row per read pair, one column per barcode
  extracted_barcodes_wide <- extracted_barcodes |>
    tidyr::pivot_wider(id_cols     = dplyr::all_of("read_name"),
                       names_from  = barcode_ID,
                       names_glue  = "{barcode_ID}",
                       values_from = barcode_seq)

  # Keep only read pairs where every barcode was found
  extracted_barcodes_wide <- extracted_barcodes_wide |> tidyr::drop_na()

  # Switch read_name for cell_barcode
  extracted_barcodes_wide <- extracted_barcodes_wide |>
    dplyr::left_join(dplyr::select(identity_amplicon_reads_forward,
                                   read_name, cell_barcode),
                     by = "read_name") |>
    dplyr::relocate(cell_barcode, .after = read_name) |>
    dplyr::select(-read_name)

  if (verbose) {
    cli::cli_process_done()
    n_reads <- nrow(extracted_barcodes_wide)
    cli::cli_alert_info("Complete sets of barcodes extracted for {n_reads} reads.")
  }

  return(extracted_barcodes_wide)
}


#' @keywords internal
group_identity_barcodes <- function(
  extracted_barcodes_tables,
  identity_barcode_info,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Grouping identity barcodes...",
                                      msg_done = "Identity barcodes grouped.")

  # Collapse across extracted_barcodes tables
  collapsed_extracted_barcodes <- purrr::reduce(extracted_barcodes_tables,
                                                dplyr::full_join,
                                                by = "cell_barcode")

  # Group by cell_barcode and each unique combination of BCs, then count occurrences
  cell_to_extracted_BCs <- collapsed_extracted_barcodes |>
    dplyr::group_by_all() |>
    dplyr::count()

  # Filter for only the rows where each BC is the correct number of bases
  cells_to_keep <- rep(TRUE, nrow(cell_to_extracted_BCs))

  for (BC_ID in identity_barcode_info[["barcode_ID"]]) {
    correct_length <- identity_barcode_info[["length"]][identity_barcode_info[["barcode_ID"]] == BC_ID]
    # Flip to FALSE if extracted barcode for BC_ID is the wrong length
    cells_to_keep <- cells_to_keep & (nchar(cell_to_extracted_BCs[[BC_ID]]) == correct_length)
  }

  cell_to_extracted_BCs <- cell_to_extracted_BCs[cells_to_keep, ]

  if (verbose) {
    cli::cli_process_done()
    n_cells <- length(unique(cell_to_extracted_BCs$cell_barcode))
    cli::cli_alert_info("Complete sets of barcodes extracted for {n_cells} cells.")
  }

  return(cell_to_extracted_BCs)
}



