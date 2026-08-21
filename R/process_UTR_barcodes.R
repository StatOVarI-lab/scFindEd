#' Process UTR Barcodes
#'
#' This function takes barcode reads data and first performs some basic QC
#' filtering. It then extracts BC1 and BC2/BC3s from provided flanking
#' sequences, before grouping cells by cell_barcode and counting reads.
#'
#' @inheritParams global_params source_dir sample_name UTR_barcode_AmpID min_reads output_dir verbose
#'
#' @return A table linking cells to extracted barcodes.
#'
#' @export
process_UTR_barcodes <- function(
  source_dir,
  sample_name,
  UTR_barcode_AmpID,
  min_reads,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import datasets --

  UTR_barcode_reads_file_name <- paste0(sample_name, "_",
                                        UTR_barcode_AmpID, "_data.tsv.gz")
  UTR_barcode_reads <- import_UTR_barcode_reads(source_dir,
                                                UTR_barcode_reads_file_name,
                                                verbose)

  UTR_barcode_info_file_name <- paste0(sample_name, "_UTR_barcode_info.csv")
  UTR_barcode_info <- import_UTR_barcode_info(source_dir,
                                              UTR_barcode_info_file_name,
                                              verbose)

  # If cell_barcode is mis-named as barcode, convert the column name
  if ("barcode" %in% colnames(UTR_barcode_reads)) {
    colnames(UTR_barcode_reads)[colnames(UTR_barcode_reads) == "barcode"] <- "cell_barcode"
  }

  # -- Filter reads according at min_reads threshold --

  invisible(list2env(filter_UTR_barcode_reads(
    UTR_barcode_reads,
    min_reads,
    verbose
  ), envir = environment()))

  # -- Extract UTR barcodes using known flanking sequences --

  extracted_barcodes_table <- extract_UTR_barcodes(
    filtered_UTR_barcode_reads_forward,
    filtered_UTR_barcode_reads_reverse,
    UTR_barcode_info,
    verbose
  )

  # -- Group UTR barcodes --

  cell_to_extracted_BCs <- group_UTR_barcodes(extracted_barcodes_table, verbose)

  # -- Optional file saving (only if output_dir is provided) --
  if (!is.null(output_dir)) {
    cell_to_extracted_BCS_file_name <- paste0(sample_name, "_cell_to_extracted_BCs.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(cell_to_extracted_BCs,
                     file = file.path(output_dir, cell_to_extracted_BCS_file_name))

    cli::cli_alert_success(
      "Grouped UTR barcodes exported to: {normalizePath(output_dir)}")
  }
}


#' @keywords internal
filter_UTR_barcode_reads <- function(
  UTR_barcode_reads,
  min_reads,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Filtering UTR barcode reads...",
                                      msg_done = "UTR barcode reads filtered.")

  # Mark and filter for barcodes that appear more than min_reads times
  barcodes.keep <-
    names(table(UTR_barcode_reads$cell_barcode))[table(UTR_barcode_reads$cell_barcode)>min_reads]
  filtered_UTR_barcode_reads <-
    UTR_barcode_reads[UTR_barcode_reads$cell_barcode%in%barcodes.keep,]

  # Subset of forward reads with a start position between 40 and 43
  filtered_UTR_barcode_reads_forward <-
    filtered_UTR_barcode_reads[filtered_UTR_barcode_reads$pos%in%40:43,]
  # Subset of reverse reads with a start position between 142 and 148
  filtered_UTR_barcode_reads_reverse <-
    filtered_UTR_barcode_reads[filtered_UTR_barcode_reads$pos%in%142:148,]

  # Removing the two reads with unusual start position for forward read
  # and are included for both forward and reverse read
  filtered_UTR_barcode_reads_reverse <- filtered_UTR_barcode_reads_reverse |>
    dplyr::filter(!duplicated(read_name))

  if (verbose) cli::cli_process_done()

  return(list(
    filtered_UTR_barcode_reads_forward = filtered_UTR_barcode_reads_forward,
    filtered_UTR_barcode_reads_reverse = filtered_UTR_barcode_reads_reverse
  ))
}


#' @keywords internal
extract_UTR_barcodes <- function(
  UTR_barcode_reads_forward,
  UTR_barcode_reads_reverse,
  UTR_barcode_info,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Extracting UTR barcodes...",
                                      msg_done = "UTR barcodes extracted.")

  reads_by_orientation <- list(forward = UTR_barcode_reads_forward,
                               reverse = UTR_barcode_reads_reverse)

  # Iterate over each barcode in UTR_barcode_info
  extracted_barcodes <- purrr::pmap_dfr(
    UTR_barcode_info,
    function(barcode_ID, amplicon_seq, start, end, length, optimal_flank_size, orientation, direction, ...) {

      # Extract actual flank sequence, ensuring direction is a valid value
      if (direction == "upstream") {
        flank_seq <- stringr::str_sub(amplicon_seq, start - optimal_flank_size, start - 1)
      } else if (direction == "downstream") {
        flank_seq <- stringr::str_sub(amplicon_seq, end + 1, end + optimal_flank_size)
      } else {
        cli::cli_abort('Each direction in UTR_barcode_info must be either "upstream" or "downstream".')
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
                       names_glue  = "BC{barcode_ID}",
                       values_from = barcode_seq)

  # Keep only read pairs where every barcode was found
  n_reads <- nrow(extracted_barcodes_wide)
  extracted_barcodes_wide <- extracted_barcodes_wide |> tidyr::drop_na()

  # Switch read_name for cell_barcode
  extracted_barcodes_wide <- extracted_barcodes_wide |>
    dplyr::left_join(dplyr::select(UTR_barcode_reads_forward,
                                   read_name, cell_barcode),
                     by = "read_name") |>
    dplyr::relocate(cell_barcode, .after = read_name) |>
    dplyr::select(-read_name)

  if (verbose) {
    cli::cli_process_done()
    n_reads_all_barcodes <- nrow(extracted_barcodes_wide)
    cli::cli_alert_info("Complete set of barcodes extracted for {n_reads_all_barcodes} out of {n_reads} total reads.")
  }

  return(extracted_barcodes_wide)
}


#' @keywords internal
group_UTR_barcodes <- function(
  extracted_barcodes_table,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Grouping UTR barcodes...",
                                      msg_done = "UTR barcodes grouped.")

  # Group by cell_barcode, BC1, BC2, and BC3, then count occurrences
  cell_to_extracted_BCs <- extracted_barcodes_table |>
    dplyr::group_by_all() |>
    dplyr::count()

  # Filter for only the rows where BC1 is 20 bases, and BC2 and BC3 are 30 bases
  cell_to_extracted_BCs <- cell_to_extracted_BCs |>
    dplyr::filter(nchar(BC1) == 20, nchar(BC2) == 30, nchar(BC3) == 30)

  if (verbose) {
    cli::cli_process_done()
    n_cells <- nrow(cell_to_extracted_BCs)
    cli::cli_alert_info("Complete set of barcodes extracted for {n_cells} cells.")
  }

  return(cell_to_extracted_BCs)
}



