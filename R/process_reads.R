#' Process Reads
#'
#' This function takes various input files (amplicon barcoded data, HDR outcomes,
#' cell to UTR barcode system, sample info, readnames FASTQ, barcodes FASTQ) and
#' performs several key processing steps, inlcuding: quality control, read
#' pairing/match/classification, and data validation.
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param sample_name String containing the file name structure for the
#' experiment sample.
#' @param amplicon_gRNA_number Integer representing the HDR to be processed.
#' @param amplicon_data_file_name String containing the file name structure for
#' the amplicon barcoded data file name.
#' @param HDR_barcoded_file_name String containing the file name structure for
#' the HDR outcomes data file name.
#' @param cell_to_UTR_BCs_file_name String containing the file name structure
#' for the cell to UTR barcode system file name.
#' @param QC_thresholds List containing the following QC parameters:
#' `min_barcode_reads`, `min_read_quality`, `min_start_positions_reads`,
#' `paired_readnames`.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return A data frame containing paired and classified (HDR/WT) read counts
#' and percentages for the provided amplicon_gRNA_number, or HDR.

#' @examples
#' process_reads(
#' source_dir = "data",
#' sample_name = "MB_barcoded",
#' amplicon_gRNA_number = 1,
#' amplicon_data_file_name = "amplicons_barcoded_data.txt",
#' HDR_barcoded_file_name = "HDR_file_barcoded_with_amplicon_sequences.csv",
#' cell_to_UTR_BCs_file_name = "cell_to_UTR_barcode_system.csv",
#' QC_thresholds = list(
#'   min_barcode_reads = 5,
#'   min_read_quality = 30,
#'   min_start_position_reads = 500,
#'   paired_readnames = TRUE
#'   ),
#' output_dir = "data",
#' verbose
#' )

#' @importFrom dplyr %>%

#' @export
process_reads <- function(
  source_dir = "",
  sample_name = "sample",
  amplicon_gRNA_number,
  amplicon_data_file_name,
  HDR_barcoded_file_name,
  cell_to_UTR_BCs_file_name,
  QC_thresholds,
  output_dir = NULL,
  verbose
) {
  cli::cli_alert_info(
    "-- Processing amplicon_{amplicon_gRNA_number} from sample: {sample_name} --"
  )

  # -- Import Data Files --

  # Import relevant data files, and load to the function environment
  imported_data <- import_process_reads_data(
    source_dir,
    sample_name,
    amplicon_gRNA_number,
    amplicon_data_file_name,
    HDR_barcoded_file_name,
    cell_to_UTR_BCs_file_name,
    verbose
  )

  invisible(list2env(imported_data, envir = environment()))

  cli::cli_alert_success("Successfully imported datasets.")

  # -- Filter for Selected Amplicon --

  mapped_reads_df_amplicon <- filter_for_selected_amplicon(
    mapped_reads_df,
    HDR_barcoded_outcomes,
    amplicon_gRNA_number,
    df_BC_system,
    verbose
  )

  # -- Mapped Reads Quality Control --

  mapped_reads_df_amplicon <- process_reads_QC(
    mapped_reads_df_amplicon,
    QC_thresholds,
    verbose
  )

  # -- Pair Reads --

  mapped_reads_df_amplicon_paired_count <- pair_reads(
    mapped_reads_df_amplicon,
    verbose
  )

  # -- Merge R1 and R2
  mapped_reads_df_amplicon_paired_count <-
    cbind(mapped_reads_df_amplicon_paired_count,
          merge_r1_r2(mapped_reads_df_amplicon_paired_count$sequence_R1,
                      mapped_reads_df_amplicon_paired_count$sequence_R2,
                      min_overlap = 2,
                      verbose))

  # -- Classify Reads --

  mapped_reads_df_amplicon_paired_count <- classify_processed_reads(
    mapped_reads_df_amplicon_paired_count,
    HDR_barcoded_outcomes,
    amplicon_gRNA_number,
    verbose
  )

  # -- Validate and Save Results --

  df_to_save_count <- export_process_reads(
    mapped_reads_df_amplicon_paired_count,
    amplicon_gRNA_number,
    sample_name,
    output_dir
  )

  return(df_to_save_count)
}


#' @keywords internal
import_process_reads_data <- function(
  source_dir = "",
  sample_name = "sample",
  amplicon_gRNA_number = "n",
  amplicon_data_file_name,
  HDR_barcoded_file_name,
  cell_to_UTR_BCs_file_name,
  verbose
) {
  # -- Read in the amplicon data --

  amplicon_data <- import_amplicons_barcoded_data(
    source_dir,
    amplicon_data_file_name,
    verbose
  )

  amplicon_data <- amplicon_data %>%
    dplyr::mutate(amplicon_gRNA = paste0(AMPLICON_NAME, "_", sgRNA_SEQUENCE)) %>%
    dplyr::filter(!is.na(sgRNA_SEQUENCE))

  # -- Read in the HDR outcome file --

  HDR_barcoded_outcomes <- import_HDR_barcoded_data(
    source_dir,
    HDR_barcoded_file_name,
    verbose
  )

  # -- Load mapped reads --

  mapped_reads_file_name <- paste0(sample_name, "_info.csv")
  mapped_reads_df <- import_sample_info_data(
    source_dir,
    mapped_reads_file_name,
    verbose
  )

  # -- Load barcodes --

  barcodes_file_name <- paste0(sample_name, "_barcodes_fastq.txt")
  barcodes <- import_barcode_reads_data(
    source_dir,
    barcodes_file_name,
    verbose
  )

  # -- Load readnames --

  readnames_file_name <- paste0(sample_name, "_R2_001_readnames_fastq.txt")
  readnames <- import_readnames_data(
    source_dir,
    readnames_file_name,
    verbose
  )

  # Map barcodes to readnames
  names(barcodes) <- readnames

  # Map mapped_reads_df to barcodes
  mapped_reads_df$barcode <- barcodes[mapped_reads_df$readname]

  # -- Import list of barcodes from UTR_BC counts --

  df_BC_system <- import_cell_to_UTR_barcode_system(
    source_dir,
    cell_to_UTR_BCs_file_name,
    verbose
  )
  df_BC_system <- df_BC_system %>%
    dplyr::rename(barcode = cell_barcode)

  # Return as list to be unpacked into the global function environment
  return(list(
    amplicon_data = amplicon_data,
    HDR_barcoded_outcomes = HDR_barcoded_outcomes,
    mapped_reads_df = mapped_reads_df,
    barcodes = barcodes,
    readnames = readnames,
    df_BC_system = df_BC_system
  ))
}


#' @keywords internal
filter_for_selected_amplicon <- function(
  mapped_reads_df,
  HDR_barcoded_outcomes,
  amplicon_gRNA_number,
  df_BC_system,
  verbose
) {
  if (verbose) cli::cli_inform("Filtering for Amplicon {amplicon_gRNA_number}...")

  # Filter for selected amplicon
  mapped_reads_df_amplicon <- mapped_reads_df[
    mapped_reads_df$amplicon_name == HDR_barcoded_outcomes$AMPLICON_NAME[amplicon_gRNA_number], ]

  # Only keep the barcodes listed form UTR_BC counts
  # (cell_to_UTR_barcode_system.csv)
  mapped_reads_df_amplicon <- mapped_reads_df_amplicon[
    mapped_reads_df_amplicon$barcode %in% df_BC_system$barcode, ]

  if (verbose) cli::cli_alert_success(
    "Successfully filtered for Amplicon {amplicon_gRNA_number}.")

  return(mapped_reads_df_amplicon)
}


#' @keywords internal
process_reads_QC <- function(
  mapped_reads_df_amplicon,
  QC_thresholds,
  verbose
) {
  if (verbose) cli::cli_inform("Performing quality control...")

  # Keep barcodes with reads >= QC_thresholds$min_barcode_reads
  barcode_counts <- table(mapped_reads_df_amplicon$barcode)
  barcodes_keep <-
    names(barcode_counts)[barcode_counts >= QC_thresholds$min_barcode_reads]
  mapped_reads_df_amplicon <-
    mapped_reads_df_amplicon[mapped_reads_df_amplicon$barcode %in% barcodes_keep, ]

  # Keep reads with quality >= QC_thresholds$min_read_quality
  mapped_reads_df_amplicon <-
    mapped_reads_df_amplicon[mapped_reads_df_amplicon$quality >= QC_thresholds$min_read_quality, ]

  # Keep start positions with >= min_start_position_reads reads
  start_positions_counts <- table(mapped_reads_df_amplicon$start)
  start_positions_keep <-
    names(start_positions_counts)[start_positions_counts > QC_thresholds$min_start_position_reads]
  mapped_reads_df_amplicon <-
    mapped_reads_df_amplicon[mapped_reads_df_amplicon$start %in% start_positions_keep, ]

  # Keep only readnames that appear exactly twice (paired)
  if(QC_thresholds$paired_readnames == TRUE) {
    readname_counts <- table(mapped_reads_df_amplicon$readname)
    reads_keep <- names(readname_counts)[readname_counts == 2]
    mapped_reads_df_amplicon <-
      mapped_reads_df_amplicon[mapped_reads_df_amplicon$readname %in% reads_keep, ]
  }

  # Order by readname
  mapped_reads_df_amplicon <-
    mapped_reads_df_amplicon[order(mapped_reads_df_amplicon$readname), ]

  if (verbose) cli::cli_alert_success("Successfully performed quality control.")

  return(mapped_reads_df_amplicon)
}


#' @keywords internal
pair_reads <- function(
  mapped_reads_df_amplicon,
  verbose
) {
  if (verbose) cli::cli_inform("Pairing reads...")

  # Pair reads (R1 and R2)
  paired_idx <- seq(1, nrow(mapped_reads_df_amplicon) - 1, by = 2)
  mapped_reads_df_amplicon_paired <- cbind(
    mapped_reads_df_amplicon[paired_idx,     c("readname", "amplicon_name",
                                               "barcode",  "start",
                                               "cigar",    "sequence")],
    mapped_reads_df_amplicon[paired_idx + 1, c("readname", "amplicon_name",
                                               "barcode", "start",
                                               "cigar",   "sequence")]
  )
  colnames(mapped_reads_df_amplicon_paired) <- c(
    "readname_R1", "amplicon_name",  "barcode",
    "start_R1",    "cigar_R1",       "sequence_R1",
    "readname_R2", "amplicon_name2", "barcode2",
    "start_R2",    "cigar_R2",       "sequence_R2"
  )

  # Keep only pairs where readnames match
  mapped_reads_df_amplicon_paired <-
    mapped_reads_df_amplicon_paired[mapped_reads_df_amplicon_paired$readname_R1 ==
                                      mapped_reads_df_amplicon_paired$readname_R2, ]

  # Remove duplicate columns
  mapped_reads_df_amplicon_paired <-
    mapped_reads_df_amplicon_paired[, c("amplicon_name", "barcode",  "start_R1",
                                        "cigar_R1",      "start_R2", "cigar_R2",
                                        "sequence_R1",   "sequence_R2")]

  # Ensure R1 is always the lower start
  swap_idx <- mapped_reads_df_amplicon_paired$start_R1 >
                mapped_reads_df_amplicon_paired$start_R2

  if (any(swap_idx)) {
    temp <- mapped_reads_df_amplicon_paired[swap_idx, ]
    temp <- temp %>%
      dplyr::mutate(
        start_R1 = start_R2,
        start_R2 = start_R1,
        sequence_R1 = sequence_R2,
        sequence_R2 = sequence_R1,
        cigar_R1 = cigar_R2,
        cigar_R2 = cigar_R1
      )
    mapped_reads_df_amplicon_paired[swap_idx, ] <- temp
  }

  # Count pairs
  mapped_reads_df_amplicon_paired_count <- mapped_reads_df_amplicon_paired %>%
    dplyr::group_by_all() %>%
    dplyr::tally(name = "n") %>%
    dplyr::ungroup()

  # Convert to data frame and apply rownames
  mapped_reads_df_amplicon_paired_count <-
    as.data.frame(mapped_reads_df_amplicon_paired_count)
  rownames(mapped_reads_df_amplicon_paired_count) <-
    seq_len(nrow(mapped_reads_df_amplicon_paired_count))

  if (verbose) cli::cli_alert_success("Successfully paired reads.")

  return(mapped_reads_df_amplicon_paired_count)
}


#' @keywords internal
classify_processed_reads <- function(
  mapped_reads_df_amplicon_paired_count,
  HDR_barcoded_outcomes,
  amplicon_gRNA_number,
  verbose
) {
  if (verbose) cli::cli_inform("Classifying processed reads...")

  # Classify WT
  mapped_reads_df_amplicon_paired_count$is_WT <-
    check_sequences_contain_string(
      mapped_reads_df_amplicon_paired_count,
      HDR_barcoded_outcomes$amplicon_seq_in_window[amplicon_gRNA_number],
      column_name = "merged_sequence",
      batch_size = 10000,
      return_indices = FALSE,
      strict_dna = TRUE
    )

  mapped_reads_df_amplicon_paired_count$HDR <- "NA"
  mapped_reads_df_amplicon_paired_count$HDR[mapped_reads_df_amplicon_paired_count$is_WT] <- "WT"

  # Classify HDR
  mapped_reads_df_amplicon_paired_count$is_HDR <-
    check_sequences_contain_string(
      mapped_reads_df_amplicon_paired_count,
      HDR_barcoded_outcomes$HDR_seq_in_window[amplicon_gRNA_number],
      column_name = "merged_sequence",
      batch_size = 10000,
      return_indices = FALSE,
      strict_dna = TRUE
    )

  HDR_name <- paste0("HDR_", amplicon_gRNA_number)
  mapped_reads_df_amplicon_paired_count$HDR[mapped_reads_df_amplicon_paired_count$is_HDR] <- HDR_name

  if (verbose) cli::cli_alert_success("Successfully Classified processed reads.")

  return(mapped_reads_df_amplicon_paired_count)
}


#' @keywords internal
export_process_reads <- function(
  mapped_reads_df_amplicon_paired_count,
  amplicon_gRNA_number,
  sample_name,
  output_dir
) {
  if (verbose) cli::cli_inform("Exporting results...")

  # Modify results table to select columns to export
  df_to_save <- mapped_reads_df_amplicon_paired_count[,-1]
  df_to_save <- df_to_save %>%
    dplyr::filter(!(is.na(merged_sequence))) %>%
    dplyr::select(barcode, HDR)

  # Calculate n (total) and pct (percentage that match classification)
  df_to_save_count <- df_to_save %>%
    dplyr::group_by(barcode, HDR) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup()

  # Rename columns for classified read summary stats
  renamed_col <- paste0("HDR_", amplicon_gRNA_number)
  names(df_to_save_count)[2] <- renamed_col

  renamed_col <- paste0("n_", amplicon_gRNA_number)
  names(df_to_save_count)[3] <- renamed_col

  renamed_col <- paste0("pct_", amplicon_gRNA_number)
  names(df_to_save_count)[4] <- renamed_col

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    file_save_name <- file.path(
      output_dir,
      paste0("reads_classified_", sample_name, "_HDR_", amplicon_gRNA_number, ".csv")
    )

    readr::write_csv(df_to_save_count, file_save_name)

    cli::cli_alert_success("Results saved to: {normalizePath(output_dir)}.")
  }

  return(df_to_save_count)
}

