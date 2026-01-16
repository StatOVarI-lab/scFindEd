#' Process Barcodes
#'
#' This function takes barcode reads data and first performs some basic QC
#' filtering. It then extracts BC1 and BC2/BC3s from provided flanking
#' sequences, before grouping cells by cell_barcode and counting reads.
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param barcode_amplicon_file_name String containing the file name structure
#' for the barcode amplicon information file.
#' @param cell_barcodes_file_name String containing the file name structure for
#' the extracted read barcodes.
#' @param readnames_file_name String containing the file name structure for the
#' extracted readnames.
#' @param min_reads QC threshold: minimum number of reads.
#' @param quality_score_threshold QC threshold: minimum quality score.
#' @param seq_before_barcode_1 BC1 flanking sequence.
#' @param seq_before_barcode_2 BC2 flanking sequence.
#' @param seq_after_barcode_3 BC3 flanking sequence.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return A table linking cells to extracted barcodes.

#' @examples
#' process_barcodes(
#' source_dir = "data",
#' barcode_amplicon_file_name = "barcode_amplicon_info.txt",
#' cell_barcodes_file_name = "MB_barcoded_barcodes_fastq.txt",
#' readnames_file_name = "MB_barcoded_R2_001_readnames_fastq.txt",
#' min_reads = 2,
#' quality_score_threshold = 30,
#' seq_before_barcode_1 = "GTACGAAGACTG",
#' seq_before_barcode_2 = "CGATACTTGTTCGATCCTTCGC",
#' seq_after_barcode_3 = "CTTGTCTCGCTCCGGG",
#' output_dir = "data",
#' verbose = TRUE)

#' @importFrom dplyr %>%

#' @export
process_barcodes <- function(
  source_dir,
  barcode_amplicon_file_name,
  cell_barcodes_file_name,
  readnames_file_name,
  min_reads,
  quality_score_threshold,
  seq_before_barcode_1,
  seq_before_barcode_2,
  seq_after_barcode_3,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import datasets --

  imported_data <- import_process_barcodes_data(
    source_dir,
    barcode_amplicon_file_name,
    cell_barcodes_file_name,
    readnames_file_name,
    verbose
  )

  invisible(list2env(imported_data, envir = environment()))

  # -- Filter amplicon data according to QC thresholds --

  filtered_amplicon_data <- filter_amplicon_data(
    barcode_amplicon_data,
    min_reads,
    quality_score_threshold,
    verbose
  )

  invisible(list2env(filtered_amplicon_data, envir = environment()))

  # -- Extract barcodes using known flanking sequences --

  readname_cell_table <- extract_barcodes(
    barcode_amplicon_data_forward,
    barcode_amplicon_data_reverse,
    seq_before_barcode_1,
    seq_before_barcode_2,
    seq_after_barcode_3,
    verbose
  )

  cell_BC1_BC2_BC3 <- group_barcodes(readname_cell_table, verbose)

  # -- Optional file saving (only if output_dir is provided) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(cell_BC1_BC2_BC3,
            file = file.path(output_dir, "cell_to_extracted_BCs.rds"))

    cli::cli_alert_success(
      "Grouped barcodes exported to: {normalizePath(output_dir)}")
  }

  # Return the grouped data invisibly
  invisible(cell_BC1_BC2_BC3)
}


#' @keywords internal
import_process_barcodes_data <- function(
  source_dir,
  barcode_amplicon_file_name,
  cell_barcodes_file_name,
  readnames_file_name,
  verbose
) {
  # -- Import barcode amplicon information --

  barcode_amplicon_data <- import_barcode_amplicon_info(
    source_dir,
    barcode_amplicon_file_name,
    verbose
  )

  # -- Import extracted read barcodes --

  cell_barcodes <- import_barcode_reads_data(
    source_dir,
    cell_barcodes_file_name,
    verbose
  )

  # -- Import extracted readnames --

  readnames <- import_readnames_data(
    source_dir,
    readnames_file_name,
    verbose
  )

  # Assign readnames as the names of the cell_barcodes vector
  names(cell_barcodes) <- readnames

  # Map cell barcodes to the amplicon data, using readname as the key
  barcode_amplicon_data$cell_barcode <-
    cell_barcodes[barcode_amplicon_data$readname]

  cli::cli_alert_success("Successfully imported datasets.")

  return(list(
    barcode_amplicon_data = barcode_amplicon_data,
    cell_barcodes = cell_barcodes,
    readnames = readnames
  ))
}


#' @keywords internal
filter_amplicon_data <- function(
  barcode_amplicon_data,
  min_reads,
  quality_score_threshold,
  verbose
) {
  # Mark and filter for barcodes that appear more than @min_reads times
  barcodes.keep <-
    names(table(barcode_amplicon_data$cell_barcode))[table(barcode_amplicon_data$cell_barcode)>min_reads]
  barcode_amplicon_data <-
    barcode_amplicon_data[barcode_amplicon_data$cell_barcode%in%barcodes.keep,]

  # Filter the reads with a quality score < @quality_score_threshold
  barcode_amplicon_data <- barcode_amplicon_data[barcode_amplicon_data$quality>=quality_score_threshold,]

  # Subset of forward reads with a start position between 40 and 43
  barcode_amplicon_data_forward <- barcode_amplicon_data[barcode_amplicon_data$start%in%40:43,]
  # Subset of reverse reads with a start position between 142 and 148
  barcode_amplicon_data_reverse <- barcode_amplicon_data[barcode_amplicon_data$start%in%142:148,]

  # Removing the two reads with unusual start position for forward read
  # and are included for both forward and reverse read
  barcode_amplicon_data_reverse <- barcode_amplicon_data_reverse %>%
    dplyr::filter(!duplicated(readname))

  if (verbose) cli::cli_alert_success("Amplicon data filtered.")

  return(list(
    barcode_amplicon_data_forward = barcode_amplicon_data_forward,
    barcode_amplicon_data_reverse = barcode_amplicon_data_reverse
  ))
}


#' @keywords internal
extract_barcodes <- function(
  barcode_amplicon_data_forward,
  barcode_amplicon_data_reverse,
  seq_before_barcode_1,
  seq_before_barcode_2,
  seq_after_barcode_3,
  verbose
) {
  # Filter FORWARD reads that contain both flanking sequences BEFORE the barcode
  barcode_amplicon_data_with_seq_before_BC_1_2 <-
    barcode_amplicon_data_forward[grepl_seq1_2(
      barcode_amplicon_data_forward$sequence,
      seq_before_barcode_1, seq_before_barcode_2), ]

  # Filter REVERSE reads that contain the flanking sequence AFTER the barcode
  barcode_amplicon_data_with_seq_after_barcode_3 <-
    barcode_amplicon_data_reverse[sapply(
      barcode_amplicon_data_reverse$sequence,
      function(x) grepl(seq_after_barcode_3, x)), ]

  # Apply extract_BC1 to each sequence and store results in BC1 column
  barcode_amplicon_data_with_seq_before_BC_1_2$BC1 <-
    unlist(sapply(barcode_amplicon_data_with_seq_before_BC_1_2$sequence,
                  function(x) extract_BC1(x, seq_before_barcode_1)))

  # Apply extract_BC2 to each sequence and store results in BC2 column
  barcode_amplicon_data_with_seq_before_BC_1_2$BC2 <-
    unlist(sapply(barcode_amplicon_data_with_seq_before_BC_1_2$sequence,
                  function(x) extract_BC2(x, seq_before_barcode_2)))

  # Create a table with selected columns: readname, cell_barcode, BC1, and BC2
  readname_cell_table <-
    barcode_amplicon_data_with_seq_before_BC_1_2 %>%
    dplyr::select(readname, cell_barcode, BC1, BC2)

  # Apply extract_BC3 to each sequence and store results in BC3 column
  barcode_amplicon_data_with_seq_after_barcode_3$BC3 <-
    sapply(barcode_amplicon_data_with_seq_after_barcode_3$sequence,
           function(x) extract_BC3(x, seq_after_barcode_3))

  # Remove duplicate readnames to ensure unique mapping
  barcode_amplicon_data_with_seq_after_barcode_3 <-
    barcode_amplicon_data_with_seq_after_barcode_3 %>%
    dplyr::filter(!(duplicated(readname)))

  readname_cell_table <- readname_cell_table %>%
    dplyr::left_join(
      barcode_amplicon_data_with_seq_after_barcode_3 %>%
        dplyr::select(readname, BC3),
      by = "readname"
    )

  # Remove rows with any missing (NA) values from the table
  readname_cell_table <-
    readname_cell_table[stats::complete.cases(readname_cell_table), ]

  if (verbose) cli::cli_alert_success("Barcodes extracted.")

  return(readname_cell_table)
}


#' @keywords internal
group_barcodes <- function(
  readname_cell_table,
  verbose
) {
  # Group by cell_barcode, BC1, BC2, and BC3, then count occurrences
  cell_BC1_BC2_BC3 <- readname_cell_table %>%
    dplyr::select(cell_barcode,BC1,BC2,BC3) %>%
    dplyr::group_by_all()%>%
    dplyr::count()

  # Remove rows with any missing values
  cell_BC1_BC2_BC3 <- cell_BC1_BC2_BC3[stats::complete.cases(cell_BC1_BC2_BC3),]

  # Filter for only the rows where BC1 is 20 bases, and BC2 and BC3 are 30 bases
  cell_BC1_BC2_BC3 <- cell_BC1_BC2_BC3 %>%
    dplyr::filter(nchar(BC1) == 20, nchar(BC2) == 30, nchar(BC3) == 30)

  if (verbose) cli::cli_alert_success("Barcodes grouped.")

  return(cell_BC1_BC2_BC3)
}



