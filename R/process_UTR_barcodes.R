#' Process UTR Barcodes
#'
#' This function takes barcode reads data and first performs some basic QC
#' filtering. It then extracts BC1 and BC2/BC3s from provided flanking
#' sequences, before grouping cells by cell_barcode and counting reads.
#'
#' @inheritParams global_params source_dir sample_name UTR_barcode_AmpID output_dir verbose
#'
#' @param min_reads Quality control threshold: minimum number of reads.
#' @param seq_before_barcode_1 Barcode 1 flanking sequence.
#' @param seq_before_barcode_2 Barcode 2 flanking sequence.
#' @param seq_after_barcode_3 Barcode 3 flanking sequence.
#'
#' @return A table linking cells to extracted barcodes.
#'
#' @importFrom dplyr %>%
#'
#' @export
process_UTR_barcodes <- function(
  source_dir,
  sample_name,
  UTR_barcode_AmpID,
  min_reads,
  seq_before_barcode_1,
  seq_before_barcode_2,
  seq_after_barcode_3,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import datasets --

  UTR_barcode_reads_file_name <- paste0(sample_name, "_",
                                        UTR_barcode_AmpID, "_data.tsv.gz")

  UTR_barcode_reads <- import_UTR_barcode_reads(
    source_dir,
    UTR_barcode_reads_file_name,
    verbose
  )

  print(colnames(UTR_barcode_reads))

  if ("barcode" %in% colnames(UTR_barcode_reads)) {
    colnames(UTR_barcode_reads)[colnames(UTR_barcode_reads) == "barcode"] <- "cell_barcode"
  }

  print(colnames(UTR_barcode_reads))

  # -- Filter reads according at min_reads threshold --

  invisible(list2env(filter_UTR_barcode_reads(
    UTR_barcode_reads,
    min_reads,
    verbose
  ), envir = environment()))

  # -- Extract UTR barcodes using known flanking sequences --

  readname_barcodes_table <- extract_UTR_barcodes(
    filtered_UTR_barcode_reads_forward,
    filtered_UTR_barcode_reads_reverse,
    seq_before_barcode_1,
    seq_before_barcode_2,
    seq_after_barcode_3,
    verbose
  )

  # -- Group UTR barcodes --

  cell_to_extracted_BCs <- group_UTR_barcodes(readname_barcodes_table, verbose)

  # -- Optional file saving (only if output_dir is provided) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(cell_to_extracted_BCs,
            file = file.path(output_dir, "cell_to_extracted_BCs.rds"))

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
  filtered_UTR_barcode_reads_reverse <- filtered_UTR_barcode_reads_reverse %>%
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
  seq_before_barcode_1,
  seq_before_barcode_2,
  seq_after_barcode_3,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Extracting UTR barcodes...",
                                      msg_done = "UTR barcodes extracted.")

  # Filter FORWARD reads that contain both flanking sequences BEFORE the barcode
  UTR_barcode_reads_with_seq_before_BC_1_2 <-
    UTR_barcode_reads_forward[grepl_seq1_2(
      UTR_barcode_reads_forward$seq,
      seq_before_barcode_1, seq_before_barcode_2), ]

  # Filter REVERSE reads that contain the flanking sequence AFTER the barcode
  UTR_barcode_reads_with_seq_after_barcode_3 <-
    UTR_barcode_reads_reverse[sapply(
      UTR_barcode_reads_reverse$seq,
      function(x) grepl(seq_after_barcode_3, x)), ]

  # Apply extract_BC1 to each sequence and store results in BC1 column
  UTR_barcode_reads_with_seq_before_BC_1_2$BC1 <-
    unlist(sapply(UTR_barcode_reads_with_seq_before_BC_1_2$seq,
                  function(x) extract_BC1(x, seq_before_barcode_1)))

  # Apply extract_BC2 to each sequence and store results in BC2 column
  UTR_barcode_reads_with_seq_before_BC_1_2$BC2 <-
    unlist(sapply(UTR_barcode_reads_with_seq_before_BC_1_2$seq,
                  function(x) extract_BC2(x, seq_before_barcode_2)))

  # Create a table with selected columns: read_name, cell_barcode, BC1, and BC2
  readname_barcodes_table <-
    UTR_barcode_reads_with_seq_before_BC_1_2 %>%
    dplyr::select(read_name, cell_barcode, BC1, BC2)

  # Apply extract_BC3 to each sequence and store results in BC3 column
  UTR_barcode_reads_with_seq_after_barcode_3$BC3 <-
    sapply(UTR_barcode_reads_with_seq_after_barcode_3$seq,
           function(x) extract_BC3(x, seq_after_barcode_3))

  # Remove duplicate readnames to ensure unique mapping
  UTR_barcode_reads_with_seq_after_barcode_3 <-
    UTR_barcode_reads_with_seq_after_barcode_3 %>%
    dplyr::filter(!(duplicated(read_name)))

  readname_barcodes_table <- readname_barcodes_table %>%
    dplyr::left_join(
      UTR_barcode_reads_with_seq_after_barcode_3 %>%
        dplyr::select(read_name, BC3),
      by = "read_name"
    )

  # Remove rows with any missing (NA) values from the table
  readname_barcodes_table <-
    readname_barcodes_table[stats::complete.cases(readname_barcodes_table), ]

  if (verbose) cli::cli_process_done()

  return(readname_barcodes_table)
}


#' @keywords internal
group_UTR_barcodes <- function(
  readname_barcodes_table,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Grouping UTR barcodes...",
                                      msg_done = "UTR barcodes grouped.")

  # Group by cell_barcode, BC1, BC2, and BC3, then count occurrences
  cell_to_extracted_BCs <- readname_barcodes_table %>%
    dplyr::select(cell_barcode, BC1, BC2, BC3) %>%
    dplyr::group_by_all()%>%
    dplyr::count()

  # Remove rows with any missing values
  cell_to_extracted_BCs <- cell_to_extracted_BCs[stats::complete.cases(cell_to_extracted_BCs),]

  # Filter for only the rows where BC1 is 20 bases, and BC2 and BC3 are 30 bases
  cell_to_extracted_BCs <- cell_to_extracted_BCs %>%
    dplyr::filter(nchar(BC1) == 20, nchar(BC2) == 30, nchar(BC3) == 30)

  if (verbose) cli::cli_process_done()

  return(cell_to_extracted_BCs)
}



