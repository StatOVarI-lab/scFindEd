#' Process Barcodes
#'
#' This function takes barcode reads data and first performs some basic QC
#' filtering. It then extracts BC1 and BC2/BC3s from provided flanking
#' sequences, before grouping cells by cell_barcode and counting reads.
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param cell_AmpID String containing the AmpID for the cell-barcode amplicon.
#' @param min_reads QC threshold: minimum number of reads.
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
#' cell_AmpID = "AMPL902654",
#' min_reads = 2,
#' seq_before_barcode_1 = "GTACGAAGACTG",
#' seq_before_barcode_2 = "CGATACTTGTTCGATCCTTCGC",
#' seq_after_barcode_3 = "CTTGTCTCGCTCCGGG",
#' output_dir = "data",
#' verbose = TRUE)

#' @importFrom dplyr %>%

#' @export
process_barcodes <- function(
  source_dir,
  cell_AmpID,
  min_reads,
  seq_before_barcode_1,
  seq_before_barcode_2,
  seq_after_barcode_3,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import datasets --

  mapped_cell_barcode_reads <- import_barcoded_mapped_reads(
    source_dir, 
    cell_AmpID, 
    verbose
  )
  if ("barcode" %in% colnames(mapped_cell_barcode_reads)) {
    colnames(mapped_cell_barcode_reads)[colnames(mapped_cell_barcode_reads) == "barcode"] <- "cell_barcode"
  }

  # -- Filter cell-barcode reads according at min_reads threshold --

  filtered_cell_barcode_reads <- filter_cell_barcode_reads(
    mapped_cell_barcode_reads,
    min_reads,
    verbose
  )

  invisible(list2env(filtered_cell_barcode_reads, envir = environment()))

  # -- Extract barcodes using known flanking sequences --

  readname_cell_table <- extract_barcodes(
    filtered_cell_barcode_reads_forward,
    filtered_cell_barcode_reads_reverse,
    seq_before_barcode_1,
    seq_before_barcode_2,
    seq_after_barcode_3,
    verbose
  )
  
  # -- Group barcodes --

  cell_BC1_BC2_BC3 <- group_barcodes(readname_cell_table, verbose)

  # -- Optional file saving (only if output_dir is provided) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(cell_BC1_BC2_BC3,
            file = file.path(output_dir, "cell_to_extracted_BCs.rds"))

    cli::cli_alert_success(
      "Grouped barcodes exported to: {normalizePath(output_dir)}")
  }
}


#' @keywords internal
filter_cell_barcode_reads <- function(
  mapped_cell_barcode_reads,
  min_reads,
  verbose
) {
  if (verbose) cli::cli_inform("Filtering cell barcode reads...")
  # Mark and filter for barcodes that appear more than @min_reads times
  barcodes.keep <-
    names(table(mapped_cell_barcode_reads$cell_barcode))[table(mapped_cell_barcode_reads$cell_barcode)>min_reads]
  filtered_cell_barcode_reads <-
    mapped_cell_barcode_reads[mapped_cell_barcode_reads$cell_barcode%in%barcodes.keep,]

  # Subset of forward reads with a start position between 40 and 43
  filtered_cell_barcode_reads_forward <- 
    filtered_cell_barcode_reads[filtered_cell_barcode_reads$start%in%40:43,]
  # Subset of reverse reads with a start position between 142 and 148
  filtered_cell_barcode_reads_reverse <- 
    filtered_cell_barcode_reads[filtered_cell_barcode_reads$start%in%142:148,]

  # Removing the two reads with unusual start position for forward read
  # and are included for both forward and reverse read
  filtered_cell_barcode_reads_reverse <- filtered_cell_barcode_reads_reverse %>%
    dplyr::filter(!duplicated(readname))

  if (verbose) cli::cli_alert_success("Cell barcode reads filtered.")

  return(list(
    filtered_cell_barcode_reads_forward = filtered_cell_barcode_reads_forward,
    filtered_cell_barcode_reads_reverse = filtered_cell_barcode_reads_reverse
  ))
}


#' @keywords internal
extract_barcodes <- function(
  cell_barcode_reads_forward,
  cell_barcode_reads_reverse,
  seq_before_barcode_1,
  seq_before_barcode_2,
  seq_after_barcode_3,
  verbose
) {
  if (verbose) cli::cli_inform("Extracting barcodes...")
  
  # Filter FORWARD reads that contain both flanking sequences BEFORE the barcode
  cell_barcode_reads_with_seq_before_BC_1_2 <-
    cell_barcode_reads_forward[grepl_seq1_2(
      cell_barcode_reads_forward$sequence,
      seq_before_barcode_1, seq_before_barcode_2), ]

  # Filter REVERSE reads that contain the flanking sequence AFTER the barcode
  cell_barcode_reads_with_seq_after_barcode_3 <-
    cell_barcode_reads_reverse[sapply(
      cell_barcode_reads_reverse$sequence,
      function(x) grepl(seq_after_barcode_3, x)), ]

  # Apply extract_BC1 to each sequence and store results in BC1 column
  cell_barcode_reads_with_seq_before_BC_1_2$BC1 <-
    unlist(sapply(cell_barcode_reads_with_seq_before_BC_1_2$sequence,
                  function(x) extract_BC1(x, seq_before_barcode_1)))

  # Apply extract_BC2 to each sequence and store results in BC2 column
  cell_barcode_reads_with_seq_before_BC_1_2$BC2 <-
    unlist(sapply(cell_barcode_reads_with_seq_before_BC_1_2$sequence,
                  function(x) extract_BC2(x, seq_before_barcode_2)))

  # Create a table with selected columns: readname, cell_barcode, BC1, and BC2
  readname_cell_table <-
    cell_barcode_reads_with_seq_before_BC_1_2 %>%
    dplyr::select(readname, cell_barcode, BC1, BC2)

  # Apply extract_BC3 to each sequence and store results in BC3 column
  cell_barcode_reads_with_seq_after_barcode_3$BC3 <-
    sapply(cell_barcode_reads_with_seq_after_barcode_3$sequence,
           function(x) extract_BC3(x, seq_after_barcode_3))

  # Remove duplicate readnames to ensure unique mapping
  cell_barcode_reads_with_seq_after_barcode_3 <-
    cell_barcode_reads_with_seq_after_barcode_3 %>%
    dplyr::filter(!(duplicated(readname)))

  readname_cell_table <- readname_cell_table %>%
    dplyr::left_join(
      cell_barcode_reads_with_seq_after_barcode_3 %>%
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
  if (verbose) cli::cli_inform("Grouping barcodes...")
  
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



