#' Process Barcodes
#'
#' This function takes barcode reads data and first performs some basic QC
#' filtering. It then extracts BC1 and BC2/BC3s from provided flanking
#' sequences, before grouping cells by cell_barcode and counting reads.
#'
#' @param barcode_amplicon_data Data frame with barcode amplicon info.
#' @param cell_barcodes Vector with extracted read barcodes.
#' @param readnames Vector with extracted readnames.
#' @param min_reads QC threshold: minimum number of reads.
#' @param quality_score_threshold QC threshold: minimum quality score.
#' @param seq_before_barcode_1 BC1 flanking sequence.
#' @param seq_before_barcode_2 BC2 flanking sequence.
#' @param seq_after_barcode_3 BC3 flanking sequence.

#' @return A table linking cells to extracted barcodes.

#' @examples
#' process_barcodes(
#' barcode_amplicon_data = barcode_amplicon_data_subset,
#' cell_barcodes = cell_barcodes_subset,
#' readnames = readnames_subset,
#' min_reads = 2,
#' quality_score_threshold = 30,
#' seq_before_barcode_1 = "GTACGAAGACTG",
#' seq_before_barcode_2 = "CGATACTTGTTCGATCCTTCGC",
#' seq_after_barcode_3 = "CTTGTCTCGCTCCGGG")
#' @importFrom dplyr %>%

#' @export
process_barcodes <- function(
  barcode_amplicon_data,
  cell_barcodes,
  readnames,
  min_reads,
  quality_score_threshold,
  seq_before_barcode_1,
  seq_before_barcode_2,
  seq_after_barcode_3,
  output_dir = NULL
) {
  # -- QC Filtering --

  message("Performing QC Filtering...")

  # Assign readnames as the names of the cell_barcodes vector
  names(cell_barcodes) <- readnames
  # Map cell barcodes to the amplicon data, using readname as the key
  barcode_amplicon_data$cell_barcode <-
    cell_barcodes[barcode_amplicon_data$readname]

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

  message("QC Filtering Complete.")

  # -- Extracting Barcodes 1, 2 and 3 --

  message("Extracting Barcodes 1, 2 and 3...")

  # Filter FORWARD reads that contain both flanking sequences BEFORE the barcode
  barcode_amplicon_data_with_seq_before_BC_1_2 <-
    barcode_amplicon_data_forward[grepl_seq1_2(
      barcode_amplicon_data_forward$sequence,
      seq_before_barcode_1,seq_before_barcode_2),]

  # Filter REVERSE reads that contain the flanking sequence AFTER the barcode
  barcode_amplicon_data_with_seq_after_barcode_3 <-
    barcode_amplicon_data_reverse[sapply(
      barcode_amplicon_data_reverse$sequence,
      function(x) grepl(seq_after_barcode_3,x)),]

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
    dplyr::select(readname,cell_barcode,BC1,BC2)

  # Apply extract_BC3 to each sequence and store results in BC3 column
  barcode_amplicon_data_with_seq_after_barcode_3$BC3 <-
    sapply(barcode_amplicon_data_with_seq_after_barcode_3$sequence,
           function(x) extract_BC3(x, seq_after_barcode_3))

  # Set rownames of above, using readname for easier lookup
  rownames(barcode_amplicon_data_with_seq_after_barcode_3) <-
    barcode_amplicon_data_with_seq_after_barcode_3$readname

  # Remove duplicate readnames to ensure unique mapping
  barcode_amplicon_data_with_seq_after_barcode_3 <-
    barcode_amplicon_data_with_seq_after_barcode_3 %>%
    dplyr::filter(!(duplicated(readname)))

  # Map BC3 values to readname_cell_table using readname as the key
  readname_cell_table$BC3 <- unlist(
    barcode_amplicon_data_with_seq_after_barcode_3[readname_cell_table$readname,"BC3"])

  # Remove rows with any missing (NA) values from the table
  readname_cell_table <-
    readname_cell_table[complete.cases(readname_cell_table),]

  message("Finished Extracting Barcodes 1, 2 and 3.")

  # -- Grouping --

  message("Grouping Barcodes...")

  # Group by cell_barcode, BC1, BC2, and BC3, then count occurrences
  cell_BC1_BC2_BC3 <- readname_cell_table %>%
    dplyr::select(cell_barcode,BC1,BC2,BC3) %>%
    dplyr::group_by_all()%>%
    dplyr::count()

  # Remove rows with any missing values
  cell_BC1_BC2_BC3 <- cell_BC1_BC2_BC3[complete.cases(cell_BC1_BC2_BC3),]

  # Filter for only the rows where BC1 is 20 bases, and BC2 and BC3 are 30 bases
  cell_BC1_BC2_BC3 <- cell_BC1_BC2_BC3 %>%
    dplyr::filter(nchar(BC1) == 20,nchar(BC2) == 30, nchar(BC3) == 30)

  message("Finished Grouping Barcodes.")

  # -- Optional file saving (only if output_dir is provided) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(cell_BC1_BC2_BC3, file = file.path(output_dir, "cell_to_extracted_BCs.rds"))
    message(glue::glue("Grouped Barcodes exported to: {normalizePath(output_dir)}"))
  }

  # Return the grouped data invisibly
  invisible(cell_BC1_BC2_BC3)
}

#' Identifying Closest Whitelisted Barcodes
#'
#' This function identifies the three shortest Hamming distances between barcode
#' and whitelisted barcode, for each of BC1 and BC2/BC3.
#'
#' @param BC1_whitelist Character vector of whitelisted BC1 barcodes.
#' @param BC2_BC3_whitelist Character vector of whitelisted BC2/BC3 barcode pairs.
#' @param cell_to_extracted_BCs A data frame linking cells to extracted barcodes.

#' @return A list of barcode matches for BC1 and a list of barcode matches for
#' BC2/BC3.

#' @examples
#' identify_whitelisted_barcodes(
#' BC1_whitelist = BC1_whitelist_subset,
#' BC2_BC3_whitelist = BC2_BC3_whitelist,
#' cell_to_extracted_BCs = cell_to_extracted_BCs_subset
#' )
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom stringdist stringdist

#' @export
identify_whitelisted_barcodes <- function(
    BC1_whitelist,
    BC2_BC3_whitelist,
    cell_to_extracted_BCs,
    output_dir = NULL
) {
  # -- Identify the Three BC1s with the Shortest Hamming Distance

  # Get unique BC1 sequences from the dataset
  cell_BC1_BC2_BC3 <- cell_to_extracted_BCs
  unique_BC_1 <- unique(cell_BC1_BC2_BC3$BC1)
  message(glue::glue("Number of Unique BC1s: {length(unique_BC_1)}"))

  # Loop through each unique BC1 sequence and find the 3 closest whitelist
  # barcodes to the current BC1 using Hamming distance
  cli_progress_bar("Identifying BC1s", total = length(unique_BC_1))
  list_BC1 <- list()
  for (j in 1:length(unique_BC_1)){
    list_BC1[[j]] <- identify_barcode_ID(unique_BC_1[j],BC1_whitelist)
    cli_progress_update() # Update progress
  }
  cli_progress_done() # Close the progress bar

  # Assign BC1 sequences as names to the list entries
  names(list_BC1) <- unique_BC_1

  # -- Identify the Three BC2/BC3s with the Shortest Hamming Distance

  # Get unique BC2 and BC3 sequences from the dataset
  unique_BC_2_3 <- unique(c(cell_BC1_BC2_BC3$BC2,cell_BC1_BC2_BC3$BC3))
  message(glue::glue("Number of Unique BC2_3s: {length(unique_BC_2_3)}"))

  # Loop through each unique BC2/BC3 sequence and find the 3 closest whitelist
  # barcodes to the current BC2/BC3 using Hamming distance
  cli_progress_bar("Identifying BC2/BC3s", total = length(unique_BC_2_3))
  list_BC2_3 <- list()
  for (j in 1:length(unique_BC_2_3)){
    list_BC2_3[[j]] <- identify_barcode_ID(unique_BC_2_3[j],BC2_BC3_whitelist)
    cli_progress_update() # Update progress
  }
  cli_progress_done() # Close the progress bar

  # Assign BC2/BC3 sequences as names to the list entries
  names(list_BC2_3) <- unique_BC_2_3


  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(list_BC1, file = file.path(output_dir, "list_BC1.rds"))
    saveRDS(list_BC2_3, file = file.path(output_dir, "list_BC2_3.rds"))
    message(glue::glue("Results saved to: {normalizePath(output_dir)}"))
  }

  invisible(list(BC1_matches = list_BC1, BC2_BC3_matches = list_BC2_3))
}




