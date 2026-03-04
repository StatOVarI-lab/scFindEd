#' Process Reads
#'
#' This function takes the mapped reads for each amplicon (.rds), and both the 
#' barcoded HDR outcomes and the cell-to-UTR barcode system, to perform the following:
#' - Filters for only barcodes found in the cell-to-UTR barcode system.
#' - Pairs and merges reads (R1 & R2).
#' - Generates a complete list of possible SNP and indel variants.
#' - Uses this list to assign a genotype to each read pair.
#' - Collects and counts read-pair genotypes across cell barcodes.
#' 
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param flank_len Integer denoting the length of the flanks used to locate and 
#' identify the editing window within a read.
#' @param editing_window_extension Integer denoting the number of bases to extend
#' each editing window by.
#' @param n_cores Integer denoting the number of cores to be utilised by the 
#' internal variant generator.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return An exported summary CSV file for each amplicon, named 
#' "reads_classified_no_extension_Element_HDR_{HDR_number}.csv"

#' @importFrom dplyr %>%

#' @export
process_reads <- function(
  source_dir,
  flank_len = 8,
  editing_window_extension = 0,
  n_cores = 1,
  output_dir = NULL,
  verbose
) {
  # -- Import Static Data Files --
  
  invisible(list2env(import_process_reads_data(
    source_dir,
    verbose
  ), envir = environment()))
  
  number_of_HDRs <- nrow(HDR_barcoded_outcomes)
  
  # -- Loop through each HDR --
  
  for (i in 1:number_of_HDRs) {
    if (verbose) {
      cli::cli_h2("Processing HDR {i} of {number_of_HDRs}")
    }
    
    AmpID <- HDR_barcoded_outcomes$AMPLICON_NAME[i]
    
    # -- Load Mapped Reads for AmpID --
    
    mapped_reads <- import_barcoded_mapped_reads(source_dir, AmpID, verbose)
    if ("cell_barcode" %in% colnames(mapped_reads)) {
      colnames(mapped_reads)[colnames(mapped_reads) == "cell_barcode"] <- "barcode"
    }
    
    # -- Filter for ONLY Barcodes Recorded in cell_to_UTR_barcode_system --
    
    mapped_reads_barcoded <- mapped_reads[mapped_reads$barcode %in% barcode_system$barcode,]
    mapped_reads_barcoded <- mapped_reads_barcoded[order(mapped_reads_barcoded$readname), ]
    
    # -- Pair Reads --
    
    mapped_reads_paired_count <- pair_reads(
      mapped_reads_barcoded,
      i,
      verbose
    )
    
    # -- Merge R1 and R2 --
    
    mapped_reads_paired_count <-
      cbind(mapped_reads_paired_count,
            merge_r1_r2(mapped_reads_paired_count$sequence_R1,
                        mapped_reads_paired_count$sequence_R2,
                        min_overlap = 2,
                        verbose))
    
    # Filter out rows with NA in $merged_sequence
    mapped_reads_paired_count <- mapped_reads_paired_count %>%
      dplyr::filter(!is.na(merged_sequence))
    
    # -- Classify Genotypes
    
    # Extract templates for HDR
    WT_editing_window <- HDR_barcoded_outcomes$amplicon_seq_in_window[i]
    HDR_editing_window <- HDR_barcoded_outcomes$HDR_seq_in_window[i]
    amplicon_sequence <- HDR_barcoded_outcomes$Amplicon.sequence[i]
    
    # Generate possible variants
    df_variants <- generate_all_variants(
      WT_editing_window,
      HDR_editing_window,
      amplicon_sequence,
      all_reads = mapped_reads_paired_count,
      flank_len,
      editing_window_extension,
      verbose
    )
    
    # Map genotypes
    mapped_reads_paired_count <- map_genotypes(
      mapped_reads_paired_count,
      df_variants,
      sequence_col = "merged_sequence",
      n_cores, 
      verbose
    )
    
    # -- Export Processed Reads --
    
    export_processed_reads(mapped_reads_paired_count, i, output_dir, verbose)
    
    # Remove datasets for HDR i from memory
    rm(mapped_reads, mapped_reads_barcoded, 
       mapped_reads_paired_count, df_variants)
    gc()
  }
}


#' @keywords internal
import_process_reads_data <- function(
  source_dir,
  verbose
) {
  # -- Read in the HDR outcome file --
  
  HDR_barcoded_outcomes <- import_HDR_barcoded_data(
    source_dir,
    "HDR_file_barcoded_with_amplicon_sequences.csv",
    verbose
  )
  
  # -- Import cell to UTR barcode system --
  
  barcode_system <- import_cell_to_UTR_barcode_system(
    source_dir,
    "cell_to_UTR_barcode_system.csv",
    verbose
  )
  barcode_system <- barcode_system %>%
    dplyr::rename(barcode = cell_barcode)
  
  return(list(
    HDR_barcoded_outcomes = HDR_barcoded_outcomes,
    barcode_system = barcode_system
  ))
}


#' @keywords internal
pair_reads <- function(
  mapped_reads_barcoded,
  HDR_no,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Pairing HDR {HDR_no} reads...",
                                      msg_done = "HDR {HDR_no} reads paired.")
  
  # Pair reads (R1 and R2)
  paired_idx <- seq(1, nrow(mapped_reads_barcoded) - 1, by = 2)
  mapped_reads_paired <- cbind(
    mapped_reads_barcoded[paired_idx,     c("readname", "amplicon_name",
                                            "barcode",  "start",
                                            "cigar",    "sequence")],
    mapped_reads_barcoded[paired_idx + 1, c("readname", "amplicon_name",
                                            "barcode", "start",
                                            "cigar",   "sequence")]
  )
  colnames(mapped_reads_paired) <- c(
    "readname_R1", "amplicon_name",  "barcode",
    "start_R1",    "cigar_R1",       "sequence_R1",
    "readname_R2", "amplicon_name2", "barcode2",
    "start_R2",    "cigar_R2",       "sequence_R2"
  )
  
  # Keep only pairs where readnames match
  mapped_reads_paired <-
    mapped_reads_paired[mapped_reads_paired$readname_R1 ==
                          mapped_reads_paired$readname_R2, ]
  
  # Remove duplicate + unnecessary columns
  mapped_reads_paired <-
    mapped_reads_paired[, c("amplicon_name", "barcode",  "start_R1",
                            "cigar_R1",      "start_R2", "cigar_R2",
                            "sequence_R1",   "sequence_R2")]
  
  # Ensure R1 is always the lower start
  swap_idx <- mapped_reads_paired$start_R1 >
    mapped_reads_paired$start_R2
  
  if (any(swap_idx)) {
    temp <- mapped_reads_paired[swap_idx, ]
    temp <- temp %>%
      dplyr::mutate(
        start_R1 = start_R2,
        start_R2 = start_R1,
        sequence_R1 = sequence_R2,
        sequence_R2 = sequence_R1,
        cigar_R1 = cigar_R2,
        cigar_R2 = cigar_R1
      )
    mapped_reads_paired[swap_idx, ] <- temp
  }
  
  # Count pairs
  mapped_reads_paired_count <- mapped_reads_paired %>%
    dplyr::group_by_all() %>%
    dplyr::tally(name = "n") %>%
    dplyr::ungroup()
  
  # Convert to data frame and apply rownames
  mapped_reads_paired_count <- as.data.frame(mapped_reads_paired_count)
  rownames(mapped_reads_paired_count) <- seq_len(nrow(mapped_reads_paired_count))
  
  if (verbose) cli::cli_process_done()
  
  return(mapped_reads_paired_count)
}


#' @keywords internal
export_processed_reads <- function(
  mapped_reads_paired_count, 
  amplicon_gRNA_number,
  output_dir, 
  verbose
) {
  df_to_export <- mapped_reads_paired_count[, -1]
  df_to_export <- df_to_export %>% 
    dplyr::filter(!(is.na(merged_sequence))) %>% 
    dplyr::select(barcode, genotype)
  
  df_to_export_counts <- df_to_export %>% 
    dplyr::group_by(barcode, genotype) %>% 
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>% 
    dplyr::mutate(pct = 100 * n / sum(n)) %>% 
    dplyr::ungroup()
  
  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    file_save_name <- file.path(output_dir, paste0(
      "reads_classified_no_extension_Element_HDR_", amplicon_gRNA_number, ".csv"))
    
    readr::write_csv(df_to_export_counts, file_save_name)
    
    cli::cli_alert_success("HDR {amplicon_gRNA_number} results saved to: {normalizePath(output_dir)}.")
  }
}

