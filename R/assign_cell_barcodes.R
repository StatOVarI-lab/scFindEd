#' Assign Cell Barcodes
#' 
#' This function takes the mapped reads for each amplicon (.tsv.gz), and both the 
#' barcodes-readnames FASTQ and amplicons list, to perform the following:
#' - Filter reads for a minimum quality.
#' - Pairs reads.
#' - Merges cell barcodes into the mapped reads data, matching by readname.
#' 
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param sample_name String containing the file name structure for the
#' experiment sample.
#' @param min_quality_threshold Integer denoting the filtration threshold for 
#' read quality.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return An exported .rds file for each amplicon, containing the mapped reads 
#' each tagged with the corresponding cell barcode, named "{AmpID}_info.rds".

#' @importFrom dplyr %>%
#' @import data.table

#' @export
assign_cell_barcodes <- function(
  source_dir, 
  sample_name,
  min_quality_threshold,
  output_dir = NULL,
  verbose
) {
  # -- Import Data Files --
  
  invisible(list2env(import_assign_cell_barcodes_data(
    source_dir,
    sample_name,
    verbose
  ), envir = environment()))
  
  # -- Loop through amplicons --
  
  for (AmpID in amplicons) {
    if (verbose) cat("\n")
    if (verbose) cli::cli_inform("Processing amplicon: {AmpID}")
    
    # Import Dataset for AmpID
    mapped_reads <- import_mapped_reads(source_dir, AmpID, 
                                        min_quality_threshold, verbose)
    
    # Check for empty file
    if (nrow(mapped_reads) == 0) {
      if (verbose) cli::cli_alert_warning(
        "No mapped reads for {AmpID}. Skipping to next amplicon.")
      next
    }
    
    # Convert mapped_reads to a data.table and assign readname as the key
    data.table::setDT(mapped_reads)
    data.table::setkey(mapped_reads, readname)
    
    # Pair mapped_reads
    if (verbose) cli::cli_inform(" - Pairing mapped reads...")
    
    keep_readnames <- mapped_reads[ , .(N = .N), by = .(readname) ][ N == 2L, readname ] 
    
    mapped_reads <- mapped_reads[readname %chin% keep_readnames]
    
    # Merge barcodes into mapped reads
    mapped_reads <- barcodes_readnames[mapped_reads, nomatch = 0]
    
    # -- Optional file saving (only if output_dir is specified) --
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      
      file_save_name <- file.path(output_dir, paste0(AmpID, "_info.rds"))
      
      saveRDS(mapped_reads, file_save_name)
      
      cli::cli_alert_success("Results saved to: {normalizePath(output_dir)}.")
    }
    
    # Remove Dataset for AmpID from Memory
    rm(mapped_reads)
    gc()
  }
}

#' @keywords internal
import_assign_cell_barcodes_data <- function(
  source_dir,
  sample_name,
  verbose
) {
  # -- Read in barcodes_readnames tsv file --
  
  barcodes_readnames_file_name <- paste0(sample_name, "_R1_001_barcodes_readnames.tsv.gz")
  barcodes_readnames <- import_barcodes_readnames(
    source_dir,
    barcodes_readnames_file_name,
    verbose
  )
  
  # Convert barcodes_readnames to a data.table and assign readname as the key
  if (verbose) cli::cli_inform(" - Converting barcodes_readnames to data.table...")
  data.table::setDT(barcodes_readnames)
  data.table::setkey(barcodes_readnames, readname)
  
  # -- Read in amplicon list --
  
  amplicon_list_file_name <- "amplicon_list.txt"
  amplicon_list <- import_amplicon_list(
    source_dir,
    amplicon_list_file_name,
    verbose
  )
  
  amplicons <- amplicon_list$X1
  
  return(list(
    barcodes_readnames = barcodes_readnames,
    amplicons = amplicons
  ))
}
