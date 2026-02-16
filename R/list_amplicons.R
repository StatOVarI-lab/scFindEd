

list_amplicons <- function(
  source_dir,
  design_summary_file_name,
  output_dir = NULL,
  verbose
) {
  # Import design summary file
  design_summary_file <- import_design_summary_info(
    source_dir, 
    design_summary_file_name, 
    verbose
  )
  
  # Filter for only rows with amplicon IDs and an amplicon sequence
  amplicon_data <- design_summary_file[
    design_summary_file$AmpID != "" & !is.na(design_summary_file$AmpID) &
    design_summary_file$`Amplicon sequence` != "" & !is.na(design_summary_file$`Amplicon sequence`),
  ]
  
  # Extract only the list of amplicon IDs
  amplicon_list <- amplicon_data$AmpID
  
  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    file_save_name <- file.path(output_dir, "amplicon_list.txt")
    
    writeLines(amplicon_list, file_save_name)
    
    cli::cli_alert_success("Amplicon list saved to: {normalizePath(output_dir)}.")
  }
  
  return(amplicon_list)
}