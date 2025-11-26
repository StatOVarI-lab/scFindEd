#' Merge HDR Barcodes with Design Summary File
#'
#' This function takes a HDR barcoded csv file path and a design summary file
#' path, matches amplicons and merges the two files according to the matched
#' amplicons.
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param HDR_barcoded_file_name String containing the file name structure for
#' the HDR barcoded file.
#' @param design_summary_file_name String containing the file name structure for
#' the design summary file.
#' @param output_dir String containing the relative path to the directory that
#' the merged data file will be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return A data frame containing the merged amplicon summary data. Also
#' exported as a .csv file to the defined output_dir.

#' @examples
#' merge_HDR_amplicon_sequences(
#' source_dir = "data",
#' HDR_barcoded_file_name = "HDR_file_barcoded.csv",
#' design_summary_file_name = "6915_design_summary_extended_information.csv",
#' output_dir = "data",
#' verbose = TRUE
#' )

#' importFrom("base", "sequence")

#' @export
merge_HDR_amplicon_sequences <- function(
  source_dir,
  HDR_barcoded_file_name,
  design_summary_file_name,
  output_dir = NULL,
  verbose
) {
  # -- Importing the Data --

  imported_data <- import_merge_HDR_data(
    source_dir,
    HDR_barcoded_file_name,
    design_summary_file_name,
    verbose
  )

  invisible(list2env(imported_data, envir = environment()))

  # -- Data Exploration --

  unique_amplicons <- amplicon_data_exploration(HDR_df, design_df, verbose)

  invisible(list2env(unique_amplicons, envir = environment()))

  # -- Matching Analysis --

  matches <- amplicon_matching_analysis(
    unique_HDR_amplicons,
    unique_design_amplicons,
    verbose
  )

  # -- Amplicon Sequence Availability --

  if (verbose) amplicon_sequence_availability(design_df, matches)

  # -- Merge Data --

  merged_df <- perform_amplicon_merge(HDR_df, design_df, verbose)

  # -- Merge Validation --

  validate_amplicon_merge(merged_df, HDR_df, verbose)

  # -- Display Summary by Amplicon Name --

  if (verbose) summarise_by_amplicon_name(merged_df)

  # -- Display Merged Data Sample --

  if (verbose) display_merged_amplicon_data_sample(merged_df)

  # -- Export Merged Data --

  export_merged_amplicon_data(merged_df, output_dir, verbose)

  # -- Summarise Outcomes --

  if (verbose) summarise_amplicon_merge_outcomes(merged_df)

  return(merged_df)
}


#' @keywords internal
import_merge_HDR_data <- function(
  source_dir,
  HDR_barcoded_file_name,
  design_summary_file_name,
  verbose
) {
  if (verbose) cli::cli_h2("DATA IMPORT")

  # -- Load the HDR barcoded file --
  HDR_df <- import_HDR_barcoded_data(
    source_dir,
    HDR_barcoded_file_name,
    verbose
  )

  if (verbose) cli::cli_alert_info(
    "HDR file shape: {nrow(HDR_df)} rows, {ncol(HDR_df)} columns.")
  if (verbose) cli::cli_alert_info(
    "HDR file columns: {paste(colnames(HDR_df), collapse = ', ')}.")

  # -- Load the design summary file --
  design_df <- import_design_summary_info(
    source_dir,
    design_summary_file_name,
    verbose
  )

  if (verbose) cli::cli_alert_info(
    "Design summary shape: {nrow(design_df)} rows, {ncol(design_df)} columns.")

  # -- Check that all required columns for merging are present --
  required_cols <- list(
    HDR_df = c("AMPLICON_NAME"),
    design_df = c("AmpID", "Amplicon sequence")
  )

  if (check_required_columns(required_cols) == FALSE) {
    cli::cli_abort("Missing columns detected - function aborted.")
  }

  if (verbose) cli::cli_text("\n")
  cli::cli_alert_success("Data successfully imported.")

  return(list(
    HDR_df = HDR_df,
    design_df = design_df
  ))
}


#' @keywords internal
amplicon_data_exploration <- function(
  HDR_df,
  design_df,
  verbose
) {
  if (verbose) cli::cli_h2("DATA EXPLORATION")

  # Assess unique amplicons in the HDR barcoded file
  unique_HDR_amplicons <- unique(HDR_df$AMPLICON_NAME)

  if (verbose) cli::cli_alert_info(glue::glue(
    "Number of unique AMPLICON_NAMEs in HDR barcoded file: ",
    "{length(unique_HDR_amplicons)}."
  ))
  if (verbose) cli::cli_alert_info(glue::glue(
    "HDR amplicons: ",
    "{paste(sort(unique_HDR_amplicons), collapse = ', ')}."
  ))
  if (verbose) cli::cli_text("\n")

  # Assess unique amplicons in the design summary file
  unique_design_amplicons <- unique(design_df$AmpID)

  if (verbose) cli::cli_alert_info(glue::glue(
    "Number of unique AmpIDs in design summary file: ",
    "{length(unique_design_amplicons)}."
  ))
  if (verbose) cli::cli_alert_info(glue::glue(
    "Design amplicons: ",
    "{paste(sort(unique_design_amplicons), collapse = ', ')}."
  ))

  return(list(
    unique_HDR_amplicons = unique_HDR_amplicons,
    unique_design_amplicons = unique_design_amplicons
  ))
}


#' @keywords internal
amplicon_matching_analysis <- function(
  unique_HDR_amplicons,
  unique_design_amplicons,
  verbose
) {
  if (verbose) cli::cli_h2("MATCHING ANALYSIS")

  # Find matches
  matches <- intersect(unique_HDR_amplicons, unique_design_amplicons)
  HDR_only <- setdiff(unique_HDR_amplicons, unique_design_amplicons)
  design_only <- setdiff(unique_design_amplicons, unique_HDR_amplicons)

  # Report matches to the user

  if (verbose) cli::cli_alert_info(
    "Number of amplicon matches found: {length(matches)}.")
  if (verbose) cli::cli_alert_info(
    "Matches: {paste(sort(matches), collapse = ', ')}.")
  if (verbose) cli::cli_text("\n")

  # Report amplicons uniquely found in HDR barcoded file
  if (verbose) cli::cli_alert_info(
    "Unique to HDR barcoded file: {length(HDR_only)}.")
  if (length(HDR_only) > 0) {
    if (verbose) cli::cli_alert_info(
      "Unique amplicons: {paste(sort(HDR_only), collapse = ', ')}.")
  }
  if (verbose) cli::cli_text("\n")

  # Report amplicons uniquely found in design summary file
  if (verbose) cli::cli_alert_info(
    "Unique to design summary file: {length(design_only)}.")
  if (length(design_only) > 0) {
    if (verbose) cli::cli_alert_info(
      "Unique amplicons: {paste(sort(design_only), collapse = ', ')}.")
  }

  return(matches)
}


#' @keywords internal
amplicon_sequence_availability <- function(
  design_df,
  matches
) {
  cli::cli_h2("AMPLICON SEQUENCE AVAILABILITY")

  # Filter design_df to only matched amplicon IDs, and select sample of columns
  matching_design_data <- design_df %>%
    dplyr::filter(AmpID %in% matches) %>%
    dplyr::select(AmpID, Name, `Amplicon sequence`)

  cli::cli_alert_info("Amplicon sequence status for matching IDs:")

  # Iterate through matched amplicons, and check for presence of a sequence
  for (amp_id in sort(matches)) {
    # Filter for iterated amplicon and extract the sequence
    amp_seq <- matching_design_data %>%
      dplyr::filter(AmpID == amp_id) %>%
      dplyr::pull(`Amplicon sequence`)

    # Checks whether the 1st element of amp_seq is a non-empty, non-NA sequence
    has_sequence <- !is.na(amp_seq[1]) && nchar(trimws(amp_seq[1])) > 0

    # Updates user on the status of each amplicon sequence
    amp_status_text <- if (has_sequence) "Has sequence" else "No sequence"

    cli::cli_alert_success("{amp_id}: {amp_status_text}.")
  }
}


#' @keywords internal
perform_amplicon_merge <- function(
  HDR_df,
  design_df,
  verbose
) {
  if (verbose) cli::cli_h2("PERFORMING MERGE")

  # Subset and reformat design_df in preparation for merge
  design_subset <- design_df %>%
    dplyr::select(AmpID, `Amplicon sequence`) %>%
    dplyr::rename(AMPLICON_NAME = AmpID)

  # Perform the merge (left join to keep all HDR entries)
  merged_df <- HDR_df %>%
    dplyr::left_join(design_subset, by = "AMPLICON_NAME")

  cli::cli_alert_success("Successfully merged files.")

  # Report merge summary to the user
  if (verbose) cli::cli_alert_info(
    "Original HDR barcoded file rows: {nrow(HDR_df)}.")
  if (verbose) cli::cli_alert_info(
    "Merged file rows: {nrow(merged_df)}.")
  if (verbose) cli::cli_alert_info(
    "Columns in merged file: {paste(colnames(merged_df), collapse = ', ')}.")

  return(merged_df)
}


#' @keywords internal
validate_amplicon_merge <- function(
  merged_df,
  HDR_df,
  verbose
) {
  if (verbose) cli::cli_h2("VALIDATING MERGE")

  # Calculation of validation statistics
  has_sequence <- !is.na(merged_df$`Amplicon sequence`)
  entries_with_sequences <- sum(has_sequence)
  entries_without_sequences <- sum(!has_sequence)
  success_rate <- (entries_with_sequences / nrow(merged_df)) * 100

  # -- Update user on validation statistics --

  if (verbose == FALSE) cli::cli_text("\n")

  cli::cli_text("Merge validation summary:")

  # Provide update on preservation of original HDR entries
  if (nrow(merged_df) == nrow(HDR_df)) {
    cli::cli_alert_success("All original HDR entries preserved.")
  } else {
    cli::cli_alert_danger("Some original HDR entries NOT preserved.")
  }

  # Provide update on any introduction of duplicate rows
  if (nrow(merged_df) == nrow(dplyr::distinct(merged_df))) {
    cli::cli_alert_success("No duplicate rows introduced.")
  } else {
    cli::cli_alert_danger("Duplicate rows introduced.")
  }

  # Provide update on the number of amplicon sequences with entries
  if (entries_with_sequences == nrow(merged_df)) {
    cli::cli_alert_success(glue::glue(
      "Entries with amplicon sequences: ",
      "{entries_with_sequences} out of {nrow(merged_df)}."
    ))
  } else if (entries_with_sequences == 0) {
    cli::cli_alert_danger(glue::glue(
      "Entries with amplicon sequences: ",
      "{entries_with_sequences} out of {nrow(merged_df)}."
    ))
  } else {
    cli::cli_alert_warning(glue::glue(
      "Entries with amplicon sequences: ",
      "{entries_with_sequences} out of {nrow(merged_df)}."
    ))
  }

  # Provide update on the number of entries without amplicon sequences
  if (entries_without_sequences == 0) {
    cli::cli_alert_success(
      "Entries without amplicon sequences: {entries_without_sequences}.")
  } else {
    cli::cli_alert_warning(
      "Entries without amplicon sequences: {entries_without_sequences}.")
  }

  # Report success rate %
  if (success_rate == 100) {
    cli::cli_alert_success("Success rate: {sprintf('%.1f', success_rate)}%.")
  } else if (success_rate == 0) {
    cli::cli_alert_danger("Success rate: {sprintf('%.1f', success_rate)}%.")
  } else {
    cli::cli_alert_warning("Success rate: {sprintf('%.1f', success_rate)}%.")
  }
}


#' @keywords internal
summarise_by_amplicon_name <- function(
  merged_df
) {
  cli::cli_h2("SEQUENCE STATUS BY AMPLICON NAME")

  # Group by amplicon name, and summarise the status of each amplicon
  sequence_status <- merged_df %>%
    dplyr::group_by(AMPLICON_NAME) %>%
    dplyr::summarise(
      total_entries = dplyr::n(),
      entries_with_sequence = sum(!is.na(`Amplicon sequence`)),
      has_sequence = any(!is.na(`Amplicon sequence`)),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(AMPLICON_NAME)

  # Print tibble of the grouped and summarised amplicons
  print(sequence_status)
}


#' @keywords internal
display_merged_amplicon_data_sample <- function(
  merged_df
) {
  cli::cli_h2("SAMPLE OF MERGED DATA")

  # Filter for first 3 results with valid amplicon sequences
  sample_with_sequences <- merged_df %>%
    dplyr::filter(!is.na(`Amplicon sequence`)) %>%
    dplyr::slice_head(n = 3)

  # As long as there exist sample sequences in the merged data...
  if (nrow(sample_with_sequences) > 0) {
    # For each sample, display the contents of each column
    for (i in 1:nrow(sample_with_sequences)) {
      row <- sample_with_sequences[i, ] # Extract entire ith sample
      amp_seq <- row$`Amplicon sequence` # Extract ith amplicon sequence

      cli::cli_text("AMPLICON_NAME: {row$AMPLICON_NAME}.")
      cli::cli_alert_info("sgRNA_SEQUENCE: {row$sgRNA_SEQUENCE}.")
      cli::cli_alert_info("SNP_change: {row$SNP_change}.")

      # Limit printing of amplicon sequence to 50 bases
      if (nchar(amp_seq) > 50) {
        cli::cli_alert_info(glue::glue(
          "Amplicon sequence (first 50 chars): ",
          "{substr(amp_seq, 1, 50)}..."
        ))
      } else {
        cli::cli_alert_info("Amplicon sequence: {amp_seq}.")
      }
      if (i < 3) cli::cli_text("\n")
    }
  }
}


#' @keywords internal
export_merged_amplicon_data <- function(
  merged_df,
  output_dir,
  verbose
) {
  if (verbose) cli::cli_h2("EXPORTING MERGED DATA")

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(merged_df, file = file.path(
      output_dir, "HDR_file_barcoded_with_amplicon_sequences.csv"))

    # Report data export outcome to the user
    cli::cli_alert_success("Merged data exported to: {normalizePath(output_dir)}.")

    if (verbose) cli::cli_alert_info("Total rows: {nrow(merged_df)}.")
    if (verbose) cli::cli_alert_info("Total columns: {ncol(merged_df)}.")

  } else {
    cli::cli_alert_danger(
      "No output directory specified - merged data could not be exported.")
  }
}


#' @keywords internal
summarise_amplicon_merge_outcomes <- function(
  merged_df
) {
  cli::cli_h2("MERGED DATA STRUCTURE")

  # Print each column name to the user
  cli::cli_ol()
  for (i in seq_along(colnames(merged_df))) {
    cli::cli_li("{colnames(merged_df)[i]}.")
  }
  cli::cli_end()
}
