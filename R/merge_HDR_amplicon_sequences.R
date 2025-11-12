#' Merge HDR Barcodes with Design Summary File
#'
#' This function takes a HDR barcoded csv file path and a design summary file
#' path, matches amplicons and merges the two files according to the matched
#' amplicons.
#'
#' @param HDR_barcoded_file_path String containing the relative path to the HDR
#' barcoded file.
#' @param design_summary_file_path String containing the relative path to the
#' design summary file.
#' @param output_dir String containing the relative path to the directory that
#' the merged data file will be exported to.

#' @return A data frame containing the merged amplicon summary data. Also
#' exported as a .csv file to the defined output_dir.

#' @examples
#' merge_HDR_amplicon_sequences(
#' HDR_barcoded_file_path = "data/HDR_file_barcoded.csv",
#' design_summary_file_path = "data/6915_design_summary_extended_information.csv",
#' output_dir = "data"
#' )

#' @importFrom readr read_csv write_csv
#' @importFrom dplyr %>% filter select pull rename left_join distinct group_by
#' summarise n arrange slice_head
#' @importFrom glue glue

#' @export
merge_HDR_amplicon_sequences <- function(
  HDR_barcoded_file_path,
  design_summary_file_path,
  output_dir = NULL
) {
  # -- Importing the Data --

  message("-- DATA IMPORT --")

  # Load the HDR barcoded file
  HDR_df <- read_csv(HDR_barcoded_file_path, show_col_types = FALSE)
  message(glue::glue(
    "HDR file shape: {nrow(HDR_df)} rows, {ncol(HDR_df)} columns."
  ))
  message(glue::glue(
    "HDR file columns: {paste(colnames(HDR_df), collapse = ', ')}."
  ))

  # Load the design summary file
  design_df <- read_csv(design_summary_file_path, show_col_types = FALSE)
  message(glue::glue(
    "Design summary shape: {nrow(design_df)} rows, {ncol(design_df)} columns."
  ))

  # Check that all required columns for merging are present --
  required_cols <- list(
    HDR_df = c("AMPLICON_NAME"),
    design_df = c("AmpID", "Amplicon sequence")
  )

  if (check_required_columns(required_cols) == FALSE) {
    stop("Missing columns detected - function aborted.")
  }

  # -- Data Exploration --

  unique_amplicons <- amplicon_data_exploration(HDR_df, design_df)
  unique_HDR_amplicons <- unique_amplicons[[1]]
  unique_design_amplicons <- unique_amplicons[[2]]

  # -- Matching Analysis --

  matches <- amplicon_matching_analysis(
    unique_HDR_amplicons, unique_design_amplicons)

  # -- Amplicon Sequence Availability --

  amplicon_sequence_availability(design_df, matches)

  # -- Merge Data --

  merged_df <- perform_amplicon_merge(HDR_df, design_df)

  # -- Merge Validation --

  validate_amplicon_merge(merged_df, HDR_df)

  # -- Display Summary by Amplicon Name --

  summarise_by_amplicon_name(merged_df)

  # -- Display Merged Data Sample --

  display_merged_amplicon_data_sample(merged_df)

  # -- Export Merged Data --

  export_merged_amplicon_data(merged_df, output_dir)

  # -- Summarise Outcomes --

  summarise_amplicon_merge_outcomes(merged_df)

  # -- Merge Completed --

  message("\n-- MERGE COMPLETE --\n")

  return(merged_df)
}



#' @keywords internal
amplicon_data_exploration <- function(
  HDR_df,
  design_df
) {
  message("\n-- DATA EXPLORATION --\n")

  # Assess unique amplicons in the HDR barcoded file
  unique_HDR_amplicons <- unique(HDR_df$AMPLICON_NAME)

  message(glue::glue(
    "Number of unique AMPLICON_NAMEs in HDR barcoded file: ",
    "{length(unique_HDR_amplicons)}."
  ))
  message(glue::glue(
    "HDR amplicons: ",
    "{paste(sort(unique_HDR_amplicons), collapse = ', ')}.\n\n"
  ))

  # Assess unique amplicons in the design summary file
  unique_design_amplicons <- unique(design_df$AmpID)

  message(glue::glue(
    "Number of unique AmpIDs in design summary file: ",
    "{length(unique_design_amplicons)}."
  ))
  message(glue::glue(
    "Design amplicons: ",
    "{paste(sort(unique_design_amplicons), collapse = ', ')}."
  ))

  return(list(unique_HDR_amplicons, unique_design_amplicons))
}


#' @keywords internal
amplicon_matching_analysis <- function(
  unique_HDR_amplicons,
  unique_design_amplicons
) {
  message(glue::glue("\n\n-- MATCHING ANALYSIS --\n\n"))

  # Find matches
  matches <- intersect(unique_HDR_amplicons, unique_design_amplicons)
  HDR_only <- setdiff(unique_HDR_amplicons, unique_design_amplicons)
  design_only <- setdiff(unique_design_amplicons, unique_HDR_amplicons)

  # Report matches to the user
  message(glue::glue("Number of amplicon matches found: {length(matches)}."))
  message(glue::glue("Matches: {paste(sort(matches), collapse = ', ')}."))

  # Report amplicons uniquely found in HDR barcoded file
  message(glue::glue("\n\nUnique to HDR barcoded file: {length(HDR_only)}."))
  if (length(HDR_only) > 0) {
    message(glue::glue(
      "Unique amplicons: {paste(sort(HDR_only), collapse = ', ')}."
    ))
  }

  # Report amplicons uniquely found in design summary file
  message(glue::glue("\n\nUnique to design summary file: {length(design_only)}."))
  if (length(design_only) > 0) {
    message(glue::glue(
      "Unique amplicons: {paste(sort(design_only), collapse = ', ')}."
    ))
  }

  return(matches)
}


#' @keywords internal
amplicon_sequence_availability <- function(
  design_df,
  matches
) {
  message("\n-- AMPLICON SEQUENCE AVAILABILITY --\n")

  # Filter design_df to only matched amplicon IDs, and select sample of columns
  matching_design_data <- design_df %>%
    dplyr::filter(AmpID %in% matches) %>%
    dplyr::select(AmpID, Name, `Amplicon sequence`)

  message("Amplicon sequence status for matching IDs:")

  # Iterate through matched amplicons, and check for presence of a sequence
  for (amp_id in sort(matches)) {
    # Filter for iterated amplicon and extract the sequence
    amp_seq <- matching_design_data %>%
      dplyr::filter(AmpID == amp_id) %>%
      dplyr::pull(`Amplicon sequence`)

    # Checks whether the 1st element of amp_seq is a non-empty, non-NA sequence
    has_sequence <- !is.na(amp_seq[1]) && nchar(trimws(amp_seq[1])) > 0

    # Updates user on the status of each amplicon sequence
    amp_status_symbol <- if (has_sequence) "\u2713" else "\u2717"
    amp_status_text <- if (has_sequence) "Has sequence" else "No sequence"
    message(glue::glue("{amp_status_symbol} {amp_id}: {amp_status_text}."))
  }
}


#' @keywords internal
perform_amplicon_merge <- function(
  HDR_df,
  design_df
) {
  message("\n-- PERFORMING MERGE --\n")

  # Subset and reformat design_df in preparation for merge
  design_subset <- design_df %>%
    dplyr::select(AmpID, `Amplicon sequence`) %>%
    dplyr::rename(AMPLICON_NAME = AmpID)

  # Perform the merge (left join to keep all HDR entries)
  merged_df <- HDR_df %>%
    dplyr::left_join(design_subset, by = "AMPLICON_NAME")

  # Report merge summary to the user
  message(glue::glue("Original HDR barcoded file rows: {nrow(HDR_df)}."))
  message(glue::glue("Merged file rows: {nrow(merged_df)}."))
  message(glue::glue(
    "Columns in merged file: {paste(colnames(merged_df), collapse = ', ')}."
  ))

  return(merged_df)
}


#' @keywords internal
validate_amplicon_merge <- function(
  merged_df,
  HDR_df
) {
  message("\n-- VALIDATING MERGE --\n")

  # Calculation of validation statistics
  has_sequence <- !is.na(merged_df$`Amplicon sequence`)
  entries_with_sequences <- sum(has_sequence)
  entries_without_sequences <- sum(!has_sequence)
  success_rate <- (entries_with_sequences / nrow(merged_df)) * 100

  # Update user on validation statistics
  message(glue::glue(
    "All original HDR entries preserved: ",
    "{ifelse(nrow(merged_df) == nrow(HDR_df), '\u2713', '\u2717')}."
  ))
  message(glue::glue(
    "No duplicate rows introduced: ",
    "{ifelse(nrow(merged_df) == nrow(dplyr::distinct(merged_df)), '\u2713', '\u2717')}."
  ))
  message(glue::glue(
    "Entries with amplicon sequences: ",
    "{entries_with_sequences} out of {nrow(merged_df)}."
  ))
  message(glue::glue(
    "Entries without amplicon sequences: {entries_without_sequences}."
  ))
  message(glue::glue("Success rate: {sprintf('%.1f', success_rate)}%."))
}


#' @keywords internal
summarise_by_amplicon_name <- function(
  merged_df
) {
  message("\n-- SEQUENCE STATUS BY AMPLICON NAME --\n")

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
  message("\n-- SAMPLE OF MERGED DATA --\n")

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

      message(glue::glue("\nAMPLICON_NAME: {row$AMPLICON_NAME}."))
      message(glue::glue("sgRNA_SEQUENCE: {row$sgRNA_SEQUENCE}."))
      message(glue::glue("SNP_change: {row$SNP_change}."))

      # Limit printing of amplicon sequence to 50 bases
      if (nchar(amp_seq) > 50) {
        message(glue::glue(
          "Amplicon sequence (first 50 chars): ",
          "{substr(amp_seq, 1, 50)}...\n\n"
        ))
      } else {
        message(glue::glue("Amplicon sequence: {amp_seq}.\n\n"))
      }
    }
  }
}


#' @keywords internal
export_merged_amplicon_data <- function(
  merged_df,
  output_dir
) {
  message("-- EXPORTING MERGED DATA --\n")

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    write_csv(merged_df, file = file.path(
      output_dir, "HDR_file_barcoded_with_amplicon_sequences.csv"))

    # Report data export outcome to the user
    message(glue::glue(
      "\u2713 Merged data exported to: {normalizePath(output_dir)}."
    ))
    message(glue::glue("\u2713 Total rows: {nrow(merged_df)}."))
    message(glue::glue("\u2713 Total columns: {ncol(merged_df)}."))

  } else {
    message(
      "\u2717 No output directory specified - merged data could not be exported."
    )
  }
}


#' @keywords internal
summarise_amplicon_merge_outcomes <- function(
  merged_df
) {
  message("\n-- MERGED DATA STRUCTURE --\n")

  # Print each column name to the user
  for (i in seq_along(colnames(merged_df))) {
    message(glue::glue("{i}. {colnames(merged_df)[i]}."))
  }
}
