#' Join Amplicon Summaries
#'
#' This function takes the amplicon summary files for each amplicon/HDR
#' (produced by process_reads), and joins them. It then performs minimum read
#' thresholding before categorising barcodes into
#' homozygous HDR/homozygous WT/heterozygous/other.
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param file_name String containing the file name structure for each amplicon's
#' summary file.
#' @param n_amplicons Integer denoting the number of amplicon summaries to be
#' joined.
#' @param minimum_count_thresholds Vector containing a list of mininum count
#' thresholds to filter the amplicon summaries at.
#' @param HOM_cutoff Integer denoting the minimum read percentage to be
#' classified homozygous - either WT or for the amplicon.
#' @param HET_low_cutoff Integer denoting the lower threshold read percentage of
#' both WT and the amplicon in order to be classified heterozygous.
#' @param HET_high_cutoff Integer denoting the upper threshold read percentage
#' of both WT and the amplicon in order to be classified heterozygous.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return Writes a CSV file containing the final joined summary for all
#' amplicons at the provided genotype thresholds, with one CSV file exported for
#' each value provided to `minimum_count_thresholds`.

#' @examples
#' join_amplicon_summaries(
#' source_dir = "data",
#' file_name = "reads_classified_MB_barcoded_HDR_",
#' n_amplicons = 9,
#' minimum_count_thresholds = c(2, 5, 10, 20),
#' HOM_cutoff = 70,
#' HET_low_cutoff = 30,
#' HET_high_cutoff = 70,
#' output_dir = "data",
#' verbose = TRUE
#' )

#' @export
join_amplicon_summaries <- function(
  source_dir,
  file_name,
  n_amplicons,
  minimum_count_thresholds,
  HOM_cutoff,
  HET_low_cutoff,
  HET_high_cutoff,
  output_dir,
  verbose = TRUE
) {
  # -- Load Datasets --

  amplicon_summary_dfs <- import_amplicon_summaries(
    source_dir,
    file_name,
    n_amplicons,
    verbose
  )

  if (verbose) {
    cat("\n")
    cli::cli_inform("Calculating summaries across Minimum Count Thresholds...")
  }

  for (count_threshold in minimum_count_thresholds){
    if (verbose & count_threshold != minimum_count_thresholds[1]) cat("\n")
    if (verbose) cli::cli_alert_info(
      "Calculating summaries for Minimum Count Threshold: {count_threshold}..."
    )

    # Initialise list to store results for each amplicon
    all_summaries <- list()

    for (amplicon_number in 1:n_amplicons) {

      # -- Summarise Genotypes --

      all_summaries[[amplicon_number]] <- summarise_genotypes(
        amplicon_summary_dfs[[amplicon_number]],
        amplicon_number,
        count_threshold,
        HOM_cutoff,
        HET_low_cutoff,
        HET_high_cutoff
      )
    }

    # Combine all amplicons into a single dataframe
    combined_summaries <- dplyr::bind_rows(all_summaries)

    # -- Optional file saving (only if output_dir is specified) --
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

      file_save_name <- file.path(
        output_dir, paste0("final_summary_minimum_",
                           count_threshold, "counts_",
                           HOM_cutoff, "_",
                           HET_low_cutoff, "_",
                           HET_high_cutoff, "_",
                           "all_amplicons.csv"
                           ))

      readr::write_csv(combined_summaries, file_save_name)

      cli::cli_alert_success(glue::glue(
        "Results for Minimum Count Threshold: {count_threshold} ",
        "saved to: {normalizePath(output_dir)}."))
    }
  }
}


#' @keywords internal
import_amplicon_summaries <- function(
  source_dir,
  file_name,
  n_amplicons,
  verbose
) {
  if (verbose) cli::cli_inform("Importing amplicon summary files...")

  amplicon_summary_dfs <- list()

  for (amplicon_number in 1:n_amplicons) {
    # Import dataset
    file_path <- file.path(source_dir, paste0(file_name, amplicon_number, ".csv"))
    tryCatch(
      df <- as.data.frame(readr::read_csv(file_path, show_col_types = FALSE)),
      error = function(e) {
        cli::cli_abort("Failed to import data from {.file {file_path}}.")
      }
    )

    # Standardise column names
    colnames(df) <- c("barcode", "HDR", "n", "pct")

    # Standardise HDR names in the HDR column
    df$HDR <- ifelse(df$HDR == paste0("HDR_", amplicon_number), "HDR", df$HDR)

    amplicon_summary_dfs[[amplicon_number]] <- df

    if (verbose) cli::cli_alert_success(
      "Successfully imported amplicon {amplicon_number} summary file.")
  }

  if (verbose) cat("\n")
  cli::cli_alert_success("Successfully imported all amplicon summary files.")

  return(amplicon_summary_dfs)
}


#' @keywords internal
summarise_genotypes <- function(
  amplicon_summary_df,
  amplicon_number,
  count_threshold,
  HOM_cutoff,
  HET_low_cutoff,
  HET_high_cutoff
) {
  # Compute barcode percentages --
  df_pct <- amplicon_summary_df %>%
    dplyr::group_by(barcode) %>%
    dplyr::summarise(
      WT_pct = sum(pct[HDR == "WT"], na.rm = TRUE),
      HDR_pct = sum(pct[HDR == "HDR"], na.rm = TRUE),
      total_count = sum(n),
      .groups = "drop"
    )

  # Assign genotypes
  df_pct_genotype <- df_pct %>%
    dplyr::mutate(
      genotype = dplyr::case_when(
        WT_pct >= HOM_cutoff ~ "WT",
        HDR_pct >= HOM_cutoff ~ "HOM",
        WT_pct > HET_low_cutoff & WT_pct < HET_high_cutoff &
          HDR_pct > HET_low_cutoff & HDR_pct < HET_high_cutoff ~ "HET",
        TRUE ~ "others"
      )
    )

  # Filter by threshold
  df_filtered <- df_pct_genotype %>%
    dplyr::filter(total_count >= count_threshold)

  # Add amplicon number
  df_filtered <- df_filtered %>%
    dplyr::mutate(amplicon = paste0("amplicon", amplicon_number))

  return(df_filtered)
}


