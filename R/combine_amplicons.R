#' Combine Amplicons
#'
#' This function takes the combined amplicon summaries file and the cell to
#' UTR barcode system file, merging them and summarising cell genotypes across
#' amplicons.
#'
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param combined_summaries_file_name String containing combined amplicon
#' summaries file name.
#' @param barcode_system_file_name String containing the cell to UTR barcode
#' system file name.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return A data frame containing summarised cell genotypes across each
#' amplicon, optionally saved to a CSV file entitled "cell_genotypes_summary".

#' @examples
#' cell_genotypes <- combine_amplicons(
#' source_dir = "data",
#' combined_summaries_file_name = "final_summary_minimum_2counts_70_30_70_all_amplicons.csv",
#' barcode_system_file_name = "cell_to_UTR_barcode_system.csv",
#' output_dir = "data",
#' verbose = TRUE
#' )

#' @importFrom dplyr %>%

#' @export
combine_amplicons <- function(
  source_dir,
  combined_summaries_file_name,
  barcode_system_file_name,
  output_dir,
  verbose = TRUE
) {
  # -- Data Import --

  imported_data <- combine_amplicons_data_import(
    source_dir,
    combined_summaries_file_name,
    barcode_system_file_name,
    verbose
  )

  invisible(list2env(imported_data, envir = environment()))

  # -- Merge data --

  merged_df <- merge_amplicon_summaries_and_barcode_system(
    combined_summaries_df,
    barcodes_file_df,
    source_dir,
    combined_summaries_file_name,
    barcode_system_file_name,
    verbose
  )

  # -- Check for Consistent Genotype among Different Amplicons --

  amplicon_genotypes_filtered <- filter_inconsistent_genotypes(merged_df,
                                                               verbose)

  # -- Summarise Cell Genotypes

  cell_genotypes <- summarise_cell_genotypes(amplicon_genotypes_filtered,
                                             verbose)

  # Report Summarised Genotypes

  if (verbose) report_summarised_genotypes(cell_genotypes)

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    file_save_name <- file.path(output_dir, "cell_genotypes_summary.csv")

    readr::write_csv(cell_genotypes, file_save_name, col_names = TRUE)

    cli::cli_alert_success(
      "Summarised cell genotypes saved to: {normalizePath(output_dir)}")
  }

  return(cell_genotypes)
}


#' @keywords internal
combine_amplicons_data_import <- function(
  source_dir,
  combined_summaries_file_name,
  barcode_system_file_name,
  verbose
) {
  if (verbose) {
    cli::cli_h2("IMPORTING DATA")
  } else {
    cli::cli_inform("Importing datasets...")
  }

  # -- Import combined amplicon summaries data --

  combined_summaries_df <- import_combined_summaries_data(
    source_dir,
    combined_summaries_file_name,
    verbose
  )

  # -- Import cell to UTR barcode system

  barcodes_file_df <- import_cell_to_UTR_barcode_system(
    source_dir,
    barcode_system_file_name,
    verbose
  )

  barcodes_file_df <- barcodes_file_df %>%
    dplyr::rename(barcode = cell_barcode)

  # Select only the cell barcodes where UTR_BC_count >= 2
  barcodes_file_df <- barcodes_file_df %>%
    dplyr::add_count(all_BCs, name = "UTR_BC_count")
  barcodes_file_df <- barcodes_file_df[barcodes_file_df$UTR_BC_count >= 2,]

  if (verbose) cli::cli_alert_info(glue::glue(
    "Number of unique UTR_BCs in raw Cell to UTR barcode system: ",
    "{length(unique(barcodes_file_df$all_BCs))}."
  ))

  if (verbose) cli::cli_text("\n")
  cli::cli_alert_success("Successfully imported data.")

  return(list(combined_summaries_df = combined_summaries_df,
              barcodes_file_df = barcodes_file_df))
}


#' @keywords internal
merge_amplicon_summaries_and_barcode_system <- function(
  combined_summaries_df,
  barcodes_file_df,
  source_dir,
  combined_summaries_file_name,
  barcode_system_file_name,
  verbose
) {
  if (verbose) {
    cli::cli_h2("PERFORMING MERGE")

    cli::cli_alert_info(glue::glue(
      "Merging {combined_summaries_file_name} and ",
      "{barcode_system_file_name} by `barcode`..."
    ))
  }

  # Merge and keep only the overlapping barcodes
  tryCatch({
    merged_df <- merge(combined_summaries_df, barcodes_file_df, by = "barcode")
    },
    error = function(e) {
      cli::cli_abort(
        "Failed to merge {.file combined_summaries_file_name} and
        {.file barcode_system_file_name} by {.val 'barcode'}.")
      }
    )

  cli::cli_alert_success("Datasets successfully merged.")

  if (verbose) {
    cli::cli_alert_info(
      "Number of unique cell barcodes: {length(unique(merged_df$barcode))}.")
    cli::cli_alert_info(
      "Number of unique UTR barcodes: {length(unique(merged_df$all_BCs))}.")
  }

  return(merged_df)
}


#' @keywords internal
filter_inconsistent_genotypes <- function(
  merged_df,
  verbose
) {
  if (verbose) cli::cli_h2("CHECKING FOR CONSISTENT GENOTYPES")

  # Check if there is consistent genotype among different amplicons
  amplicon_genotypes <- merged_df %>%
    dplyr::group_by(all_BCs, amplicon) %>%
    dplyr::summarise(
      n_WT        = sum(genotype == "WT"),
      n_HET       = sum(genotype == "HET"),
      n_HOM       = sum(genotype == "HOM"),
      n_others    = sum(genotype == "others"),
      # Number of rows - number of genotypes within the same amplicon (and all BCs)
      n_genotypes = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      genotype = dplyr::case_when(
        n_WT     == n_genotypes ~ "WT",
        n_HET    == n_genotypes ~ "HET",
        n_HOM    == n_genotypes ~ "HOM",
        n_others == n_genotypes ~ "others",
        TRUE ~ "mixed"
      )
    )

  # Remove the inconsistent amplicon genotypes within the same UTR_BCs
  amplicon_genotypes_filtered <- amplicon_genotypes[amplicon_genotypes$genotype!= "mixed",]

  if (verbose) {
    cli::cli_alert_success("Inconsistent amplicon genotypes removed.")

    cli::cli_alert_info("Dimensions: {dim(amplicon_genotypes_filtered)}.")
    print(utils::head(amplicon_genotypes_filtered))
  }

  return(amplicon_genotypes_filtered)
}


#' @keywords internal
summarise_cell_genotypes <- function(
  amplicon_genotypes_filtered,
  verbose
) {
  if (verbose) cli::cli_h2("SUMMARISING CELL GENOTYPES")

  # Check the genotype for each cell (UTR_BC)
  cell_genotypes <- amplicon_genotypes_filtered %>%
    dplyr::group_by(all_BCs) %>%
    dplyr::summarise(
      # Amplicon name for the non-WT genotype (if any)
      amplicon = paste(unique(amplicon[genotype != "WT"]), collapse = ","),
      n_amplicons = dplyr::n(),
      n_WT = sum(genotype == "WT"),
      n_HET = sum(genotype == "HET"),
      n_HOM = sum(genotype == "HOM"),
      n_others = sum(genotype == "others"),
      n_mutations = n_amplicons - n_WT,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      genotype = dplyr::case_when(
        n_WT == n_amplicons ~ "WT",
        n_WT == n_amplicons - 1 & n_HET == 1 ~ "HET",
        n_WT == n_amplicons - 1 & n_HOM == 1 ~ "HOM",
        n_WT == n_amplicons - 1 & n_others == 1 ~ "others",
        TRUE ~ "mixed"
      ),
      # Keep amplicon info only for non-WT genotypes
      amplicon = dplyr::if_else(
        genotype %in% c("HET", "HOM", "others"), amplicon, NA_character_)
    )

  if (verbose) {
    cli::cli_alert_success("Cell genotypes summarised.")

    cli::cli_alert_info("Dimensions: {dim(cell_genotypes)}.")
    cli::cli_alert_info(
      "Number of unique UTR barcodes: {length(unique(cell_genotypes$all_BCs))}."
    )
  }

return(cell_genotypes)
}


#' @keywords internal
report_summarised_genotypes <- function(
  cell_genotypes
) {
  cli::cli_h2("SUMMARY STATISTICS")

  cli::cli_alert_info("Summary of genotype counts:")
  print(table(cell_genotypes$genotype))

  cli::cli_alert_info("Summary of n_mutation counts:")
  print(table(cell_genotypes$n_mutations))

  cli::cli_alert_info("Number of records with n_mutations == 0 or 1: ")
  print(nrow(cell_genotypes[cell_genotypes$n_mutations == 0 |
                              cell_genotypes$n_mutations == 1,]))

  cli::cli_alert_info("Count of HOMOZYGOUS genotypes for each amplicon:")
  print(table(cell_genotypes[cell_genotypes$genotype =="HOM",]$amplicon))

  cli::cli_alert_info("Count of HETEROZYGOUS genotypes for each amplicon:")
  print(table(cell_genotypes[cell_genotypes$genotype =="HET",]$amplicon))

  cli::cli_alert_info("Count of 'OTHERS' genotypes for each amplicon:")
  print(table(cell_genotypes[cell_genotypes$genotype =="others",]$amplicon))
}





