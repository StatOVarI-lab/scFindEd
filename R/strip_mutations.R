#' Strip Specific Mutations from Variant Summary
#'
#' This function is used to trim out certain mutations from each variant,
#' specifically those that were not removed during the function
#' `identify_true_mutations()`. The most probable reason for doing so would be
#' to remove a common, benign variant from each variant it occurs in, in order
#' to reconsider the reads as wildtype at that base position. This can make it
#' easier to call consistent genotypes downstream.
#'
#' The user is required to provide a specific `AmpID` from which variants will
#' be stripped from, as well as a single character or character vector of
#' `mutations` that should be stripped.
#'
#' @param AmpID String that designates which amplicons true_mutations file the
#' mutations should be stripped from.
#' @param mutations Character or character vector that contains each mutation
#' that should be stripped from each variant.
#'
#' @inheritParams global_params source_dir sample_name output_dir verbose
#'
#' @return Re-exports the true mutations file, containing the stripped HGVSc
#' mutation codes that uniquely identifies each variant and its count, for each
#' cell barcode.
#'
#' @export
strip_mutations <- function(
  source_dir,
  sample_name,
  AmpID,
  mutations,
  output_dir = NULL,
  verbose = TRUE
) {

  # Import true mutations
  true_mutations_file_name <- paste0(sample_name, "_", AmpID, "_true_mutations.tsv")
  true_mutations <- import_true_mutations(source_dir,
                                          true_mutations_file_name,
                                          verbose)


  if (verbose) cli::cli_process_start(msg = "Stripping designated mutations from {AmpID}...",
                                      msg_done = "Designated mutations stripped from {AmpID}.")

  true_mutations <- true_mutations |>
    dplyr::mutate(row_id = dplyr::row_number())

  # Keep only the mutations that are NOT in the noise list, re-collapsed per row
  true_mutations_stripped <- true_mutations |>
    dplyr::mutate(individual_mutation = strsplit(as.character(merged_HGVSc), ";")) |>
    tidyr::unnest(individual_mutation) |>
    dplyr::filter(!is.na(individual_mutation),
                  !individual_mutation %in% mutations) |>
    dplyr::group_by(row_id) |>
    dplyr::summarise(
      merged_HGVSc_stripped = paste(unique(individual_mutation), collapse = ";"),
      .groups = "drop"
    )

  true_mutations_output <- true_mutations |>
    dplyr::select(-merged_HGVSc) |>
    dplyr::left_join(true_mutations_stripped, by = "row_id") |>
    dplyr::mutate(
      # No mutations left after stripping -> WT
      merged_HGVSc_stripped = dplyr::if_else(
        is.na(merged_HGVSc_stripped) | merged_HGVSc_stripped == "",
        "WT",
        merged_HGVSc_stripped
      )
    ) |>
    dplyr::rename(merged_HGVSc = merged_HGVSc_stripped) |>
    dplyr::group_by(dplyr::across(dplyr::any_of(c(
      "sample_id", "cell_barcode", "amp_id", "identity_barcode", "cell_count", "merged_HGVSc")))) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::across(dplyr::any_of(c("identity_barcode", "cell_barcode"))),
                   dplyr::desc(count))

  if (verbose) cli::cli_process_done()

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    true_mutations_output_file_name <- paste0(sample_name, "_", AmpID,
                                              "_true_mutations.tsv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_tsv(true_mutations_output, file = file.path(
      output_dir, true_mutations_output_file_name))

    cli::cli_alert_success(
      "{AmpID} variants with stripped mutations saved to: {normalizePath(output_dir)}.")
  }
}

