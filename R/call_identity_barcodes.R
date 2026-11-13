#' Calling Each Identity Barcode's Genotype
#'
#' This function is used to collate the information on the most probable haploid
#' genotype across cells that share the same identity barcode, and provide a
#' consensus call for each identity barcode's overall genotype.
#'
#' Each mutation that exists within all of the constituent cells' most probable
#' genotypes are considered independently for whether they meet the threshold
#' for being assigned as part of the overall genotype.
#'
#' The `consistency_threshold` argument determines how many cells need to contain
#' a mutation in its most probable genotype for the mutation to be included in
#' the identity barcode's genotype. For example, at the default value of `0.5`,
#' explicitly more than 50% of cells must contain the mutation.
#'
#' Currently, the `lower_CI_assignment_threshold` parameter is used exclusively
#' for identity barcodes with only one cell. These barcodes require more
#' confidence to be called, and as such the cell itself must have a lower
#' credible interval bound greater than or equal to
#' `lower_CI_assignment_threshold`.
#'
#' An identity barcode is a "Confident WT" if every cell is itself a confident
#' wildtype. An identity barcode is a "Non-Confident WT" if instead more than
#' 50% of cells are a confident wildtype. If less than 50% of cells are confident
#' wildtype, and no other mutation is called confidently, the identity barcode
#' will be assigned as "Non-Confident Mutant".
#'
#' Lastly, an additional check is performed against the experiment's design
#' summary file. All assigned genotypes are checked against the intended edits
#' laid out in the design summary; if a genotype matches an edit on the
#' respective amplicon, then the column `contains_edit` is marked `TRUE`.
#'
#' @param consistency_threshold Float that sets the threshold for the minimum
#' number of cells in an identity barcode that must agree for a mutation to be
#' accepted. Must be in the range `[0, 1]`.
#'
#' @inheritParams global_params source_dir sample_name panel_ID ploidy lower_CI_assignment_threshold output_dir verbose
#'
#' @return Exports an identity barcode genotypes file, which contains a row for
#' each unique identity barcode, and describes the number of cells assigned to
#' the barcode, the barcode's consensus genotype, and whether the genotype
#' matches one of the experiment's intended edits.
#'
#' @export
call_identity_barcodes <- function(
  source_dir,
  sample_name,
  panel_ID,
  ploidy = "haploid",
  lower_CI_assignment_threshold = 0.7,
  consistency_threshold = 0.5,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Ensure all provided arguments are valid --
  validate_call_identity_barcodes_args(source_dir,
                                       sample_name,
                                       panel_ID,
                                       ploidy,
                                       lower_CI_assignment_threshold,
                                       consistency_threshold,
                                       output_dir,
                                       verbose)

  # -- Import SNP Summary --

  design_summary_file_name <- paste0(panel_ID, "_design_summary_joined.csv")
  design_summary_file <- import_design_summary_info(source_dir,
                                                    design_summary_file_name,
                                                    verbose)

  # -- Import genotype summary --

  genotype_summary_file_name <- paste0(sample_name, "_haploid_genotypes_summary.csv")
  genotypes_summary <- import_haploid_genotypes_summary(source_dir,
                                                        genotype_summary_file_name,
                                                        verbose)

  if (verbose) cli::cli_process_start(msg = "Summarising mutations & edits...",
                                      msg_done = "Mutations & edits summarised.")


  # Extract a list of edit genotypes
  edit_list <- design_summary_to_genotypes(design_summary_file)

  # Summarise wildtype
  WT_summaries <- genotypes_summary |>
    dplyr::group_by(identity_barcode) |>
    dplyr::summarise(
      n_cells        = dplyr::n(),
      n_WT_cells     = sum(is.na(all_mutations)),
      n_confident_WT = sum(confident_WT, na.rm = TRUE),
      .groups = "drop"
    )

  # One row per (cell, mutation); cells with no mutations drop out here
  mutations_long <- genotypes_summary |>
    dplyr::select(cell_barcode, identity_barcode, all_mutations) |>
    dplyr::filter(!is.na(all_mutations)) |>
    tidyr::separate_longer_delim(all_mutations, delim = " | ") |>
    dplyr::rename(mutation = all_mutations) |>
    dplyr::mutate(AmpID = sub(":.*", "", mutation))

  # Vote on each mutation separately
  mutations_consistency <- mutations_long |>
    dplyr::count(identity_barcode, AmpID, mutation, name = "n_carrying") |>
    dplyr::left_join(WT_summaries, by = "identity_barcode") |>
    dplyr::mutate(
      proportion = n_carrying / n_cells,
      is_edit    = mutation %in% edit_list
    ) |>
    dplyr::filter(proportion > consistency_threshold)

  if (verbose) cli::cli_process_done()

  # -- Slightly different protocol for identity barcodes with only one cell --

  one_cell_only_mutations <- mutations_consistency |>
    dplyr::filter(n_cells == 1)

  if (nrow(one_cell_only_mutations) >= 1) {
    if (verbose) cli::cli_alert_info("Single cell identity barcodes detected. Genotyping stringently...")

    # Import the complete set of haploid genotype summaries in order to look-
    # up CI lower values
    amplicon_list <- unique(design_summary_file[["AmpID"]])
    AmpID_genotype_summaries <- list()
    for (AmpID in amplicon_list) {
      AmpID_genotype_summary_file_name <- paste0(sample_name, "_", AmpID,
                                                 "_haploid_genotypes.tsv")
      AmpID_genotype_summaries[[AmpID]] <- import_AmpID_genotypes_summary(
        source_dir,
        AmpID_genotype_summary_file_name,
        verbose
      )
    }
    combined_genotype_summaries <- dplyr::bind_rows(AmpID_genotype_summaries)

    one_cell_only_mutations <- one_cell_only_mutations |>
      dplyr::mutate(variant = stringi::stri_replace_first_regex(mutation, "^.*?:", "")) |>
      dplyr::left_join(combined_genotype_summaries |> dplyr::select(
        identity_barcode, amp_id, variant_genotype, CI_lower
      ), by = c(
        "identity_barcode" = "identity_barcode",
        "AmpID" = "amp_id",
        "variant" = "variant_genotype"
      ))

    mutations_consistency <- mutations_consistency |>
      dplyr::left_join(one_cell_only_mutations |> dplyr::select(
        identity_barcode, AmpID, mutation, CI_lower
      ), by = c(
        "identity_barcode" = "identity_barcode",
        "AmpID" = "AmpID",
        "mutation" = "mutation"
      )) |>
      dplyr::filter(n_cells >= 2 | (n_cells == 1 & CI_lower >= lower_CI_assignment_threshold))
  }

  if (verbose) cli::cli_process_start(msg = "Calling identity barcode genotypes...",
                                      msg_done = "Identity barcode genotypes called.")

  # Combine passing mutations per barcode, starting from ALL barcodes
  barcode_calls <- WT_summaries |>
    dplyr::left_join(
      mutations_consistency |>
        dplyr::group_by(identity_barcode) |>
        dplyr::summarise(
          contains_edit = sum(is_edit) >= 1,
          mutations = paste(sort(mutation), collapse = " | "),
          .groups = "drop"
        ),
      by = "identity_barcode"
    )

  # Assign genotypes
  barcode_calls <- barcode_calls |>
    dplyr::mutate(
      genotype = dplyr::case_when(
        !is.na(mutations)                ~ mutations,
        n_confident_WT == n_cells        ~ "Confident WT",
        n_confident_WT / n_cells >= 0.5  ~ "Non-Confident WT",
        n_confident_WT / n_cells <  0.5  ~ "Non-Confident Mutant",
        .default                         = NA
      )
    ) |>
    dplyr::select(identity_barcode, n_cells, genotype, contains_edit)

  if (verbose) {
    cli::cli_process_done()

    total_identity_barcodes <- nrow(barcode_calls)
    n_confident_WTs <- sum(barcode_calls[["genotype"]] == "Confident WT")
    n_mutants <- sum(barcode_calls[["genotype"]] != "Confident WT" &
                       barcode_calls[["genotype"]] != "Non-Confident WT" &
                       barcode_calls[["genotype"]] != "Non-Confident Mutant" &
                       !is.na(barcode_calls[["genotype"]]))
    n_contains_edit <- sum(barcode_calls[["contains_edit"]], na.rm = TRUE)

    cli::cli_alert_info("Mutated identity barcodes: {n_mutants} / {total_identity_barcodes}")
    cli::cli_alert_info("Edted identity barcodes: {n_contains_edit} / {total_identity_barcodes}")
    cli::cli_alert_info("Confident WT identity barcodes: {n_confident_WTs} / {total_identity_barcodes}")
  }

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    called_identity_barcodes_file_name <- paste0(sample_name, "_", ploidy,
                                                 "_identity_barcode_genotypes.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(barcode_calls, file = file.path(
      output_dir, called_identity_barcodes_file_name))

    cli::cli_alert_success(
      "Called identity barcode genotypes saved to: {normalizePath(output_dir)}.")
  }
}



