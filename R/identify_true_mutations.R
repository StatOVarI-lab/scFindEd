#' Identifying True Mutations in Read Variants
#'
#' This function is used to trim the 'noisy' mutations from each variant,
#' returning variant codes that contain only the likely 'true' mutations. A
#' 'noisy' mutation is one that appears in less than `true_mutation_proportion`
#' of the reads in a cell.
#'
#' This function also filters variants for a minimum number of reads per cell
#' - `min_reads`, and optionally for a minimum number of cells per identity
#' barcode combination - `min_cells`.
#'
#' @param min_reads Minimum number of reads for each cell. Cells with fewer than
#' `min_reads` reads will be filtered out.
#' @param true_mutation_proportion Minimum proportion of reads that a mutation
#' must appear in to be considered a 'true' mutation within it's cell.
#' @param min_cells Minimum number of cells per identity barcode combination.
#' Identity barcode combinations with fewer than `min_cells` cells will be
#' filtered out. Required if `barcoded = TRUE`.
#'
#' @inheritParams global_params source_dir sample_name panel_ID barcoded output_dir verbose
#'
#' @return For each amplicon, exports the true mutations file, containing the
#' stripped HGVSc mutation codes that uniquely identifies each variant and its
#' count, for each cell barcode.
#'
#' @export
identify_true_mutations <- function(
  source_dir,
  sample_name,
  panel_ID,
  min_reads = 10,
  true_mutation_proportion = 0.1,
  barcoded = FALSE,
  min_cells = NULL,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import SNP Summary --

  design_summary_file_name <- paste0(panel_ID, "_design_summary_joined.csv")
  design_summary_file <- import_design_summary_info(source_dir,
                                                    design_summary_file_name,
                                                    verbose)

  # Extract unique amplicon IDs
  amplicon_list <- unique(design_summary_file$AmpID)
  n_amplicons <- length(amplicon_list)

  # Create a list for storing complete set of retained cell barcodes
  retained_cell_barcodes <- vector("list", n_amplicons)

  for (i in 1:n_amplicons) {
    AmpID <- amplicon_list[i]

    # -- Import variant summaries --

    if (verbose) cli::cli_h2(AmpID)

    variant_summary_file_name <- paste0(sample_name, "_", AmpID,
                                        "_variant_summary.tsv")
    variant_summary <- import_variant_summary(source_dir,
                                              variant_summary_file_name,
                                              verbose)

    # -- Filter out discordant reads and count reads per cell --

    if (verbose) cli::cli_process_start(msg = "Removing discordant reads...",
                                        msg_done = "Discordant reads removed.")

    variant_summary <- variant_summary |>
      dplyr::filter(merged_HGVSc != "Discordance in overlapping region") |>
      dplyr::group_by(cell_barcode) |>
      dplyr::mutate(cell_count = sum(count, na.rm = TRUE)) |>
      dplyr::ungroup()

    if (verbose) {
      cli::cli_process_done()
      cli::cli_alert_info("Unique cell barcodes: {length(unique(variant_summary$cell_barcode))}")
    }

    # -- Filter by min_reads --

    if (verbose) cli::cli_process_start(msg = "Removing cells with fewer than {min_reads} reads...",
                                        msg_done = "Cells with fewer than {min_reads} removed.")

    variant_summary <- variant_summary |>
      dplyr::filter(cell_count >= min_reads) |>
      dplyr::mutate(row_id = dplyr::row_number())

    if (verbose) {
      cli::cli_process_done()
      cli::cli_alert_info("Unique cell barcodes: {length(unique(variant_summary$cell_barcode))}")
    }

    # -- Filter by min_cells (optional) --

    if (barcoded) {
      if (is.null(min_cells)) cli::cli_abort("In barcoded experiments, you must assign a value to min_cells.")
      if (verbose) cli::cli_process_start(msg = "Removing identity barcode combinations with fewer than {min_cells} cells...",
                                          msg_done = "Identity barcode combinations with fewer than {min_cells} removed.")

      variant_summary <- variant_summary |>
        dplyr::group_by(identity_barcode) |>
        dplyr::filter(dplyr::n_distinct(cell_barcode) >= min_cells) |>
        dplyr::ungroup()

      if (verbose) {
        cli::cli_process_done()
        cli::cli_alert_info("Unique cell barcodes: {length(unique(variant_summary$cell_barcode))}")
      }
    }

    # -- Identify 'true' mutations --

    if (verbose) cli::cli_process_start(msg = "Identifying 'true' mutations...",
                                        msg_done = "'True' mutations identified.")

    true_mutations <- variant_summary |>
      dplyr::mutate(individual_mutation = strsplit(as.character(merged_HGVSc), ";")) |>
      tidyr::unnest(individual_mutation) |>
      dplyr::group_by(cell_barcode, cell_count, individual_mutation) |>
      dplyr::summarise(mutation_count = sum(count, na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(mutation_count / cell_count > true_mutation_proportion) |>
      dplyr::arrange(cell_barcode, dplyr::desc(mutation_count))

    if (verbose) cli::cli_process_done()

    # -- Remove 'noise' mutations --

    if (verbose) cli::cli_process_start(msg = "Removing 'noise' mutations...",
                                        msg_done = "'Noise' mutations removed.")

    AMPL_variants_denoised <- variant_summary |>
      dplyr::mutate(individual_mutation = strsplit(as.character(merged_HGVSc), ";")) |>
      tidyr::unnest(individual_mutation) |>
      dplyr::semi_join(
        true_mutations,
        by = c("cell_barcode", "individual_mutation")
      ) |>
      dplyr::group_by(row_id) |>
      dplyr::summarise(
        merged_HGVSc_denoised = paste(unique(individual_mutation), collapse = ";"),
        .groups = "drop"
      )

    true_mutations <- variant_summary |>
      dplyr::select(-merged_HGVSc) |>  # Drop old column to avoid duplication
      dplyr::left_join(AMPL_variants_denoised, by = "row_id") |>
      dplyr::mutate(
        # If it's NA (meaning no mutations survived the filter), change it to "WT"
        merged_HGVSc_denoised = dplyr::if_else(
          is.na(merged_HGVSc_denoised) | merged_HGVSc_denoised == "",
          "WT",
          merged_HGVSc_denoised
        )
      ) |>
      dplyr::rename(merged_HGVSc = merged_HGVSc_denoised) |>
      dplyr::group_by(dplyr::across(dplyr::any_of(c(
        "sample_id", "cell_barcode", "amp_id", "identity_barcode", "cell_count", "merged_HGVSc")))) |>
      dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::across(dplyr::any_of(c("identity_barcode", "cell_barcode"))), dplyr::desc(count))

    if (verbose) cli::cli_process_done()

    # Note the retained cell barcodes for amplicon i
    retained_cell_barcodes[[i]] <- unique(true_mutations$cell_barcode)

    # -- Optional file saving (only if output_dir is specified) --

    if (!is.null(output_dir)) {
      true_mutations_output_file_name <- paste0(sample_name, "_", AmpID,
                                                "_true_mutations.tsv")

      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      readr::write_tsv(true_mutations, file = file.path(
        output_dir, true_mutations_output_file_name))

      cli::cli_alert_success(
        "{AmpID} true mutation variants saved to: {normalizePath(output_dir)}.")
    }
  }

  # Report total retained cell barcodes
  if (verbose) {
    cli::cli_h2("Summary")
    # Unique cell barcodes retained in at least one amplicon
    total_cells <- length(unique(unlist(retained_cell_barcodes)))
    cli::cli_alert_info("Total unique cell barcodes across all amplicons: {total_cells}.")
    # Unique cell barcodes retained in all amplicons
    shared_cells <- length(Reduce(intersect, retained_cell_barcodes))
    cli::cli_alert_info("Total unique cell barcodes across all amplicons: {shared_cells}.")
  }
}


