#' Identifying Frameshift Indels
#'
#' This function is an optional extension to `call_identity_barcodes()`.
#'
#' This function takes the output file from `call_identity_barcodes()` - the
#' identity barcode genotypes file - and appends an additional column that
#' states whether there is a frameshift indel or in-frame indel contained in the
#' consensus genotype. It can be used to identify identity barcodes that could
#' be used as labelled knockouts in downstream scRNA analysis.
#'
#' @inheritParams global_params source_dir sample_name ploidy output_dir verbose
#'
#' @return Appends an additional column to
#' "\{sample_name\}_\{ploidy\}_genotypes_summary.csv" that details whether the
#' genotype contains an `"In-Frame"` or `"Frameshift indel"`. If the identity
#' barcode's genotype contains no indel, returns `NA`. If the identity barcode
#' is genotyped on multiple amplicons, the output frameshift column will be
#' pipe-separated in the same order as the genotype column.
#'
#' @export
identify_frameshift_indels <- function(
  source_dir,
  sample_name,
  ploidy = "haploid",
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import identity barcode genotypes --

  identity_barcode_genotypes_file_name <- paste0(sample_name, "_", ploidy,
                                                 "_identity_barcode_genotypes.csv")
  identity_barcode_genotypes <- import_identity_barcode_genotypes(
    source_dir,
    identity_barcode_genotypes_file_name,
    verbose
  )

  if (verbose) cli::cli_process_start(msg = "Identifying and categorising indel mutations...",
                                      msg_done = "Indel mutations identified and categorised.")

  identity_barcode_genotypes <- identity_barcode_genotypes |>
    # Separate out genotypes with multiple genotypes on different amplicons
    tidyr::separate_longer_delim(genotype, delim = " | ") |>
    # Detect in each genotype's code, and track the AmpID of origin
    dplyr::mutate(AmpID = dplyr::case_when(
      stringi::stri_detect_fixed(genotype, ":")  ~ stringi::stri_extract_first_regex(genotype, "^[^:]*"),
      .default                                   = NA
    )) |>
    # If AmpID detected, remove the AmpID prefix
    dplyr::mutate(mutations = dplyr::case_when(
      !is.na(AmpID) ~ stringi::stri_replace_first_regex(genotype, "^.*?:", ""),
      .default      = NA
    )) |>
    # Now separate out each mutation on genotypes made up of multiple mutations
    tidyr::separate_longer_delim(mutations, delim = ";") |>
    # For each mutation, use custom function to calculate size of indel
    # NOTE: substitutions return NA
    dplyr::mutate(indel_length = dplyr::case_when(
      !is.na(mutations)  ~ calculate_indel_size(mutations),
      .default           = NA
    )) |>
    dplyr::summarise(
      # Summarise frameshift across all mutations on the same single AmpID genotype
      frameshift = dplyr::case_when(
        # If multiple indels on the same amplicon, do not attempt to classify
        sum(!is.na(indel_length))  > 1  ~ "Multiple Indels",
        # If the single indel's length is divisible by 3:     In-Frame
        sum(!is.na(indel_length)) == 1 & sum(indel_length, na.rm = TRUE) %% 3 == 0  ~ "In-Frame",
        # If the single indel's length is not divisible by 3: Frameshift
        sum(!is.na(indel_length)) == 1 & sum(indel_length, na.rm = TRUE) %% 3 != 0  ~ "Frameshift",
        # Not indel
        .default                        = NA_character_
      ), .by = c(identity_barcode, n_cells, genotype, contains_edit)
    ) |>
    # Re-collapse frameshift identification for barcodes with multiple genotypes
    dplyr::summarise(
      genotype   = paste(genotype,   collapse = " | "),
      frameshift = paste(frameshift, collapse = " | "),
      .by       = c(identity_barcode, n_cells, contains_edit)
    )

  if (verbose) cli::cli_process_done()

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    identity_barcode_genotypes_file_name <- paste0(sample_name, "_", ploidy,
                                                   "_identity_barcode_genotypes.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(identity_barcode_genotypes, file = file.path(
      output_dir, identity_barcode_genotypes_file_name))

    cli::cli_alert_success(
      "Frameshift-labelled identity barcode genotypes saved to: {normalizePath(output_dir)}.")
  }
}
