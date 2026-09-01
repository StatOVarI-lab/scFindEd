#' Extract Identity Barcode Information
#'
#' This function takes the `ref_amplicons_sorted.fa` file produced by scFindEd-nf
#' and extracts Identity barcode information from the sequence for the provided
#' `identity_AmpID`.
#'
#' @inheritParams global_params source_dir sample_name identity_AmpID output_dir verbose
#'
#' @return A table describing the start position, length and end position of each
#' barcode in the identity barcode amplicon.
#'
#' @export
extract_identity_barcode_info <- function(
  source_dir,
  sample_name,
  identity_AmpID,
  output_dir = NULL,
  verbose
) {
  # -- Import Reference Amplicons --

  reference_amplicons <- import_reference_amplicons(source_dir, verbose)

  if (!all(identity_AmpID %in% names(reference_amplicons))) {
    cli::cli_abort('All AmpIDs in identity_AmpID must be present in "ref_amplicons_sorted.fa".')
  }

  identity_amplicons <- reference_amplicons[identity_AmpID]

  # -- Set Up Output Dataframe --

  identity_barcode_info <- data.frame(AmpID = character(),
                                      amplicon_seq = character(),
                                      barcode_ID = integer(),
                                      start = integer(),
                                      length = integer(),
                                      end = integer(),
                                      stringsAsFactors = FALSE)

  # -- Extract N info --

  for (AmpID in identity_AmpID) {
    if (verbose) cli::cli_process_start(msg = "Extracting identity barcode info from {AmpID}...",
                                        msg_done = "Identity barcode info from {AmpID} extracted.")

    # Extract amplicon sequence for AmpID
    identity_amplicon <- identity_amplicons[[AmpID]]

    # Extract all 'N' bases
    individual_Ns <- Biostrings::matchPattern("N", identity_amplicon, fixed = TRUE)

    # Collapse adjacent 'N' positions into respective chunks
    N_chunks <- IRanges::reduce(IRanges::ranges(individual_Ns))

    # Populate into data frame and append to bottom of identity_barcode_info
    AmpID_barcode_info <- data.frame(
      AmpID = AmpID,
      amplicon_seq = as.character(identity_amplicon),
      barcode_ID = NA,
      start = IRanges::start(N_chunks),
      length = IRanges::width(N_chunks),
      end = IRanges::end(N_chunks)
    )

    identity_barcode_info <- dplyr::bind_rows(identity_barcode_info,
                                              AmpID_barcode_info)

    if (verbose) cli::cli_process_done()
  }

  # Assign unique barcode ID to each extracted barcode
  identity_barcode_info[["barcode_ID"]] <- seq(1, nrow(identity_barcode_info))

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    identity_barcode_info_file_name <- paste0(sample_name, "_identity_barcode_info.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(identity_barcode_info,
                     file = file.path(output_dir, identity_barcode_info_file_name))

    cli::cli_alert_success("Identity barcode information saved to: {normalizePath(output_dir)}.")
  }
}


