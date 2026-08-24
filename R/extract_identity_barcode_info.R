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

  identity_amplicon <- reference_amplicons[identity_AmpID][[1]]

  # -- Extract N info --

  if (verbose) cli::cli_process_start(msg = "Extracting identity barcode info from {identity_AmpID}...",
                                      msg_done = "Identity barcode info from {identity_AmpID} extracted.")

  # Extract all 'N' bases
  individual_Ns <- Biostrings::matchPattern("N", identity_amplicon, fixed = TRUE)

  # Collapse adjacent 'N' positions into respective chunks
  N_chunks <- IRanges::reduce(IRanges::ranges(individual_Ns))

  identity_barcode_info <- data.frame(
    amp_ID = identity_AmpID,
    amplicon_seq = as.character(identity_amplicon),
    barcode_ID = seq(1, length(N_chunks)),
    start = IRanges::start(N_chunks),
    length = IRanges::width(N_chunks),
    end = IRanges::end(N_chunks)
  )

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    identity_barcode_info_file_name <- paste0(sample_name, "_identity_barcode_info.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(identity_barcode_info,
                     file = file.path(output_dir, identity_barcode_info_file_name))

    cli::cli_alert_success("Identity barcode information saved to: {normalizePath(output_dir)}.")
  }
}


