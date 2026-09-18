#' Extract Identity Barcode Information
#'
#' This function takes the `ref_amplicons_sorted.fa` file produced by scFindEd-nf
#' and extracts Identity barcode information from the sequence for the provided
#' `identity_AmpID`.
#'
#' @inheritParams global_params source_dir sample_name identity_AmpID output_dir plots_output_dir verbose
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
  plots_output_dir = NULL,
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
  identity_barcode_info[["barcode_ID"]] <- paste0("BC", seq(1, nrow(identity_barcode_info)))

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    identity_barcode_info_file_name <- paste0(sample_name, "_identity_barcode_info.csv")

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(identity_barcode_info,
                     file = file.path(output_dir, identity_barcode_info_file_name))

    cli::cli_alert_success("Identity barcode information saved to: {normalizePath(output_dir)}.")
  }

  # -- Identity Barcode Visualisation (only if plots_output_dir is specified) --
  if (!is.null(plots_output_dir)) {
    visualisation_plot <- identity_barcode_visualisation(identity_barcode_info)

    visualisation_plot_file_name <- paste0(sample_name, "_identity_barcode_visualisation.png")

    # Optional file saving (only if plots_output_dir is provided)
    if (!dir.exists(plots_output_dir)) dir.create(plots_output_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(plots_output_dir, visualisation_plot_file_name),
                    plot = visualisation_plot,
                    width = 8, height = 1.5 * length(identity_AmpID) + 0.5, dpi = 300)

    cli::cli_alert_success("Identity barcode visualisation exported to: {normalizePath(plots_output_dir)}")
  }
}


#' Produces a visualisation of the identity barcodes on the identity amplicons
#' @keywords internal
identity_barcode_visualisation <- function(identity_barcode_info) {
  # Add amplicon length, plus additional bar size details
  identity_barcode_info <- identity_barcode_info |>
    dplyr::mutate(amplicon_length = nchar(amplicon_seq),
                  ymin = -0.7, ymax = 0.7, mid = (start + end) / 2)

  # Lock amplicons in a stable, predictable order
  identity_AmpIDs <- unique(identity_barcode_info[["AmpID"]])
  identity_barcode_info[["AmpID"]] <- factor(identity_barcode_info[["AmpID"]],
                                             levels = identity_AmpIDs)

  # Extract a secondary table with unique rows for each amplicon
  amp_lengths <- identity_barcode_info |>
    dplyr::distinct(AmpID, amplicon_seq, amplicon_length)

  plot <- ggplot2::ggplot() +
    # backbone: full-length amplicon shaded grey; barcode boxes are drawn
    # on top of this, so only the non-barcode remainder shows as grey
    ggplot2::geom_rect(
      data = amp_lengths,
      ggplot2::aes(xmin = 1, xmax = amplicon_length, ymin = -0.7, ymax = 0.7),
      fill = "grey85", colour = "grey60", linewidth = 0.3
    ) +
    # barcode regions as coloured boxes
    ggplot2::geom_rect(
      data = identity_barcode_info,
      ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = barcode_ID),
      colour = "black", alpha = 0.85
    ) +
    # barcode labels, sitting above each box
    ggplot2::geom_text(
      data = identity_barcode_info,
      ggplot2::aes(x = mid, y = ymax + 0.3, label = barcode_ID),
      size = 3.5, fontface = "bold"
    ) +
    # position labels, sitting below each box
    ggplot2::geom_text(
      data = identity_barcode_info,
      ggplot2::aes(x = mid, y = ymin - 0.3, label = paste0(start, "-", end)),
      size = 3, colour = "grey20"
    ) +
    ggplot2::facet_wrap(~AmpID, ncol = 1, scales = "fixed") +
    ggplot2::scale_x_continuous(name = "Amplicon Position (bp)") +
    ggplot2::scale_y_continuous(limits = c(-1.2, 1.2), breaks = NULL, name = NULL) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Barcode") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      aspect.ratio       = 0.1,
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.spacing.y    = ggplot2::unit(1.2, "lines"),
      strip.text         = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position    = "none"
    ) +
    ggplot2::labs(title = "Identity Barcode Visualisation")

  return(plot)
}


