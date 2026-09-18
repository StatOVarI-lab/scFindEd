#' Identifying Closest Whitelisted Barcodes
#'
#' This function identifies the shortest Hamming distances between the extracted
#' barcodes and the whitelist of barcodes for each of the identity barcodes in
#' the experiment.
#'
#' @inheritParams global_params source_dir sample_name output_dir plots_output_dir verbose
#'
#' @param whitelist_barcode_files A named list pairing each identity barcode ID
#' with .csv file name from which to draw the whitelisted barcode from.
#' @param n_threads An integer value that confirms how many threads/cores the
#' function should utilise in constructing the barcode distance matrices.
#'
#' @return A .rds list of the shortest distance whitelist barcodes for each of
#' the identity barcodes in the experiment. Also produces a Hamming
#' distance-based quality control plot for each identity barcode, showing the
#' distributions of distances between the extracted barcodes and the first and
#' second closest whitelist barcodes.
#'
#' @export
identify_whitelisted_barcodes <- function(
  source_dir,
  sample_name,
  whitelist_barcode_files = NULL,
  n_threads = 1,
  output_dir = NULL,
  plots_output_dir = NULL,
  verbose = TRUE
) {
  # -- Import cell to extracted BCs --

  cell_to_extracted_BCs_file_name <- paste0(sample_name, "_cell_to_extracted_BCs.csv")
  cell_to_extracted_BCs <- import_cell_to_extracted_BCs(source_dir,
                                                        cell_to_extracted_BCs_file_name,
                                                        verbose)

  # Ensure each barcode in cell_to_extracted_BCs has an associated whitelist_barcode_file
  if (!is.list(whitelist_barcode_files)) {
    cli::cli_abort("whitelist_barcode_files must be a named list of whitelist barcode file names.")
  }
  extracted_BC_IDs = setdiff(colnames(cell_to_extracted_BCs), c("cell_barcode", "n"))
  whitelist_BC_IDs = names(whitelist_barcode_files)

  if (!identical(sort(extracted_BC_IDs), sort(whitelist_BC_IDs))) {
    cli::cli_abort('The key names for whitelist_barcode_files must exactly match the barcode IDs in "cell_to_extracted_BCs.csv".')
  }

  # -- Compare to whitelists --

  for (extracted_BC_ID in extracted_BC_IDs) {
    cli::cli_h2(extracted_BC_ID)

    # Import barcode whitelist
    BC_ID_whitelist_file_name <- whitelist_barcode_files[[extracted_BC_ID]]
    BC_ID_whitelist <- import_barcode_whitelist(source_dir,
                                                BC_ID_whitelist_file_name,
                                                verbose)
    names(BC_ID_whitelist) <- paste0(extracted_BC_ID, "_", 1:length(BC_ID_whitelist))

    # -- Extract closest whitelist barcodes --

    no_closest_barcodes <- 3
    BC_ID_closest_barcodes <- closest_whitelist_barcodes(cell_to_extracted_BCs,
                                                         extracted_BC_ID,
                                                         BC_ID_whitelist,
                                                         no_closest_barcodes,
                                                         n_threads,
                                                         verbose)

    # -- Optional file saving (only if output_dir is provided) --
    if (!is.null(output_dir)) {
      BC_ID_closest_barcodes_file_name <- paste0(sample_name, "_", extracted_BC_ID,
                                                 "_closest_barcodes.rds")

      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      saveRDS(BC_ID_closest_barcodes,
              file = file.path(output_dir, BC_ID_closest_barcodes_file_name))

      cli::cli_alert_success(
        "Closest whitelist barcodes exported to: {normalizePath(output_dir)}")
    }

    # -- Hamming distance QC plot --

    if (!is.null(plots_output_dir)) {
      BC_ID_hamming_plot <- barcode_hamming_plot(
        BC_ID_closest_barcodes,
        paste0(extracted_BC_ID, ": Hamming Distance Comparison"))

      BC_ID_hamming_plot_file_name <- paste0(sample_name, "_", extracted_BC_ID,
                                             "_hamming_plot.png")

      # Optional file saving (only if plots_output_dir is provided)
      if (!dir.exists(plots_output_dir)) dir.create(plots_output_dir, recursive = TRUE)
      ggplot2::ggsave(file.path(output_dir, BC_ID_hamming_plot_file_name),
                      plot = BC_ID_hamming_plot,
                      width = 8, height = 5, dpi = 300)

      cli::cli_alert_success("QC plot for {extracted_BC_ID} exported to: {normalizePath(output_dir)}")
    }
  }
}


#' @keywords internal
closest_whitelist_barcodes <- function(
  cell_to_extracted_BCs,
  extracted_BC_ID,
  barcode_whitelist,
  no_closest_barcodes,
  n_threads,
  verbose
) {
  # Get list of unique extracted barcodes
  unique_extracted_barcodes <- unique(cell_to_extracted_BCs[[extracted_BC_ID]])
  no_unique_extracted_barcodes <- length(unique_extracted_barcodes)

  # Get list of whitelist names
  whitelist_names <- names(barcode_whitelist)

  if (verbose) cli::cli_alert_info("{length(unique_extracted_barcodes)} unique extracted barcodes x {length(barcode_whitelist)} whitelist barcodes.")

  # Batch the distance matrix calculations to reduce memory overhead
  chunk_size <- 1000L
  BC_ID_closest_barcodes <- vector("list", no_unique_extracted_barcodes)

  cli::cli_progress_bar("Computing whitelist distances",
                        total = no_unique_extracted_barcodes)

  for (start in seq(1L, no_unique_extracted_barcodes, by = chunk_size)) {
    # Set start and end indices for batch
    end <- min(start + chunk_size - 1L, no_unique_extracted_barcodes)
    idx <- start:end

    # Calculate distance matrix between barcode whitelist and batch of extracted barcodes
    dist_batch <- stringdist::stringdistmatrix(
      unique_extracted_barcodes[idx],
      barcode_whitelist,
      method = "hamming",
      useBytes = TRUE,
      nthread = n_threads # Exploit multi-threading
    )

    # Extract top X distances per row using partial quick-select sorting
    BC_ID_closest_barcodes[idx] <- lapply(seq_len(nrow(dist_batch)), function(i) {
      row_dist <- dist_batch[i, ]
      top_idxs <- order(row_dist)[1:no_closest_barcodes] # Pull out the top ordered vals
      output <- row_dist[top_idxs]
      names(output) <- whitelist_names[top_idxs] # Re-attach whitelist barcode IDs
      output
    })

    cli::cli_progress_update(inc = (end - start + 1L))
  }

  cli::cli_alert_success("Whitelist barcode distances computed.")

  # Assign barcode sequences as names to the list entries
  names(BC_ID_closest_barcodes) <- unique_extracted_barcodes

  return(BC_ID_closest_barcodes)
}


#' Produces a histogram of Hamming distances to the nearest whitelist barcode.
#' @keywords internal
barcode_hamming_plot <- function(
  closest_barcodes,
  title
) {
  # Convert list to matrix once
  dist_mat <- do.call(rbind, closest_barcodes)

  # Build the long data frame directly
  hamming_df <- data.frame(
    distance = c(dist_mat[, 1], dist_mat[, 2]),
    type = factor(rep(c("Closest", "Second Closest"), each = nrow(dist_mat)))
  )

  plot <- ggplot2::ggplot(hamming_df, ggplot2::aes(x = distance, fill = type)) +
    ggplot2::geom_histogram(position = "dodge", binwidth = 0.5, color = "black", alpha = 0.8) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = title,
      subtitle = "Hamming distances between extracted and closest whitelist barcodes",
      x = "Hamming Distance",
      y = "Number of Cells",
      fill = "Barcode Match"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.y = ggplot2::element_line(color = "grey95", linewidth = 0.2)
    )

  return(plot)
}
