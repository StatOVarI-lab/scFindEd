#' Plot Variant Summary
#'
#' This function is used to visualise the distribution of variants across a
#' sample. It can be used to generate plot(s) for all variants in each amplicon,
#' or in a specific amplicon. It can also be used to generate plot(s) for all
#' variants that belong to a specific cell barcode / UTR barcode combination,
#' also for either every amplicon or a specific amplicon.
#'
#' @param amplicon A string containing the amplicon ID (AmpID), to be used when
#' generating a plot for a specific amplicon.
#' @param barcode A string containing the cell barcode / UTR barcode combination,
#' to be used when generating plot(s) only for a specific cell / UTR barcode
#' combination subset.
#'
#' @inheritParams global_params source_dir sample_name panel_ID output_dir verbose
#'
#' @return By default, exports a bulked variant summary plot for every amplicon.
#' Can also export a specific amplicon, a specific cell barcode / UTR barcode
#' combination, or both.
#'
#' @export
plot_variant_summary <- function(
  source_dir,
  sample_name,
  panel_ID,
  amplicon = NULL,
  barcode = NULL,
  output_dir = NULL,
  verbose = TRUE
) {
  # Import SNP summary
  design_summary_file_name <- paste0(panel_ID, "_design_summary_joined.csv")
  design_summary_file <- import_design_summary_info(source_dir,
                                                    design_summary_file_name,
                                                    verbose)

  # If amplicon IS defined, filter for only amplicon
  if (!is.null(amplicon)) {
    if (amplicon %in% design_summary_file$AmpID) {
      design_summary_file <- dplyr::filter(design_summary_file, AmpID == amplicon)
    } else {
      cli::cli_abort("{amplicon} is not included in {design_summary_file_name}.")
    }
  }

  # Extract unique amplicon IDs
  amplicon_list <- unique(design_summary_file$AmpID)
  n_amplicons <- length(amplicon_list)

  # Extract unique reference amplicon sequences
  reference_seqs <- unique(design_summary_file$amplicon_seq)

  # -- Repeat for each amplicon... --

  for (i in 1:n_amplicons) {
    AmpID <- amplicon_list[i]
    reference_seq <- reference_seqs[i]

    if (verbose) cli::cli_h2(AmpID)

    # Import true mutations
    true_mutations_file_name <- paste0(sample_name, "_", AmpID,
                                       "_true_mutations.tsv")
    true_mutations <- import_true_mutations(source_dir,
                                            true_mutations_file_name,
                                            verbose)

    # Check for presence of cell_barcode / UTR_barcode
    if (!is.null(barcode)) {
      if (!barcode %in% true_mutations$cell_barcode &
          !barcode %in% true_mutations$UTR_barcode) {
        # Skip to next amplicon if not present
        cli::cli_alert_warning("{barcode} is not present in {AmpID}.")
        next
      }
    }

    # -- Select Variants --

    variant_data <- select_variants(true_mutations, AmpID, barcode, verbose)

    # -- Structure Alignment Levels --

    # Define isolated reference labels and the spacer lane
    standalone_ref <- "WT_Reference"
    gap_separator <- "GAP_Spacer"

    # Construct the y-axis hierarchy
    y_axis_levels <- c(standalone_ref, gap_separator, variant_data$variant)

    # Create a clean lookup vector to format the y-axis labels later
    axis_label_lookup <- c(
      "Reference Amplicon",
      "",
      variant_data$display_label
    )
    names(axis_label_lookup) <- y_axis_levels

    # -- Compile Alignments --

    compiled_alignments <- compile_alignments(variant_data,
                                              reference_seq,
                                              y_axis_levels,
                                              axis_label_lookup,
                                              verbose)

    # -- Calculate Informative Window --

    informative_window_coordinates <- calculate_informative_window(compiled_alignments,
                                                                   reference_seq,
                                                                   design_summary_file,
                                                                   AmpID)

    compiled_alignments <- compiled_alignments |>
      dplyr::filter(sequence_pos >= informative_window_coordinates[["window_start"]] &
                    sequence_pos <= informative_window_coordinates[["window_end"]])

    # -- Construct Variant Summary Plot --

    variant_summary_plot <- construct_variant_summary_plot(compiled_alignments,
                                                           AmpID,
                                                           barcode,
                                                           y_axis_levels,
                                                           axis_label_lookup,
                                                           informative_window_coordinates,
                                                           verbose)

    # -- Optional file saving (only if output_dir is specified) --

    if (!is.null(output_dir)) {
      variant_summary_plot_file_name <- paste0(AmpID, "_variant_summary.png")

      # Create sub-directory only if barcode is defined
      if (!is.null(barcode)) {
        if (!dir.exists(paste0(output_dir, "/", barcode))) {
          dir.create(paste0(output_dir, "/", barcode), recursive = TRUE)
        }
        variant_summary_plot_file_name <- paste0(
          barcode, "/", variant_summary_plot_file_name)
      }

      # Dynamically adjust plot height based on the number of variants
      plot_height <- 1.5 + (nrow(variant_data) / 5)
      if ((nrow(dplyr::filter(compiled_alignments, is_variant == TRUE))) != 0) {
        plot_height <- plot_height + 0.5
      }

      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      ggplot2::ggsave(
        file.path(output_dir, variant_summary_plot_file_name),
        plot = variant_summary_plot,
        width = 10,
        height = plot_height,
        dpi = 300)
    }
  }

  cli::cli_alert_success(
    "Variant summary plot(s) saved to: {normalizePath(output_dir)}.")
}


#' @keywords internal
select_variants <- function(
  true_mutations,
  AmpID,
  barcode,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Filtering and selecting variants...",
                                      msg_done = "Variants filtered and selected.")

  if (is.null(barcode)) {
    # Remove wildtypes, only if not filtering by barcode
    true_mutations <- true_mutations |>
      dplyr::filter(combined_HGVSc_denoised != "WT")
  }

  # If barcode IS defined, filter for only cell_barcode / UTR_barcode == barcode
  if (!is.null(barcode)) {
    if (barcode %in% true_mutations$cell_barcode) {
      true_mutations <- dplyr::filter(true_mutations, cell_barcode == barcode)
    } else if (barcode %in% true_mutations$UTR_barcode) {
      true_mutations <- dplyr::filter(true_mutations, UTR_barcode == barcode)
    }
  }

  # Calculate counts and proportions of each variant, as well as a total count
  true_mutations <- true_mutations |>
    dplyr::group_by(sample_id, amp_id, combined_HGVSc_denoised) |>
    dplyr::summarise(variant_copy = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::ungroup() |>
    dplyr::mutate(total_count = sum(variant_copy)) |>
    dplyr::mutate(variant_proportion = variant_copy / total_count) |>
    dplyr::arrange(desc(variant_proportion), desc(variant_copy))

  # Limit the number of printed variants to 20 rows
  variant_data <- true_mutations |>
    dplyr::slice_head(n = 20) |>
    dplyr::select(
      variant    = combined_HGVSc_denoised,
      count      = variant_copy,
      proportion = variant_proportion
    ) |>
    dplyr::mutate(
      display_label = paste0(round(proportion * 100, 1), "% [", count, " reads]"))

  if (verbose) cli::cli_process_done()

  return(variant_data)
}


#' @keywords internal
compile_alignments <- function(
  variant_data,
  reference_seq,
  y_axis_levels,
  axis_label_lookup,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Structuring variant alignments...",
                                      msg_done = "Variant alignments structured.")

  # Build a comprehensive tracking map including the isolated reference lane
  mutation_mapping <- c("WT", variant_data$variant)
  names(mutation_mapping) <- c(y_axis_levels[1], variant_data$variant)

  # -- Convert Mutation Codes to Plottable Data --

  compiled_alignments <- lapply(names(mutation_mapping), function(label) {
    variant_code <- mutation_mapping[label]
    alignments_df <- mutation_code_alignment_tracking(variant_code, reference_seq)
    alignments_df$sequence_pos <- 1:nrow(alignments_df)
    alignments_df$sample <- label
    return(alignments_df)
  }) %>%
    dplyr::bind_rows()

  # Assign layout factor order explicitly
  compiled_alignments$sample <- factor(compiled_alignments$sample, levels = rev(y_axis_levels))

  if (verbose) cli::cli_process_done()

  return(compiled_alignments)
}


#' @keywords internal
construct_variant_summary_plot <- function(
  compiled_alignments,
  AmpID,
  barcode,
  y_axis_levels,
  axis_label_lookup,
  informative_window_coordinates,
  verbose
) {
  # -- Construct plot --

  if (verbose) cli::cli_process_start(msg = "Constructing variant summary plot...",
                                      msg_done = "Variant summary plot constructed.")

  # Variably define plot title
  plot_title <- AmpID
  if (!is.null(barcode)) {
    plot_title <- paste0(plot_title, " - ", barcode)
  }

  # Grab colour palettes
  base_colours <- get_base_colours()
  border_colours <- get_border_colours()

  # Isolate rows with mutations
  mutation_boxes <- compiled_alignments |>
    dplyr::filter(is_variant == TRUE)
  n_variants <- nrow(mutation_boxes)

  # Calculate edit bounding box coordinates
  if (n_variants > 0) {
    mutation_boxes <- mutation_boxes |>
      dplyr::group_by(sample, event_ID) |>
      dplyr::summarise(
        xmin = min(sequence_pos) - 0.5,
        xmax = max(sequence_pos) + 0.5,
        ymin = first(as.numeric(factor(sample, levels = rev(y_axis_levels)))) - 0.5,
        ymax = first(as.numeric(factor(sample, levels = rev(y_axis_levels)))) + 0.5,
        event_type = first(event_type),
        .groups = 'drop'
      )
  } else {
    # Empty fallback to prevent Inf errors when there are no variants
    mutation_boxes <- data.frame(
      xmin = numeric(0), xmax = numeric(0),
      ymin = numeric(0), ymax = numeric(0),
      event_type = character(0), sample = character(0)
    )
  }

  # Construct ggplot
  variant_summary_plot <- ggplot(compiled_alignments, aes(x = sequence_pos, y = sample)) +
    geom_tile(aes(fill = base), colour = "white", lwd = 0.2, show.legend = FALSE) +

    geom_text(aes(label = base), colour = "black", size = 2.6, fontface = "plain") +

    scale_y_discrete(
      drop = FALSE,
      limits = rev(y_axis_levels),
      labels = axis_label_lookup) +

    scale_fill_manual(values = base_colours) +

    scale_x_continuous(
      breaks = seq(informative_window_coordinates[["window_start"]],
                   informative_window_coordinates[["window_end"]],
                   by = 5),
      limits = c(informative_window_coordinates[["window_start"]] - 0.5,
                 informative_window_coordinates[["window_end"]] + 0.5),
      expand = c(0, 0)) +

    labs(title = plot_title,
         x = "Sequence Position Relative to Reference Amplicon",
         y = NULL)

  if (n_variants > 0) {
    # Only plot mutation boxes if there *are* mutations
    variant_summary_plot <- variant_summary_plot +
      geom_rect(data = mutation_boxes, aes(xmin = xmin, xmax = xmax,
                                           ymin = ymin, ymax = ymax,
                                           colour = event_type),
                fill = NA, lwd = 0.5, inherit.aes = FALSE) +

      scale_colour_manual(
        values = border_colours,
        breaks = c("insertion", "deletion", "substitution"),
        labels = c("Insertion", "Deletion", "Substitution"))
  }

  # Add theme for the plot *after* conditionally adding mutation boxes
  variant_summary_plot <- variant_summary_plot +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      axis.text.y = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

  if (verbose) cli::cli_process_done()

  return(variant_summary_plot)
}


