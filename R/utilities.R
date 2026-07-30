#' Checks if both seq1 and seq2 are present in each element of vector v
#' @keywords internal
grepl_seq1_2 <- function(
  v,
  seq1,
  seq2
) {
  # Check for presence of seq1
  temp_1 <- sapply(v, function(x) grepl(seq1, x, fixed = TRUE))

  # Check for presence of seq2
  temp_2 <- sapply(v, function(x) grepl(seq2, x, fixed = TRUE))

  # Return TRUE when both temp_1 and temp_2 are TRUE
  return(temp_1 & temp_2)
}


#' Extracts up to 20 bases following seq_before_barcode_1
#' @keywords internal
extract_BC1 <- function(
  x,
  seq_before_barcode_1
) {
  # Find seq_before_barcode_1 within x
  barcode_check <- regexpr(seq_before_barcode_1, x, fixed = TRUE)

  # Report error if seq_before_barcode_1 is not present inside of x
  if (barcode_check == -1) {
    cli::cli_abort("extract_BC1(): flanking sequence not found.
                   This likely indicates an upstream filtering issue.")
  }

  barcode_start <- barcode_check + attr(barcode_check, "match.length")

  # Extracts up to 20 bases following seq_before_barcode_1
  substr(x, barcode_start, barcode_start + 19)
}


#' Extracts up to 30 bases following seq_before_barcode_2
#' @keywords internal
extract_BC2 <- function(
  x,
  seq_before_barcode_2
) {
  # Find seq_before_barcode_2 within x
  barcode_check <- regexpr(seq_before_barcode_2, x, fixed = TRUE)

  # Report error if seq_before_barcode_2 is not present inside of x
  if (barcode_check == -1) {
    cli::cli_abort("extract_BC2(): flanking sequence not found.
                   This likely indicates an upstream filtering issue.")
  }

  barcode_start <- barcode_check + attr(barcode_check, "match.length")

  # Extracts up to 30 bases following seq_before_barcode_2
  substr(x, barcode_start, barcode_start + 29)
}


#' Extracts up to 30 bases preceding seq_after_barcode_3
#' @keywords internal
extract_BC3 <- function(
  x,
  seq_after_barcode_3
) {
  # Find seq_after_barcode_3 within x
  barcode_check <- regexpr(seq_after_barcode_3, x, fixed = TRUE)

  # Report error if seq_after_barcode_3 is not present inside of x
  if (barcode_check == -1) {
    cli::cli_abort("extract_BC3(): flanking sequence not found.
                   This likely indicates an upstream filtering issue.")
  }

  # Get position of the barcode start from the position of the barcode end
  barcode_end <- barcode_check - 1
  barcode_start <- max(1, barcode_end - 29)

  # Extracts up to 30 bases preceding seq_after_barcode_3
  substr(x, barcode_start, barcode_end)
}


#' Returns the n_barcodes whitelist barcodes with the smallest distance
#' @keywords internal
extract_closest_barcodes <- function(
  i,
  dist_matrix,
  whitelist_names,
  n_barcodes
) {
  distances <- dist_matrix[i, ]

  # Attach the barcode names to the distances
  names(distances) <- whitelist_names

  # Return the lowest 3 distances
  return(sort(distances)[1:3])
}


#' Corrects a barcode by selecting the closest match, as long as the closest
#' match has its distance below a defined threshold.
#' @keywords internal
correct_BC <- function(distances, thresh){
  # If the distance < the threshold, select the barcode as a match
  if (distances[1] < thresh){
    output <- names(distances)[1]
  }
  # Otherwise, return no output
  else {
    output <- ""
  }
  return(output)
}


#' Converts a list of barcode distances into a data frame, better suited for
#' building the downstream QC ggplot histogram
#' @keywords internal
prepare_hamming_plot_dataframe <- function(
  BC_list,
  BC_name,
  verbose
) {
  # Distance from nearest whitelist barcode
  closest <- as.numeric(unlist(lapply(BC_list, function(x) x[1])))

  # Distance from second nearest whitelist barcode
  second_closest <- as.numeric(unlist(lapply(BC_list, function(x) x[2])))

  # Combine into a data frame
  BC_df <- data.frame(
    distance = c(closest, second_closest),
    type = factor(c(rep("Closest", length(closest)),
                    rep("Second Closest", length(second_closest))))
  )

  return(BC_df)
}


#' Produces a histogram of Hamming distances to the nearest whitelist barcode.
#' @keywords internal
barcode_hamming_plot <- function(
  data,
  title
) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = distance, fill = type)) +
    ggplot2::geom_histogram(position = "dodge", binwidth = 0.5, color = "black", alpha = 0.8) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = title,
      subtitle = "Distribution of distances by barcode match type",
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


#' Calculate a read's end position from the start position and the read's CIGAR string
#' @keywords internal
calculate_read_end_position <- function(pos, cigar) {
  # Break the cigar strings into their constituent chunks
  cigar_ops <- cigarillo::cigars_as_RleList(cigar)

  # Extract chunk operations and chunk lengths separately
  ops  <- S4Vectors::runValue(cigar_ops)
  lens <- S4Vectors::runLength(cigar_ops)

  # Create a mask for only ops with value "M" (match) or "D" (deletion)
  M_or_D <- ops == "M" | ops == "D"
  lens_filtered <- lens * M_or_D

  # Sum the lengths of these "M" and "D" chunks
  total_alignment_lengths <- sum(lens_filtered)

  return (pos + total_alignment_lengths - 1)
}


#' Calculates a coordinate window around the bases that contain mutations in
#' the variants in compiled_alignments.
#' @keywords internal
calculate_informative_window <- function(
  compiled_alignments,
  reference_seq,
  design_summary_file,
  amplicon
) {
  # Extract basic reference coordinate data
  all_mutation_positions <- compiled_alignments |>
    dplyr::filter(is_variant == TRUE) |>
    dplyr::pull(sequence_pos)
  # In order to define the informative window
  if (length(all_mutation_positions) > 0) {
    window_start <- max(1, min(all_mutation_positions) - 5)
    window_end   <- min(nchar(reference_seq), max(all_mutation_positions) + 5)
  } else {
    edit_positions <- design_summary_file |>
      dplyr::filter(AmpID == amplicon) |>
      dplyr::select(start, length)
    window_start <- min(edit_positions$start) - 5
    window_end   <- max(edit_positions$start + edit_positions$length) + 4
  }

  return(list(
    window_start = window_start,
    window_end = window_end
  ))
}


#' Used to convert a variant_code to an aligned track of bases relative to a
#' reference sequence.
#' @keywords internal
mutation_code_alignment_tracking <- function(
  variant_code,
  reference_seq
) {
  # Initialise tracking vectors based on the reference sequence
  seq_bases       <- strsplit(reference_seq, "")[[1]]
  is_base_variant <- rep(FALSE, length(seq_bases))
  is_base_sub     <- rep(FALSE, length(seq_bases))
  base_event_ID   <- rep(NA_integer_, length(seq_bases))
  base_event_type <- rep(NA_character_, length(seq_bases))

  # Process variants if the sequence is not wildtype
  if (variant_code != "WT") {
    mutations <- strsplit(variant_code, ";")[[1]]

    # Sort mutations right-to-left (descending) to safeguard index alignments
    mutation_positions <- as.numeric(stringr::str_extract(mutations, "\\d+"))
    mutations          <- mutations[order(mutation_positions, decreasing = TRUE)]

    event_counter <- 1
    for (mutation in mutations) {
      clean_mutation <- sub("^c\\.", "", mutation)

      # -- Substitions --
      if (grepl(">", clean_mutation)) {
        pos <- as.numeric(stringr::str_extract(clean_mutation, "^\\d+"))

        seq_bases[pos]       <- stringr::str_extract(clean_mutation, "[A-Z]$")
        is_base_sub[pos]     <- TRUE
        is_base_variant[pos] <- TRUE
        base_event_ID[pos]   <- event_counter
        base_event_type[pos] <- "substitution"

        # -- Deletions --
      } else if (grepl("del", clean_mutation)) {
        range_vals <- as.numeric(stringr::str_extract_all(clean_mutation, "\\d+")[[1]])
        idx_range  <- range_vals[1]:tail(range_vals, 1)

        seq_bases[idx_range]       <- "-"
        is_base_variant[idx_range] <- TRUE
        base_event_ID[idx_range]   <- event_counter
        base_event_type[idx_range] <- "deletion"

        # -- Insertions --
      } else if (grepl("ins", clean_mutation)) {
        pos          <- as.numeric(stringr::str_extract(clean_mutation, "^\\d+"))
        insert_bases <- strsplit(stringr::str_extract(clean_mutation, "[A-Z]+$"), "")[[1]]
        ins_len      <- length(insert_bases)

        # Capture tracking vector lengths safely before expansion
        current_len <- length(seq_bases)
        idx_left    <- 1:pos
        idx_right   <- (pos + 1):current_len

        # Expand tracking arrays using precise boundary pieces
        seq_bases       <- c(seq_bases[idx_left],       insert_bases,                seq_bases[idx_right])
        is_base_variant <- c(is_base_variant[idx_left], rep(TRUE, ins_len),          is_base_variant[idx_right])
        base_event_ID   <- c(base_event_ID[idx_left],   rep(event_counter, ins_len), base_event_ID[idx_right])
        base_event_type <- c(base_event_type[idx_left], rep("insertion", ins_len),   base_event_type[idx_right])
        is_base_sub     <- c(is_base_sub[idx_left],     rep(FALSE, ins_len),         is_base_sub[idx_right])
      }

      event_counter <- event_counter + 1
    }
  }

  # 3. Unified return statement
  return(data.frame(
    base       = seq_bases,
    is_variant = is_base_variant,
    is_sub     = is_base_sub,
    event_ID   = base_event_ID,
    event_type = base_event_type,
    stringsAsFactors = FALSE
  ))
}

#' Used to grab pre-defined colours for each base in producing plots.
#' @keywords internal
get_base_colours <- function() {
  return (c("A" = "#adf7b6",
            "T" = "#ffc09f",
            "C" = "#ffee93",
            "G" = "#a0ced9",
            "-" = "#f6f6f6"))
}

#' Used to grab pre-defined colours for each mutation type when constructing
#' borders in producing variant summary plots.
#' @keywords internal
get_border_colours <- function() {
  return(c("insertion"    = "#1B5E20",
           "deletion"     = "#B71C1C",
           "substitution" = "#0D47A1"))
}

