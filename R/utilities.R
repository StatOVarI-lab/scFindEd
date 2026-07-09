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


#' Takes a list of data frame names with vectors of column names to check for in
#' each, and returns TRUE/FALSE dependent on whether every required column is
#' present.
#' @keywords internal
check_required_columns <- function(
  required_cols
) {
  columns_present = TRUE

  for (df_name in names(required_cols)) {
    check_df <- get(df_name, envir = parent.frame())
    cols_to_check <- required_cols[[df_name]]

    # Finds columns present in required_cols but not in check_df
    missing_cols <- setdiff(cols_to_check, colnames(check_df))

    # Update the user if columns are missing
    if(length(missing_cols) > 0) {
      cli::cli_abort(
        "Missing column(s) in {df_name}: {paste(missing_cols, collapse = ', ')}")
      columns_present = FALSE
    }
  }
  return(columns_present)
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



