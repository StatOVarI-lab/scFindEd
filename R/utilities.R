#' @keywords internal
# Checks if both seq1 and seq2 are present in each element of vector v
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


#' @keywords internal
# Extracts up to 20 bases following seq_before_barcode_1
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


#' @keywords internal
# Extracts up to 30 bases following seq_before_barcode_2
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


#' @keywords internal
# Extracts up to 30 bases preceding seq_after_barcode_3
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


#' @keywords internal
#' Returns the n_barcodes whitelist barcodes with the smallest distance
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


#' @keywords internal
#' Corrects a barcode by selecting the closest match, as long as the closest
#' match has its distance below a defined threshold.
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


#' @keywords internal
#' Converts a list of barcode distances into a data frame, better suited for
#' building the downstream QC ggplot histogram
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


#' @keywords internal
#' Produces a histogram of Hamming distances to the nearest whitelist barcode.
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


#' @keywords internal
#' Takes a list of data frame names with vectors of column names to check for in
#' each, and returns TRUE/FALSE dependent on whether every required column is
#' present.
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


# Inputs: two character vectors of equal length (NAs allowed).
# Returns data.frame with merged_sequence (NA if no overlap) and overlap_length (0 if none).
#' @keywords internal
merge_r1_r2 <- function(
  r1,
  r2,
  min_overlap = 4L,
  verbose
) {
  if (verbose) cli::cli_process_start(msg = "Merging R1 and R2...",
                                      msg_done = "Sucessfully merged R1 and R2.")
  
  # Ensure both read vectors are the same length
  if (length(r1) != length(r2)) {
    cli::cli_abort("Reads r1 and r2 are not of the same length.")
  }
  
  # Convert both r1 and r2 to characters
  # Important in case they are passed as factors
  r1 <- as.character(r1)
  r2 <- as.character(r2)
  n <- length(r1)
  
  # Record lengths of every read for both r1 and r2
  len1 <- ifelse(is.na(r1), 0L, nchar(r1))
  len2 <- ifelse(is.na(r2), 0L, nchar(r2))
  
  # Calculates the minimum overlap possible between every read pair
  max_pair_overlap <- pmin(len1, len2)
  overlap <- integer(n)
  
  # Only consider pairs with potential overlap greater than the minimum overlap
  eligible <- which(max_pair_overlap >= min_overlap)
  
  # Only consider if there are eligible read pairs
  if (length(eligible) > 0L) {
    # Record the highest overlap across all read pairs
    maxLEN <- max(max_pair_overlap[eligible])
    
    if (verbose) cli::cli_progress_bar(
      "Calculating overlaps",
      total = length(seq.int(maxLEN, min_overlap, -1L))
    )
    
    # Iterate from the highest overlap to the minimum overlap parameter
    # in descending order
    for (LEN in seq.int(maxLEN, min_overlap, -1L)) {
      # Return indices of eligible pairs that haven't yet found an overlap,
      # and can support an overlap of length LEN
      cand <- eligible[overlap[eligible] == 0L & max_pair_overlap[eligible] >= LEN]
      
      # Skip to next value of LEN, if no candidate pairs exist
      if (length(cand) == 0L) next
      
      # Extract suffix and prefix of the read pairs
      s1_starts <- len1[cand] - LEN + 1L
      s1_sub <- substring(r1[cand], s1_starts, len1[cand])
      s2_sub <- substring(r2[cand], 1L, LEN)
      
      # Return indices of candidate pairs that match exactly across length LEN
      matches <- which(!is.na(s1_sub) & !is.na(s2_sub) & s1_sub == s2_sub)
      
      # If there are matches for length LEN, record their overlap lengths
      if (length(matches) > 0L) overlap[cand[matches]] <- as.integer(LEN)
      
      if (verbose) cli::cli_progress_update() # Update progress bar
      
      # If all eligible pairs have a found valid overlap, break function early
      if (all(overlap[eligible] > 0L)) break
    }
    
    if (verbose) cli::cli_progress_done() # Clear progress bar
  }
  
  # Create character vector of same length as the read vectors
  merged <- rep(NA_character_, n)
  
  # Record indices of read pairs that had overlap
  has_ol <- which(overlap > 0L)
  
  # Only consider if there are any read pairs that have overlap
  if (length(has_ol) > 0L) {
    # Determine start point for extracting the overlap tail from r2
    tail_start <- overlap[has_ol] + 1L
    
    # Extracts the non-overlapping suffix of r2 for each overlapping pair
    tails <- substring(r2[has_ol], tail_start, len2[has_ol])
    
    # Merges r1 with the non-overlapping suffix of r2 for each overlapping pair
    merged[has_ol] <- paste0(r1[has_ol], tails)
  }
  
  if (verbose) cli::cli_process_done()
  
  return(data.frame(
    merged_sequence = merged,
    overlap_length = overlap,
    stringsAsFactors = FALSE
  ))
}


#' @keywords internal
#' Searches for substring matches of a target string in DNA sequences.
#' All sequences and search strings are converted to uppercase for easier
#' comparison.
#' Returns: Logical vector indicating which sequences contain the target string
#' as a substring, or list with logical vector and indices if return_indices =
#' TRUE
check_sequences_contain_string <- function(
  data,
  target_string,
  column_name = "merged_sequence",
  batch_size = 10000,
  return_indices = FALSE,
  strict_dna = TRUE
) {
  # Input validation
  if (!column_name %in% colnames(data)) {
    cli::cli_abort("Column {columns_name} not found in data.")
  }
  
  if (is.null(target_string) ||
      length(target_string) == 0 ||
      nchar(target_string) == 0) {
    cli::cli_abort("Input target_string cannot be an empty data object.")
  }
  
  # Capitalise target_string
  target_string <- toupper(target_string)
  
  # Extract sequences and calculate total count
  sequences <- data[[column_name]]
  n_sequences <- length(sequences)
  
  # Calculate n_batches from the n_sequences and batch size
  n_batches <- ceiling(n_sequences / batch_size)
  
  # Initialise objects for storing match results
  matches <- logical(n_sequences)
  match_indices <- integer(0)
  
  # Process in batches
  for (batch_idx in seq_len(n_batches)) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, n_sequences)
    
    # Extract batch
    batch_sequences <- sequences[start_idx:end_idx]
    
    # Substring search using grepl
    batch_matches <- grepl(target_string, batch_sequences, fixed = TRUE)
    batch_matches[is.na(batch_matches)] <- FALSE
    
    # Store batch results
    matches[start_idx:end_idx] <- batch_matches
    
    # Collect indices if needed
    if (return_indices) {
      batch_match_indices <- which(batch_matches) + start_idx - 1
      match_indices <- c(match_indices, batch_match_indices)
    }
  }
  
  if (return_indices) {
    return(list(matches = matches, indices = match_indices))
  } else {
    return(matches)
  }
}




