#' @keywords internal
#' Imports data from a file
#' Expected file name: "amplicons_barcoded_data.txt"
import_amplicons_barcoded_data <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- readr::read_tsv(file_path, show_col_types = FALSE),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}


#' @keywords internal
#' Imports data from a file
#' Expected file name: HDR_file_barcoded_with_amplicon_sequences.csv"
import_HDR_barcoded_data <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- suppressMessages(readr::read_csv(file_path, show_col_types = FALSE)),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}


#' @keywords internal
#' Imports data from a file
#' Expected file name: "sample_name_info.csv"
import_sample_info_data <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- data.table::fread(file_path, showProgress = verbose),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )
  data.table::setnames(data, c("readname", "amplicon_name", "start",
                               "quality", "cigar", "sequence"))

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}


#' @keywords internal
#' Imports data from a file
#' Expected file name: "sample_name_barcodes_fastq.txt"
import_barcode_reads_data <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- readLines(file_path),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}


#' @keywords internal
#' Imports data from a file
#' Expected file name: "sample_name_R2_001_readnames_fastq.txt"
import_readnames_data <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- data.table::fread(
      file_path, header = FALSE, showProgress = verbose)$V1,
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}


#' @keywords internal
#' Imports data from a file
#' Expected file name: "cell_to_UTR_barcode_system.csv"
import_cell_to_UTR_barcode_system <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- utils::read.csv(file_path, row.names = 1),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )
  data <- tibble::as_tibble(data)

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
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
  if (verbose) cli::cli_inform("Merging R1 and R2...")

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

  if (verbose) cli::cli_alert_success("Successfully merged R1 and R2.")

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












