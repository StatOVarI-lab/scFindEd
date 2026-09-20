#' Corrects a barcode by selecting the closest match, as long as the closest
#' match has its distance below a defined threshold.
#' @keywords internal
correct_BC <- function(distances, threshold){
  # If the distance < the threshold, select the barcode as a match
  if (distances[1] < threshold){
    output <- names(distances)[1]
  }
  # Otherwise, return no output
  else {
    output <- ""
  }
  return(output)
}


#' Cumulative sum that resets at the start of each contiguous group, and allows
#' data.table to vectorise cleanly
#' @keywords internal
vectorised_cumsum <- function(ref_step, row_ID) {
  cs <- cumsum(as.numeric(ref_step))
  row_ID_start <- c(TRUE, row_ID[-1] != row_ID[-length(row_ID)])
  row_ID_sizes <- rle(row_ID)$lengths
  offset <- rep(cs[row_ID_start] - ref_step[row_ID_start], row_ID_sizes)

  return(cs - offset)
}


#' Chunk CIGAR strings and flatten them into a single, long data table
#' @keywords internal
construct_cigar_chunks_dt <- function(cigars, row_IDs, batch_size = 1000000) {
  # Construct batches
  n <- length(cigars)
  batch_starts <- seq.int(1, n, by = batch_size)
  batch_ends <- pmin(batch_starts + batch_size - 1, n)

  batched_dt <- lapply(seq_along(batch_starts), function(i) {
    # Get batch indices
    batch_idx <- batch_starts[i]:batch_ends[i]

    # Break the CIGAR strings into their constituent chunks
    cigar_ops <- cigarillo::cigars_as_RleList(cigars[batch_idx])
    # Extract chunk operations and chunk lengths separately
    ops  <- S4Vectors::runValue(cigar_ops)
    lens <- S4Vectors::runLength(cigar_ops)

    # Flatten all CIGAR chunks into a long data table
    data.table::data.table(
      row_ID = rep(row_IDs[batch_idx], lengths(ops)),
      op     = unlist(ops, use.names = FALSE),
      len    = unlist(lens, use.names = FALSE)
    )
  })

  return(data.table::rbindlist(batched_dt, use.names = FALSE))
}


#' Calculate a read's end position from the start position and the read's CIGAR string
#' @keywords internal
calculate_read_end_position <- function(pos, cigars, batch_size = 1000000) {
  # Construct batches
  n <- length(cigars)
  batch_starts <- seq.int(1, n, by = batch_size)
  batch_ends <- pmin(batch_starts + batch_size - 1, n)

  cigar_extents <- unlist(lapply(seq_along(batch_starts), function(i) {
    # Get batch indices
    batch_idx <- batch_starts[i]:batch_ends[i]

    # Calculate distance CIGAR extends along the reference
    cigarillo::cigar_extent_along_ref(cigars[batch_idx])
  }), use.names = FALSE)

  return(pos + cigar_extents - 1)
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

