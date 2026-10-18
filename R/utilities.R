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


#' Used to convert design summary table to a list of all edits' HGVSc-style
#' mutation codes
#' @keywords internal
design_summary_to_genotypes <- function(design_summary) {
  AmpIDs <- unique(design_summary[["AmpID"]])

  genotypes <- character()

  for (ampID in AmpIDs) {
    AmpID_edits <- design_summary |>
      dplyr::select(AmpID, SNP, block_index, start, length, ref, change) |>
      dplyr::filter(AmpID == ampID) |>
      dplyr::arrange(length, start)  # controls the order genotypes appear in

    for (edit in unique(AmpID_edits$SNP)) {
      blocks <- AmpID_edits[AmpID_edits$SNP == edit, ]
      blocks <- blocks[order(blocks$block_index), ]

      edit_string <- character()
      for (i in seq_len(nrow(blocks))) {
        ref_bases <- strsplit(blocks$ref[i], "")[[1]]
        alt_bases <- strsplit(blocks$change[i], "")[[1]]
        pos       <- blocks$start[i] + seq_along(ref_bases) - 1

        edit_string <- c(edit_string, paste0("c.", pos, ref_bases, ">", alt_bases))
      }

      genotypes <- c(genotypes, paste0(ampID, ":", paste(edit_string, collapse = ";")))
    }
  }

  return(genotypes)
}


# -- Synthetic Experiment Generation --

# Transition matrix modeling human genomic sequential nucleotide frequencies
# (CpG depleted)
#' @keywords internal
markov_DNA_generation <- function(sequence_length = 250) {
  # Rows: From (A, C, G, T), Cols: To (A, C, G, T)
  P <- matrix(c(
    0.30, 0.20, 0.20, 0.30,  # A -> A, C, G, T
    0.32, 0.28, 0.08, 0.32,  # C -> A, C, G, T (CpG depleted)
    0.26, 0.22, 0.26, 0.26,  # G -> A, C, G, T
    0.22, 0.18, 0.28, 0.32   # T -> A, C, G, T
  ), nrow = 4, byrow = TRUE, dimnames = list(BASES, BASES))

  amplicon_seq <- character(sequence_length)
  amplicon_seq[1] <- sample(BASES, 1)

  for (i in seq_len(sequence_length)[-1]) {
    previous_base <- amplicon_seq[i-1]
    amplicon_seq[i] <- sample(BASES, 1, prob = P[previous_base, ])
  }

  return(paste0(amplicon_seq, collapse = ""))
}


#' Used to reverse complement a vector of sequences
#' @keywords internal
reverse_conversion <- function(seqs) {
  # Ensure conversion to character
  seqs <- as.character(seqs)

  # Set hashes for base complements
  old <- "ACGTacgt"
  new <- "TGCAtgca"

  # Due to high number of duplicate seqs, only perform conversion on unique sequences
  unique_seqs <- unique(seqs)
  unique_rev_comp <- stringi::stri_reverse(chartr(old, new, unique_seqs))
  reverse_complement_seqs <- unique_rev_comp[match(seqs, unique_seqs)]

  return(reverse_complement_seqs)
}


#' Used to insert mutations into a vector of sequences at a rate of 1 in 1000 bases (by default)
#' mutation_rate = per-base
#' p_sub = substitution proportion; p_ins = insertion proportion; p_del = deletion proportion
#' @keywords internal
mutate_seqs <- function(
  seqs,
  mutation_rate = 0.001,
  p_sub = 0.8,
  p_ins = 0.1,
  p_del = 0.1,
  indel_p = 0.5 # indel length ~ 1 + Geometric(indel_p); mean 1/indel_p, 1 = always single-base
) {
  # Convert to biostrings DNA set
  DNA_seqs     <- if (is.character(seqs)) Biostrings::DNAStringSet(seqs) else seqs
  n_seqs       <- length(DNA_seqs)
  seq_lengths  <- as.numeric(IRanges::width(DNA_seqs))
  total_length <- sum(seq_lengths)

  # Generate the total number of mutations based on total bases and mutation rate
  total_mutations <- rbinom(1L, total_length, mutation_rate)
  if (total_mutations == 0L) return(seqs)

  # Sample positions to mutate at, and determine the offsets
  pos     <- sort(sample.int(total_length, size = total_mutations))
  offsets <- c(0, cumsum(seq_lengths)[-n_seqs])
  seq_ID  <- findInterval(pos - 1, offsets)
  local   <- as.integer(pos - offsets[seq_ID])

  # Sample mutation type based on respective probabilities
  type  <- sample(c("sub", "ins", "del"), total_mutations, replace = TRUE,
                  prob = c(p_sub, p_ins, p_del))

  # -- Substitutions --

  # Substitutions must differ from the original base
  # Uses integer conversions and matching for efficiency
  orig <- as.character(Biostrings::subseq(DNA_seqs[seq_ID], start = local, width = 1L))
  idx  <- match(orig, BASES); idx[is.na(idx)] <- 0L
  alt  <- BASES[(idx - 1L + sample.int(3L, total_mutations, replace = TRUE)) %% 4L + 1L]

  # -- Indel Lengths --

  indel_len <- 1L + rgeom(total_mutations, indel_p)
  is_ins <- type == "ins"
  is_del <- type == "del"

  # Set base replacement values for substitutions (alt) and deletions (nothing)
  value <- alt
  value[is_del] <- ""

  # -- Insertions --

  # For insertions, draw a string of the sampled length (vectorised)
  if (any(is_ins)) {
    value[is_ins] <- vapply(indel_len[is_ins], markov_DNA_generation, character(1))
  }

  # -- Indel Range Widths --

  # ins = 0, sub = 1, del = sampled length (clamped so it can't run off the end)
  mutation_widths <- rep.int(1L, total_mutations)
  mutation_widths[is_ins] <- 0L
  mutation_widths[is_del] <- as.integer(pmin(indel_len[is_del],
                                             seq_lengths[seq_ID[is_del]] - local[is_del] + 1))

  # -- Drop mutations that fall inside an earlier mutation's span --

  # Sets a global end coordinates (clamping keeps mutation within its sequence)
  global_end_coord <- pos + mutation_widths - 1
  # Determines if a mutation event overlaps with another
  keep <- pos > c(-Inf, cummax(global_end_coord)[-total_mutations])
  # Applies the above filter
  seq_ID          <- seq_ID[keep]
  local           <- local[keep]
  value           <- value[keep]
  mutation_widths <- mutation_widths[keep]

  # Ensures a sequence whose only mutations were dropped are removed
  hit           <- unique(seq_ID)
  seq_ID_factor <- factor(seq_ID, levels = hit)

  # Generate IRanges for removed strings and their replacement string
  at  <- unname(S4Vectors::splitAsList(IRanges::IRanges(start = local, width = mutation_widths), seq_ID_factor))
  val <- unname(S4Vectors::splitAsList(value, seq_ID_factor))

  # Combines all of the above to insert every replacement string while removing correct bases
  # i.e. inserts all mutations
  DNA_seqs[hit] <- Biostrings::replaceAt(DNA_seqs[hit], at, val)

  return(as.character(DNA_seqs))
}


# Simulates a per-base sequencing error rate of 0.03% (by default)
#' @keywords internal
insert_sequencing_errors <- function(
  seqs,
  error_rate = 0.0003
) {
  # Convert to Biostrings DNA set
  DNA_seqs     <- if (is.character(seqs)) Biostrings::DNAStringSet(seqs) else seqs
  n_seqs       <- length(DNA_seqs)
  seq_lengths  <- as.numeric(IRanges::width(DNA_seqs))
  total_length <- sum(seq_lengths)

  # Generate the total number of errors based on total bases and error rate
  total_errors <- rbinom(1L, total_length, error_rate)
  if (total_errors == 0L) return(seqs)

  # Sample positions to mutate at, and determine the offsets
  pos     <- sort(sample.int(total_length, size = total_errors))
  offsets <- c(0, cumsum(seq_lengths)[-n_seqs])
  seq_ID  <- findInterval(pos - 1, offsets)
  local   <- as.integer(pos - offsets[seq_ID])

  # Errors must differ from the original base
  # Uses integer conversions and matching for efficiency
  orig <- as.character(Biostrings::subseq(DNA_seqs[seq_ID], start = local, width = 1L))
  idx  <- match(orig, BASES); idx[is.na(idx)] <- 0L
  alt  <- BASES[(idx - 1L + sample.int(3L, total_errors, replace = TRUE)) %% 4L + 1L]

  # Generate IRanges for the single base being replaced, and its replacement
  hit           <- unique(seq_ID)
  seq_ID_factor <- factor(seq_ID, levels = hit)
  at  <- unname(S4Vectors::splitAsList(IRanges::IRanges(start = local, width = 1L), seq_ID_factor))
  val <- unname(S4Vectors::splitAsList(alt, seq_ID_factor))

  # Combines all of the above to insert every replacement base while removing ref bases
  # i.e. inserts all substitutions
  DNA_seqs[hit] <- Biostrings::replaceAt(DNA_seqs[hit], at, val)

  return(as.character(DNA_seqs))
}


# Inserts identity barcodes into a template identity amplicon
#' @keywords internal
insert_identity_barcodes <- function(
  amplicon_seq,
  barcodes
) {
  empty_barcodes <- gregexpr("N+", amplicon_seq)[[1]]
  empty_barcode_starts <- as.integer(empty_barcodes)
  empty_barcode_ends   <- empty_barcode_starts + attr(empty_barcodes, "match.length") - 1L
  n_identity_barcodes  <- length(empty_barcode_starts)

  # Extract the sections around the empty barcodes
  seq_sections <- stringr::str_sub(amplicon_seq,
                                   c(1L, empty_barcode_ends + 1L),
                                   c(empty_barcode_starts - 1L, nchar(amplicon_seq)))

  inserted_barcodes <- seq_sections[1]
  for (i in seq_len(n_identity_barcodes)) {
    inserted_barcodes <- paste0(inserted_barcodes, barcodes[[paste0("BC", i)]], seq_sections[i + 1])
  }

  return(inserted_barcodes)
}


# Generates a vetor of read lengths using a truncated mixture model
# Has a tight peak around paired_read_length, with broad low-probability tails
#' @keywords internal
mixture_model_read_lengths <- function(
  total_reads,
  amplicon_length = 250,
  paired_read_length = 200,
  peak_prob = 0.99
) {
  # Draw from binomial to decide how many reads will have length drawn from tight peak/broad tail
  n_tight_peak <- rbinom(n = 1, size = total_reads, prob = peak_prob)
  n_broad_tail <- total_reads - n_tight_peak

  # Sample read lengths from (A) tight peak or (B) broad tail
  tight_peak <- truncnorm::rtruncnorm(n = n_tight_peak, a = 120, b = amplicon_length,
                                      mean = paired_read_length, sd = 0.5)
  broad_tail <- truncnorm::rtruncnorm(n = n_broad_tail, a = 120, b = amplicon_length,
                                      mean = paired_read_length, sd = 10)

  # Combine, round, shuffle
  read_lengths <- sample(as.integer(round(c(tight_peak, broad_tail))))

  return(read_lengths)
}


#' Used to determine the length of an insertion/deletion in a single mutation
#' code, or a vector of single mutation codes.
#' Currently UNUSED.
#' @keywords internal
calculate_indel_size <- function(mutation) {
  # Insertion - size is the length of the inserted sequence
  # e.g. c.100_101insGAGCTGT = 7 bases
  ins_seq <- stringi::stri_match_first_regex(
    mutation, "^c\\.\\d+_\\d+ins([ACGTN]+)$")[, 2]

  # Deletion - size is the difference between the base positions
  # e.g. c.146del     = 1 base
  # e.g. c.180_188del = 9 bases
  del_match <- stringi::stri_match_first_regex(
    mutation, "^c\\.(\\d+)(?:_(\\d+))?del[ACGTN]*$")
  del_start <- as.integer(del_match[, 2])
  del_end   <- as.integer(del_match[, 3])

  # Determine length
  indel_size <- dplyr::coalesce(
    stringi::stri_length(ins_seq),
    dplyr::case_when(
      is.na(del_start) ~ NA_integer_,
      is.na(del_end)   ~ 1L,
      .default         = del_end - del_start + 1L
    )
  )

  return(indel_size)
}
