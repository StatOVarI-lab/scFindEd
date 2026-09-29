#' Summarising Variants for Each Amplicon
#'
#' This function calculates the mutation codes HGVSc code for both the R1 and R2
#' reads on each amplicon in HGVSc nomenclature, using the data provided in the
#' amplicon's read fragments file. It exploits read pos, cigar, nm, md and seq
#' information to determine a single code that describes the difference between
#' each read and the reference amplicon.
#' The two HGVSc codes are then merged - if they are concordant in the overlapping
#' region - to produce a single code that describes a unique variant.
#' These variants are then collapsed per cell barcode, and the cell barcode is
#' matched to its corresponding identity barcode combination.
#'
#' @inheritParams global_params source_dir sample_name panel_ID barcoded output_dir verbose
#'
#' @return Exports the variant summary file, containing the HGVSc mutation code
#' that uniquely identifies each variant and its count, for each cell barcode.
#'
#' @importFrom data.table :=
#'
#' @export
summarise_variants <- function(
  source_dir,
  sample_name,
  panel_ID,
  barcoded = FALSE,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import SNP Summary --

  design_summary_file_name <- paste0(panel_ID, "_design_summary_joined.csv")
  design_summary_file <- import_design_summary_info(source_dir,
                                                    design_summary_file_name,
                                                    verbose)

  # -- Import Barcode System (optional) --

  if (barcoded) {
    barcode_system_file_name <- paste0(sample_name, "_barcode_system.csv")
    barcode_system <- import_barcode_system(source_dir,
                                            barcode_system_file_name,
                                            verbose)
  }

  # Extract unique amplicon IDs
  amplicon_list <- unique(design_summary_file$AmpID)
  n_amplicons <- length(amplicon_list)

  # Extract unique reference amplicon sequences
  reference_seqs <- unique(design_summary_file$amplicon_seq)

  # -- Repeat for each amplicon... --

  for (i in 1:n_amplicons) {
    AmpID <- amplicon_list[i]
    reference <- reference_seqs[i]

    # -- Import fragment mutations --

    if (verbose) cli::cli_h2(AmpID)

    fragment_mutations_file_name <- paste0(sample_name, "_", AmpID,
                                           "_fragment_mutations.tsv.gz")
    fragment_mutations <- import_fragment_mutations(paste0(source_dir, "/fragment_mutations"),
                                                    fragment_mutations_file_name,
                                                    verbose)

    # -- Filter out rows with NA or "" CIGAR strings --

    fragment_mutations <- fragment_mutations |>
      dplyr::filter(r1_cigar != "", !is.na(r1_cigar),
                    r2_cigar != "", !is.na(r2_cigar),
                    r1_pos   != "", !is.na(r1_pos),
                    r2_pos   != "", !is.na(r2_pos))

    # -- Calculate individual HGVSc codes for R1 reads --

    if (verbose) cli::cli_process_start(msg = "Calculating R1 HGVSc codes...",
                                        msg_done = "R1 HGVSc codes calculated.")

    fragment_mutations$r1_HGVSc <- calculate_HGVSc_codes(
      reference,
      fragment_mutations$r1_pos,
      fragment_mutations$r1_cigar,
      fragment_mutations$r1_nm,
      fragment_mutations$r1_md,
      fragment_mutations$r1_seq
    )

    # -- Calculate individual HGVSc codes for R2 reads --

    if (verbose) {
      cli::cli_process_done()
      cli::cli_process_start(msg = "Calculating R2 HGVSc codes...",
                             msg_done = "R2 HGVSc codes calculated.")
    }

    fragment_mutations$r2_HGVSc <- calculate_HGVSc_codes(
      reference,
      fragment_mutations$r2_pos,
      fragment_mutations$r2_cigar,
      fragment_mutations$r2_nm,
      fragment_mutations$r2_md,
      fragment_mutations$r2_seq
    )

    # -- Merge R1 and R2 HGVSc codes --

    if (verbose) {
      cli::cli_process_done()
      cli::cli_process_start(msg = "Merging R1 & R2 HGVSc codes...",
                             msg_done = "HGVSc codes merged.")
    }

    fragment_mutations$merged_HGVSc <- merge_HGVSc_codes(
      fragment_mutations$r1_HGVSc,
      fragment_mutations$r1_pos,
      fragment_mutations$r1_cigar,
      fragment_mutations$r2_HGVSc,
      fragment_mutations$r2_pos
    )

    if (verbose) cli::cli_process_done()

    # -- Collapse by Mutation Code --

    if (verbose) {
      cli::cli_process_done()
      cli::cli_process_start(msg = "Collapsing reads by shared HGVSc...",
                             msg_done = "Reads collapsed by shared HGVSc.")
    }

    # Group all reads in each cell that share the same merged_HGVSc and sum
    # the total fragments
    variant_summary <- fragment_mutations |>
      dplyr::group_by(sample_id, cell_barcode, amp_id, merged_HGVSc) |>
      dplyr::summarise(count = sum(n_fragments, na.rm = TRUE), .groups = "drop")

    # -- Assign Identity Barcode Combinations (optional) --

    if (barcoded) {
      if (verbose) {
        cli::cli_process_done()
        cli::cli_process_start(msg = "Linking cell barcodes to identity barcode combination...",
                               msg_done = "Cells linked to identity barcode combination.")
      }

      n_cells_pre_filter <- length(unique(variant_summary$cell_barcode))

      # Filter for only cell barcodes that have an assigned identity barcode combination
      # then match cell barcodes to attach the identity barcode column to variant_summary
      variant_summary <- variant_summary |>
        dplyr::filter(cell_barcode %in% barcode_system$cell_barcode) |>
        dplyr::left_join(barcode_system |>
                           dplyr::select(cell_barcode, all_BCs),
                         by = "cell_barcode") |>
        dplyr::rename(identity_barcode = all_BCs)

      n_cells_post_filter <- length(unique(variant_summary$cell_barcode))

      if (verbose) {
        cli::cli_process_done()
        cli::cli_alert_info("Number of cells retained after assigning identity barcodes: {n_cells_post_filter} out of {n_cells_pre_filter}.")
      }
    }

    # -- Optional file saving (only if output_dir is specified) --

    if (!is.null(output_dir)) {
      variant_summary_output_file_name <- paste0(sample_name, "_", AmpID,
                                                    "_variant_summary.tsv")

      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      readr::write_tsv(variant_summary, file = file.path(
        output_dir, variant_summary_output_file_name))

      cli::cli_alert_success(
        "{AmpID} variants with mutation codes saved to: {normalizePath(output_dir)}.")
    }
  }
}


#' Calculating HGVSc Codes for a Column Vector of Reads
#'
#' This function is used to generate a column vector of HGVSc codes that describe
#' the differences between observed reads and the reference amplicon sequence,
#' calculated from the follow data points: read sequence, start position, CIGAR
#' string, nm tag and md tag.
#'
#' @param reference_seq A singlular character containing the reference amplicon
#' sequence.
#' @param pos A column vector of read start positions.
#' @param cigar A column vector of read CIGAR strings.
#' @param nm A column vector of read nm (number of mismatches) tags.
#' @param md A column vector of read md (mismatching positions/bases) tags.
#' @param read_seq A column vector of read sequences.
#'
#' @return A column vector of calculated HGVSc codes.

#' @keywords internal
calculate_HGVSc_codes <- function(
  reference_seq,
  pos,
  cigar,
  nm,
  md,
  read_seq
) {
  # -- Set up data.table output --

  # Mute data.table progress bars for this function call only
  old_opts <- options(datatable.showProgress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  # Allocate a vector and populate with WT reads
  HGVSc_codes <- rep("WT", length(read_seq))

  # Filter for reads that actually contain mutations to avoid processing WTs
  mutation_indices <- which(nm > 0 & !is.na(cigar) & cigar != "")

  # Preventative fail-safe in case all rows are complete WTs
  if (length(mutation_indices) == 0) return(HGVSc_codes)

  # Build a data.table for only the rows with mutations
  dt_mutations <- data.table::data.table(
    row_ID   = mutation_indices,
    pos      = pos[mutation_indices],
    cigar    = cigar[mutation_indices],
    md       = md[mutation_indices],
    read_seq = read_seq[mutation_indices]
  )

  # -- MD Parsing --

  # Breaks down md into a vector of its respective chunks
  # Only break down unique mds, to reduce redundant regex calls
  unique_md <- unique(dt_mutations$md)
  md_matches_unique <- stringi::stri_extract_all_regex(unique_md, "[0-9]+|\\^[A-Z]+|[A-Z]")
  names(md_matches_unique) <- unique_md
  # Re-expand unique md matches out by lookup
  md_matches <- md_matches_unique[dt_mutations$md]
  lens_md    <- lengths(md_matches)

  md_parsed <- data.table::data.table(row_ID = rep(dt_mutations$row_ID, lens_md),
                                      pos    = rep(dt_mutations$pos, lens_md),
                                      chunk  = unlist(md_matches, use.names = FALSE))

  md_parsed[, first_char := substr(chunk, 1, 1)]
  # Checks for lengths (numbers)
  md_parsed[, is_num     := first_char %in% c("0","1","2","3","4","5","6","7","8","9")]
  # Checks for deletions ("^" presence)
  md_parsed[, is_del     := first_char == "^"]

  # Set reference step size
  md_parsed[, ref_step := 1L]
  md_parsed[is_num == TRUE, ref_step := as.integer(chunk)] # if number, step is the length
  md_parsed[is_del == TRUE, ref_step := nchar(chunk) - 1L] # if deletion, step is size of deletion

  # Track the movement along the reference position by cumsumming ref_steps
  md_parsed[, current_md_ref_pos := pos + vectorised_cumsum(ref_step, row_ID) - ref_step]

  # Create a separate data table for mismatches
  dt_mismatches <- md_parsed[is_num == FALSE & is_del == FALSE,
                             .(row_ID, target_ref_pos = current_md_ref_pos, ref_base = chunk)]

  # -- CIGAR Parsing --

  # Chunk CIGAR strings and flatten them into a single, long data table
  cigar_long <- construct_cigar_chunks_dt(dt_mutations$cigar, dt_mutations$row_ID, batch_size = 1000000)

  # Count the number of bases consumed relative to read / reference position for each chunk operation
  cigar_long[, `:=`(
    # Increase read position for "M" (match), "I" (insertion) and "S" (soft-clipping) chunks
    read_step = data.table::fifelse(op %in% c("M", "I", "S"), len, 0L),
    # Increase reference position for "M" (match) and "D" (deletion) chunks
    ref_step  = data.table::fifelse(op %in% c("M", "D"), len, 0L)
  )]

  # Calculate chunk start and end positions
  cigar_long[, `:=`(
    read_chunk_start = vectorised_cumsum(read_step, row_ID) - read_step + 1L,
    ref_chunk_start  = vectorised_cumsum(ref_step, row_ID) - ref_step
  )]

  # Merge the chunk positions into the mutation table
  cigar_long[dt_mutations, on = "row_ID", pos := i.pos]
  cigar_long[, `:=`(
    ref_chunk_start = ref_chunk_start + pos,
    read_chunk_end  = read_chunk_start + read_step - 1L
  )]
  cigar_long[, ref_chunk_end := ref_chunk_start + ref_step - 1L]

  # -- Calculate HGVSc codes for *mismatches* --

  # Filter out "M" (match) chunks and find which mismatches exist within these chunks
  M_chunks <- cigar_long[op == "M",
                         .(row_ID, read_chunk_start,
                           ref_chunk_start, ref_chunk_end,
                           lo = ref_chunk_start, hi = ref_chunk_end)]
  mismatch_mutations <- M_chunks[dt_mismatches,
                                 on = .(row_ID, lo <= target_ref_pos, hi >= target_ref_pos),
                                 nomatch = NULL]
  mismatch_mutations[, target_ref_pos := lo]
  mismatch_mutations[, c("lo", "hi") := NULL] # Clear intermediary columns

  # Calculate target_read_pos relative to the chunk
  mismatch_mutations[dt_mutations, on = "row_ID", read_seq := i.read_seq]
  mismatch_mutations[, target_read_pos := read_chunk_start + (target_ref_pos - ref_chunk_start)]
  mismatch_mutations[, read_base := substr(read_seq, target_read_pos, target_read_pos)]
  mismatch_mutations[, HGVSc := paste0("c.", target_ref_pos, ref_base, ">", read_base)]
  mismatch_HGVSc <- mismatch_mutations[, .(row_ID, HGVSc)]

  # -- Calculate HGVSc codes for *indels* --

  # "I" (insertions)
  insertion_mutations <- cigar_long[op == "I"]
  insertion_mutations[dt_mutations, on = "row_ID", read_seq := i.read_seq]
  insertion_mutations[, HGVSc := paste0("c.", (ref_chunk_start - 1L), "_", ref_chunk_start,
                                        "ins", substr(read_seq, read_chunk_start, read_chunk_end))]
  insertion_HGVSc <- insertion_mutations[, .(row_ID, HGVSc)]

  # "D" (deletions)
  deletion_mutations <- cigar_long[op == "D"]
  deletion_mutations[, HGVSc := data.table::fifelse(
    # Notation for deletion codes slightly different for 1 base and >1 base deletions
    len == 1L,
    paste0("c.", ref_chunk_start, "del"),
    paste0("c.", ref_chunk_start, "_", ref_chunk_end, "del")
  )]
  deletion_HGVSc <- deletion_mutations[, .(row_ID, HGVSc)]

  # -- Aggregate mutations --

  all_HGVSc <- data.table::rbindlist(list(mismatch_HGVSc, insertion_HGVSc, deletion_HGVSc))

  # Preventative fail-safe in case no mutations identified
  if (nrow(all_HGVSc) == 0) return(HGVSc_codes)

  # Track how many mutations exist per row_ID
  all_HGVSc[, n_mutations := .N, by = row_ID]

  # Assign HGVSc codes for reads with only 1 mutation
  singles <- all_HGVSc[n_mutations == 1]
  HGVSc_codes[singles$row_ID] <- singles$HGVSc

  # Assign HGVSc codes for reads with only >1 mutation (reduces number of non-vectorised paste() calls)
  multis <- all_HGVSc[n_mutations > 1]
  if (nrow(multis) > 0) {
    collapsed_multis <- multis[, .(merged_code = paste(HGVSc, collapse = ";")), by = row_ID]
    HGVSc_codes[collapsed_multis$row_ID] <- collapsed_multis$merged_code
  }

  return(HGVSc_codes)
}


#' Merging HGVSc Codes for Paired End Reads
#'
#' This function is used to generate a column vector of HGVSc codes that describe
#' the differences between observed read pairs and the reference amplicon sequence,
#' calculated by merging the individual HGVSc codes for the R1 and R2 of
#' paired-end reads. If the paired-end reads disagree in the overlapping region,
#' the function returns "Discordance in overlapping region" for that read pair.
#'
#' @param r1_HGVSc A column vector of calculated HGVSc codes for the R1 reads.
#' @param r1_pos A column vector of start positions for the R1 reads.
#' @param r1_cigar A column vector of CIGAR strings for the R1 reads.
#' @param r2_HGVSc A column vector of calculated HGVSc codes for the R2 reads.
#' @param r2_pos A column vector of start positions for the R2 reads.
#'
#' @return A column vector of merged HGVSc codes.

#' @keywords internal
merge_HGVSc_codes <- function(
  r1_HGVSc,
  r1_pos,
  r1_cigar,
  r2_HGVSc,
  r2_pos
) {
  # -- Set up data.table output --

  # Mute data.table progress bars for this function call only
  old_opts <- options(datatable.showProgress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  # Allocate a vector and populate with default WT status
  combined_HGVSc_codes <- rep("WT", length(r1_HGVSc))

  # Identify rows where at least one read has a mutation
  mutation_rows <- which(r1_HGVSc != "WT" | r2_HGVSc != "WT")
  if (length(mutation_rows) == 0) return(combined_HGVSc_codes)

  # Build a data.table for only the rows with mutations
  dt_mutations <- data.table::data.table(
    row_ID   = mutation_rows,
    r1_HGVSc = r1_HGVSc[mutation_rows],
    r1_pos   = r1_pos[mutation_rows],
    r1_cigar = r1_cigar[mutation_rows],
    r2_HGVSc = r2_HGVSc[mutation_rows],
    r2_pos   = r2_pos[mutation_rows]
  )

  # Calculate r1_end values
  dt_mutations[, r1_end := calculate_read_end_position(r1_pos, r1_cigar, batch_size = 1000000)]

  # -- Split HGVSc codes --

  # Split HGVSc codes for r1 reads
  chunks_r1 <- strsplit(dt_mutations$r1_HGVSc, ";", fixed = TRUE)
  lens_r1   <- lengths(chunks_r1)

  # Convert r1 mutations into long format (1 row per mutation in r1_HGVSc)
  r1_long   <- data.table::data.table(
    row_ID   = rep(dt_mutations$row_ID, lens_r1),
    r2_pos   = rep(dt_mutations$r2_pos, lens_r1),
    mutation = unlist(chunks_r1, use.names = FALSE),
    source   = "r1"
  )[mutation != "WT"]

  # Split HGVSc codes for r2 reads
  chunks_r2 <- strsplit(dt_mutations$r2_HGVSc, ";", fixed = TRUE)
  lens_r2   <- lengths(chunks_r2)

  # Convert r2 mutations into long format (1 row per mutation in r2_HGVSc)
  r2_long   <- data.table::data.table(
    row_ID   = rep(dt_mutations$row_ID, lens_r2),
    r1_end   = rep(dt_mutations$r1_end, lens_r2),
    mutation = unlist(chunks_r2, use.names = FALSE),
    source   = "r2"
  )[mutation != "WT"]

  # Combine long tables together
  mutations_long <- data.table::rbindlist(list(r1_long, r2_long), fill = TRUE)

  # Preventative fail-safe in case no mutations identified
  if (nrow(mutations_long) == 0) return(combined_HGVSc_codes)

  # -- Extract mutation positions --

  unique_mutations <- unique(mutations_long$mutation)
  pos_matrix_unique <- stringi::stri_match_first_regex(unique_mutations, "c\\.(\\d+)(?:_(\\d+))?")
  # Map each unique mutation row to mutations_long using index lookup
  match_idx <- match(mutations_long$mutation, unique_mutations)
  pos_matrix <- pos_matrix_unique[match_idx, , drop = FALSE]

  mutations_long[, mutation_start := as.integer(pos_matrix[, 2])]
  mutations_long[, mutation_end   := data.table::fifelse(!is.na(pos_matrix[, 3]),
                                                         as.integer(pos_matrix[, 3]),
                                                         mutation_start)]

  # -- Check for mutation presence in both read pairs --

  # Scans the table forward and backwards for mutations present in both pairs
  mutations_long[, in_opposite := duplicated(mutations_long, by = c("row_ID", "mutation")) |
                   duplicated(mutations_long, by = c("row_ID", "mutation"), fromLast = TRUE)]

  # -- Evaluate concordance --

  # Validates reads on r1 *before* the overlapping region, and reads on r2 *after* the overlapping region
  mutations_long[, is_valid := data.table::fcase(
    source == "r1" & (mutation_end < r2_pos | in_opposite), TRUE,
    source == "r2" & (mutation_start > r1_end | in_opposite), TRUE,
    default = FALSE
  )]

  # -- Aggregate tables --

  # Directly tag all rows with disagreeing mutations in overlapping region
  discordant_ids <- unique(mutations_long[is_valid == FALSE, row_ID])
  combined_HGVSc_codes[discordant_ids] <- "Discordance in overlapping region"

  # Filter to look exclusively at completely valid mutations
  valid_mutations <- mutations_long[!(row_ID %in% discordant_ids)]

  if (nrow(valid_mutations) > 0) {
    # Deduplicate mutations per read pair (mutations that agree will appear on both r1/r2)
    dt_valid <- unique(valid_mutations, by = c("row_ID", "mutation"))

    # Track how many mutations exist per row_ID
    dt_valid[, n_mutations := .N, by = row_ID]

    # Assign HGVSc codes for rows with only 1 mutation
    singles <- dt_valid[n_mutations == 1]
    combined_HGVSc_codes[singles$row_ID] <- singles$mutation

    # Assign HGVSc codes for rows with only >1 mutation (reduces number of paste() calls)
    multis <- dt_valid[n_mutations > 1]
    if (nrow(multis) > 0) {
      multi_collapsed <- multis[, .(combined_code = paste(mutation, collapse = ";")), by = row_ID]
      combined_HGVSc_codes[multi_collapsed$row_ID] <- multi_collapsed$combined_code
    }
  }

  return(combined_HGVSc_codes)
}




