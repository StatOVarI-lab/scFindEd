#' Generate Synthetic Experiment
#'
#' This function generates a complete synthetic scSNV-seq experiment, returning
#' they key files for scFindEd-nf processing (design summary, SNP list, R1 FASTQ,
#' R2 FASTQ), as well as a ground-truth table and a cell_barcode--readname
#' look-up table for validating results.
#'
#' *Complete description of intricacies needed.
#'
#' @inheritParams global_params output_dir verbose
#'
#' @param reads_per_amplicon Integer that fixes an approximate number of reads
#' on each amplicon.
#' @param total_cell_barcodes Integer that fixes the number of true cells.
#' @param total_whitelist_barcodes Integer that fixes the number of unique
#' barcodes in each identity barcode whitelist.
#' @param total_identity_barcode_combinations Integer that fixes the number of
#' unique identity barcode combinations that are spread across the true cells.
#' @param soup_fraction Float that defines a proportion of additional reads that
#' are simulate the 'soup' - reads that assign to both real and non-real cell
#' barcodes but don't represent the ground truth for that cell.
#' @param n_amplicons Integer that fixes the number of amplicons that are
#' generated in the synthetic design summary. Does not include the identity
#' amplicon.
#' @param n_edits Integer that fixes the number of edits that are generated in
#' thge synthetic SNP list. All edits are single nulceotide variants (SNVs).
#' @param total_identity_barcodes Integer that fixes the number of identity
#' barcodes inserted into the identity amplicon.
#' @param amplicon_length Integer that fixes the length (in bases) of each
#' amplicon sequence. Note that for simplicity, all amplicons will have the
#' same length.
#' @param paired_read_length Integer that fixes the average length (in bases) of
#' each of the forward and reverse reads in the simulated paired-end reads.
#' Lengths are drawn from a mixture model around mean `paired_read_length`. Note
#' that for simplicity, both the forward and reverse reads will draw from the
#' same average.
#' @param identity_barcode_length Integer that fixes the length (in bases) of
#' each of the inserted identity barcodes. Note that for simplicity, each
#' inserted barcode will have the same length.
#' @param random_seed Integer value used to set a consistent random seed.
#'
#' @return Six exported files: a synthetic design summary .csv, a synthetic SNP
#' list .csv, the R1 and R2 reads .fastq.gz files, and both a ground-truth
#' summary table as well as a cell_barcode--readname look-up table for results
#' validation.
#'
#' @export
generate_synthetic_experiment <- function(
  reads_per_amplicon = 1000000,
  total_cell_barcodes = 10000,
  total_whitelist_barcodes = 1000,
  total_identity_barcode_combinations = 500,
  soup_fraction = 0.2,
  n_amplicons = 5,
  n_edits = 20,
  total_identity_barcodes = 3,
  amplicon_length = 250,
  paired_read_length = 200,
  identity_barcode_length = 25,
  output_dir = NULL,
  random_seed = 42,
  verbose = TRUE
) {
  total_reads <- reads_per_amplicon * (n_amplicons + 1)

  # -- Ensure all provided arguments are valid --
  validate_generate_synthetic_experiment_args(total_reads,
                                              reads_per_amplicon,
                                              total_cell_barcodes,
                                              total_whitelist_barcodes,
                                              total_identity_barcode_combinations,
                                              soup_fraction,
                                              n_amplicons,
                                              n_edits,
                                              total_identity_barcodes,
                                              amplicon_length,
                                              paired_read_length,
                                              identity_barcode_length,
                                              output_dir,
                                              random_seed,
                                              verbose)

  # Prevent defaulting to scientific notation
  if (verbose) options(scipen = 999)

  # Set random seed
  set.seed(random_seed)

  # -- Generate Readnames --
  if (verbose) cli::cli_process_start(msg = "Generating {total_reads} readnames...",
                                      msg_done = "{total_reads} readnames generated.")
  readnames <- generate_readnames(total_reads)
  if (verbose) cli::cli_process_done()

  # -- Generate Cell Barcodes --

  if (verbose) cli::cli_process_start(msg = "Generating {total_cell_barcodes} cell barcodes...",
                                      msg_done = "{total_cell_barcodes} cell barcodes generated.")
  cell_barcodes <- generate_cell_barcodes(total_cell_barcodes, cell_barcode_length = 18)
  if (verbose) cli::cli_process_done()

  # -- Generate Identity Amplicon --

  if (verbose) cli::cli_process_start(msg = "Generating identity amplicon sequence...",
                                      msg_done = "Identity amplicon sequence generated.")
  identity_amplicon_seq <- generate_identity_amplicon(total_identity_barcodes,
                                                      amplicon_length,
                                                      identity_barcode_length)
  if (verbose) cli::cli_process_done()

  # -- Generate Design Summary --

  if (verbose) cli::cli_process_start(msg = "Generating design summary with {n_amplicons} amplicons...",
                                      msg_done = "Design summary with {n_amplicons} amplicons generated.")
  design_summary <- generate_design_summary(identity_amplicon_seq,
                                            n_amplicons,
                                            amplicon_length)
  if (verbose) cli::cli_process_done()

  # -- Generate SNP List --

  if (verbose) cli::cli_process_start(msg = "Generating SNP list with {n_edits} edits...",
                                      msg_done = "SNP list with {n_edits} edits generated.")
  SNP_list <- generate_SNP_list(design_summary,
                                n_amplicons,
                                n_edits)
  if (verbose) cli::cli_process_done()

  # -- Generate Identity Barcode Whitelists --

  identity_barcode_whitelists <- list()
  for (i in seq_len(total_identity_barcodes)) {
    if (verbose) cli::cli_process_start(msg = "Generating BC{i} whitelist with {total_whitelist_barcodes} barcodes...",
                                        msg_done = "BC{i} whitelist with {total_whitelist_barcodes} barcodes generated.")
    identity_barcode_whitelist <- generate_identity_barcode_whitelist(total_whitelist_barcodes,
                                                                      identity_barcode_length)

    identity_barcode_whitelists[[paste0("BC", i)]] <- identity_barcode_whitelist
    if (verbose) cli::cli_process_done()
  }

  # -- Generate Ground Truth Table --

  if (verbose) cli::cli_process_start(msg = "Generating ground truth table...",
                                      msg_done = "Ground truth table generated.")
  ground_truth_table <- generate_ground_truth_table(readnames,
                                                    cell_barcodes,
                                                    design_summary,
                                                    SNP_list,
                                                    identity_barcode_whitelists,
                                                    reads_per_amplicon,
                                                    total_whitelist_barcodes,
                                                    total_identity_barcode_combinations,
                                                    total_identity_barcodes)
  if (verbose) cli::cli_process_done()

  # -- Generate baseline amplicon sequences --

  if (verbose) cli::cli_process_start(msg = "Generating ground truth amplicon sequences for each identity barcode...",
                                      msg_done = "Ground truth amplicon sequences for each identity barcode generated.")
  baseline_amplicon_sequences <- generate_baseline_amplicon_sequences(design_summary,
                                                                      SNP_list,
                                                                      ground_truth_table)
  if (verbose) cli::cli_process_done()

  # -- Generate Soup reads --

  n_soup_reads <- pmax(1L, as.integer(round(rnorm(1, mean = total_reads * soup_fraction, sd = total_reads / 100))))
  soup_readnames <- generate_readnames(n_soup_reads, start = total_reads + 1L)

  if (verbose) cli::cli_process_start(msg = "Generating {n_soup_reads} additional soup reads...",
                                      msg_done = "{n_soup_reads} additional soup reads generated.")
  soup_reads <- generate_soup_reads(soup_readnames,
                                    cell_barcodes,
                                    design_summary,
                                    SNP_list,
                                    identity_barcode_whitelists,
                                    total_identity_barcodes,
                                    soup_edit_rate = 0.05)
  if (verbose) cli::cli_process_done()

  # -- Generate FASTQ reads --

  if (verbose) cli::cli_process_start(msg = "Generating complete FASTQ reads...",
                                      msg_done = "Complete FASTQ reads generated.")
  FASTQ_reads <- generate_FASTQ_reads(ground_truth_table,
                                      baseline_amplicon_sequences,
                                      soup_reads,
                                      amplicon_length,
                                      paired_read_length)
  if (verbose) cli::cli_process_done()

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    # Set output file names
    design_summary_file_name        <- "SYNTHETIC_design_summary.csv"
    SNP_list_file_name              <- "SYNTHETIC_SNP_list.csv"
    ground_truth_table_file_name    <- "ground_truth_table.csv"
    cell_barcode_readname_file_name <- "cell_barcode_readname_lookup.tsv"
    forward_FASTQ_file_name         <- "SYNTHETIC_R1_reads.fastq.gz"
    reverse_FASTQ_file_name         <- "SYNTHETIC_R2_reads.fastq.gz"

    # Export design summary file
    readr::write_csv(design_summary,
                     file = file.path(output_dir, design_summary_file_name))
    cli::cli_alert_success("Synthetic design summary saved to: {normalizePath(output_dir)}.")

    # Export SNP list file
    readr::write_csv(SNP_list |> dplyr::select(ID, Gene, SNP, gRNA, HDR_template),
                     file = file.path(output_dir, SNP_list_file_name))
    cli::cli_alert_success("Synthetic SNP list saved to: {normalizePath(output_dir)}.")

    # Export Barcode Whitelists
    for (i in seq_len(total_identity_barcodes)) {
      identity_barcode_whitelist_file_name <- paste0("BC", i, "_whitelist_barcodes.csv")
      readr::write_csv(tibble::tibble(identity_barcode_whitelists[[paste0("BC", i)]]),
                       file = file.path(output_dir, identity_barcode_whitelist_file_name),
                       col_names = FALSE)
      cli::cli_alert_success("BC{i} whitelist barcodes saved to: {normalizePath(output_dir)}.")
    }

    # Export ground truth table
    readr::write_csv(ground_truth_table |> dplyr::select(-readname) |> dplyr::distinct(),
                     file = file.path(output_dir, ground_truth_table_file_name))
    cli::cli_alert_success("Ground truth table saved to: {normalizePath(output_dir)}.")

    # Export cell_barcode--readname lookup table
    readr::write_tsv(ground_truth_table |> dplyr::select(readname, cell_barcode),
                     file = file.path(output_dir, cell_barcode_readname_file_name))
    cli::cli_alert_success("cell_barcode--readname look-up table saved to: {normalizePath(output_dir)}.")

    # Export forward FASTQ reads
    data.table::fwrite(list(as.vector(t(as.matrix(
      FASTQ_reads |> dplyr::select(FASTQ_header, fwd_seq, separator, fwd_quality)
    )))), file = file.path(output_dir, forward_FASTQ_file_name), col.names = FALSE)
    cli::cli_alert_success("Forward FASTQ reads saved to: {normalizePath(output_dir)}.")

    # Export reverse FASTQ reads
    data.table::fwrite(list(as.vector(t(as.matrix(
      FASTQ_reads |> dplyr::select(FASTQ_header, rev_seq, separator, rev_quality)
    )))), file = file.path(output_dir, reverse_FASTQ_file_name), col.names = FALSE)
    cli::cli_alert_success("Reverse FASTQ reads saved to: {normalizePath(output_dir)}.")
  }
}


#' @keywords internal
generate_readnames <- function(
  total_reads,
  start = 1
) {
  # Generate x & y coordinates
  coords <- seq(from = start, length.out = total_reads)
  y_coords <- coords %% 10000L
  x_coords <- as.integer(ceiling(coords / 10000))

  # Bind values into valid readnames
  readnames <- stringi::stri_join(
    "INSTRUMENT", "RUN", "FLOWCELL", "LANE", "TILE",
    stringi::stri_pad_left(x_coords, width = 5, pad = "0"),
    stringi::stri_pad_left(y_coords, width = 5, pad = "0"),
    sep = ":"
  )

  return(readnames)
}


#' @keywords internal
generate_cell_barcodes <- function(
  total_cell_barcodes = 10000,
  cell_barcode_length = 18
) {
  cell_barcodes <- stringi::stri_rand_strings(
    n = total_cell_barcodes, length = cell_barcode_length,
    pattern = paste0("[", paste(BASES, collapse = ""), "]")
  ) |> stringi::stri_join("-1")

  return(cell_barcodes)
}


#' @keywords internal
generate_identity_amplicon <- function(
  total_identity_barcodes = 3,
  amplicon_length = 250,
  identity_barcode_length = 25
) {
  # Generate initial sequence
  identity_amplicon_seq <- markov_DNA_generation(amplicon_length)

  for (i in 1:total_identity_barcodes) {
    # Calculate start and end positions of barcode(s)
    insert_start <- 21 + (i-1)*(identity_barcode_length + 20)
    insert_end <- insert_start + identity_barcode_length - 1

    # Generate and insert empty barcode(s) - Ns
    empty_barcode <- stringi::stri_dup("N", identity_barcode_length)
    stringi::stri_sub(identity_amplicon_seq, insert_start, insert_end) <- empty_barcode
  }

  return(identity_amplicon_seq)
}


#' @keywords internal
generate_design_summary <- function(
  identity_amplicon_seq,
  n_amplicons = 5,
  amplicon_length = 250
) {
  # Generate position metadata
  first_AmpID_no <- sample(100000:999000, 1)
  chrs <- paste0("chr", sort(sample(1:20, n_amplicons, replace = TRUE)))
  amplicon_starts <- sample(10000000:99999000, n_amplicons)
  amplicon_ends <- amplicon_starts + amplicon_length - 1
  insert_starts <- amplicon_starts + 30
  insert_ends <- amplicon_ends - 30

  # Generate amplicon sequences
  amplicon_seqs <- sapply(1:n_amplicons, function(x) markov_DNA_generation(amplicon_length))

  # Construct design summary tibble
  design_summary <- dplyr::tibble(
    AmpID          = paste0("AMPL", seq(from = first_AmpID_no, to = first_AmpID_no + n_amplicons - 1)),
    chr            = chrs,
    amplicon_start = amplicon_starts,
    insert_start   = insert_starts,
    insert_end     = insert_ends,
    amplicon_end   = amplicon_ends,
    Ref_Genome     = rep("hg38_GR", n_amplicons),
    fwd_seq        = stringi::stri_sub(amplicon_seqs, 1, 30),
    rev_seq        = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(
      stringi::stri_sub(amplicon_seqs, amplicon_length - 30 + 1, amplicon_length)
    ))),
    amplicon_seq   = amplicon_seqs,
    review         = rep(0, n_amplicons),
    Targets        = paste0(chr, ":", insert_start, "-", insert_end)
  )

  # Add Identity Amplicon Row
  identity_amplicon_start <- sample(250:750, 1)
  identity_amplicon_end   <- identity_amplicon_start + amplicon_length - 1
  identity_insert_start   <- identity_amplicon_start + 20
  identity_insert_end     <- identity_amplicon_end - 20

  identity_amplicon_row <- dplyr::tibble(
    AmpID          = "AMPL999999",
    chr            = "chrBarcode",
    amplicon_start = identity_amplicon_start,
    insert_start   = identity_insert_start,
    insert_end     = identity_insert_end,
    amplicon_end   = identity_amplicon_end,
    Ref_Genome     = "hg38_GR",
    fwd_seq        = stringi::stri_sub(identity_amplicon_seq, 1, 20),
    rev_seq        = as.character(Biostrings::reverseComplement(Biostrings::DNAString(
      stringi::stri_sub(identity_amplicon_seq, amplicon_length - 20 + 1, amplicon_length)
    ))),
    amplicon_seq   = identity_amplicon_seq,
    review         = 0,
    Targets        = paste0("chrBarcode:", identity_insert_start, "-", identity_insert_end)
  )

  return(dplyr::bind_rows(design_summary, identity_amplicon_row))
}


#' @keywords internal
generate_SNP_list <- function(
  design_summary,
  n_amplicons = 5,
  n_edits = 20
) {
  gene_IDs <- paste0("Gene", seq_len(n_amplicons))
  design_summary[["Gene"]] <- c(gene_IDs, NA)

  # Sample Gene IDs and generate SNP positions
  SNP_list <- dplyr::tibble(
    Gene    = sample(gene_IDs, n_edits, replace = TRUE),
    SNP_loc = sample(75:175, n_edits)
  ) |>
    dplyr::arrange(as.numeric(stringi::stri_sub(Gene, 5)), SNP_loc) |>
    dplyr::mutate(ID = seq_len(n_edits)) |>
    # Bring across amplicon_seq from design summary
    dplyr::left_join(
      design_summary |> dplyr::select(Gene, AmpID, amplicon_seq),
      by = "Gene") |>
    # Extract ref base and generate alt base
    dplyr::mutate(ref_base = stringi::stri_sub(amplicon_seq, from = SNP_loc, length = 1)) |>
    dplyr::rowwise() |>
    dplyr::mutate(alt_base = sample(setdiff(BASES, ref_base), 1)) |>
    dplyr::ungroup() |>
    # Insert edit to amplicon_seq
    dplyr::mutate(edited_seq = stringi::stri_sub_replace(amplicon_seq,
                                                         from = SNP_loc,
                                                         to = SNP_loc,
                                                         replacement = alt_base)) |>
    # Cut out section of edited_seq for HDR template
    dplyr::mutate(HDR_template = stringi::stri_sub(edited_seq,
                                                   from = SNP_loc - 50,
                                                   to = SNP_loc + 49)) |>
    # Cut out section of edited_seq for gRNA
    dplyr::mutate(gRNA = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(
      stringi::stri_sub(edited_seq, from = SNP_loc - 25, to = SNP_loc - 6)
    )))) |>
    # Add unique SNP identifier
    dplyr::group_by(Gene) |>
    dplyr::mutate(SNP = paste0("SNP", dplyr::row_number())) |>
    dplyr::ungroup() |>
    # Select and re-order necessary columns
    dplyr::select(ID, Gene, SNP, AmpID, gRNA, HDR_template, SNP_loc,
                  ref_base, alt_base, amplicon_seq, edited_seq)

  return(SNP_list)
}


#' @keywords internal
generate_identity_barcode_whitelist <- function(
  total_whitelist_barcodes = 1000,
  identity_barcode_length = 25
) {
  barcode_whitelist <- character(0)

  while (length(barcode_whitelist) < total_whitelist_barcodes) {
    needed <- total_whitelist_barcodes - length(barcode_whitelist)
    new_barcodes <- sapply(1:needed, function(x) markov_DNA_generation(identity_barcode_length))
    barcode_whitelist <- unique(c(barcode_whitelist, new_barcodes))
  }

  return(sort(barcode_whitelist))
}


#' @keywords internal
generate_ground_truth_table <- function(
  readnames,
  cell_barcodes,
  design_summary,
  SNP_list,
  identity_barcode_whitelists,
  reads_per_amplicon = 100000,
  total_whitelist_barcodes = 1000,
  total_identity_barcode_combinations = 500,
  total_identity_barcodes = 3
) {
  # -- Generate identity barcode combinations --

  identity_barcode_combinations <- generate_identity_barcode_combinations(
    identity_barcode_whitelists,
    total_whitelist_barcodes,
    total_identity_barcode_combinations,
    total_identity_barcodes
  )

  # -- Generate lookup table between all_BCs and edit --

  identity_barcode_edit_assignment <- generate_identity_barcode_edit_assignment(
    identity_barcode_combinations[["all_BCs"]],
    SNP_list
  )

  # -- Generate lookup table between cell_barcode and identity_barcode (all_BCs) --
  cell_identity_barcode_assignment <- skewed_cell_identity_barcode_assignment(
    cell_barcodes,
    identity_barcode_combinations[["all_BCs"]],
    skewness = 0.4
  )

  # -- Generate lookup table between readnames and cell_barcode --
  readname_cell_barcode_assignment <- skewed_readname_cell_barcode_assignment(
    readnames,
    cell_barcodes,
    skewness = 0.7
  ) |>
    dplyr::mutate(AmpID = rep(design_summary[["AmpID"]], each = reads_per_amplicon))

  # Construct identity barcode side of the ground truth table
  ground_truth_table <- readname_cell_barcode_assignment |>
    dplyr::left_join(cell_identity_barcode_assignment, by = "cell_barcode") |>
    dplyr::left_join(identity_barcode_combinations, by = "all_BCs") |>
    dplyr::left_join(identity_barcode_edit_assignment, by = "all_BCs")

  return(ground_truth_table)
}


#' @keywords internal
generate_identity_barcode_combinations <- function(
  identity_barcode_whitelists,
  total_whitelist_barcodes = 1000,
  total_identity_barcode_combinations = 500,
  total_identity_barcodes = 3
) {
  # Create empty table of correct length for storing true barcode combinations
  identity_barcode_combinations <- dplyr::tibble(.rows = total_identity_barcode_combinations)

  # For each identity barcode, sample a vector of barcodes from respective whitelist
  for (i in 1:total_identity_barcodes) {
    barcode_ID <- paste0("BC", i)

    sampled_barcode_vals <- sample(1:total_whitelist_barcodes,
                                   total_identity_barcode_combinations,
                                   replace = TRUE)
    sampled_barcodes <- identity_barcode_whitelists[[barcode_ID]][sampled_barcode_vals]

    identity_barcode_combinations <- identity_barcode_combinations |>
      dplyr::mutate(!!rlang::sym(barcode_ID) := sampled_barcodes)

    # Generate a single code (all_BCs) that details the identity barcode combination
    if (i == 1) {
      all_BCs <- paste0("BC", i, "_", sampled_barcode_vals)
    } else {
      all_BCs <- paste0(all_BCs, "_BC", i, "_", sampled_barcode_vals)
    }
  }

  identity_barcode_combinations <- dplyr::bind_cols(all_BCs = all_BCs,
                                                    identity_barcode_combinations)

  return(identity_barcode_combinations)
}


#' @keywords internal
generate_identity_barcode_edit_assignment <- function(
  identity_barcodes,
  SNP_list
) {
  n_identity_barcodes <- length(identity_barcodes)
  n_wildtype          <- floor(n_identity_barcodes / 2)
  n_edited            <- n_identity_barcodes - n_wildtype
  SNP_IDs <- sample(SNP_list[["ID"]])

  # Randomised list of SNP ID genotypes, with half as wildtype
  assignments <- c(rep("WT", n_wildtype),
                   as.character(sample(SNP_IDs, n_edited, replace = TRUE)))

  # Assign SNP IDs to identity_barcodes at random
  identity_barcode_edit_assignment <- dplyr::tibble(all_BCs = identity_barcodes) |>
    dplyr::mutate(SNP_ID = sample(assignments))

  return(identity_barcode_edit_assignment)
}


#' @keywords internal
skewed_cell_identity_barcode_assignment <- function(
  cell_barcodes,
  identity_barcodes,
  skewness = 0.4
) {
  n_cell_barcodes     <- length(cell_barcodes)
  n_identity_barcodes <- length(identity_barcodes)

  # Shuffle cell_barcodes to ensure randomness
  shuffled_cells <- sample(cell_barcodes)

  guaranteed_assignments <- identity_barcodes
  remaining_cells <- shuffled_cells[(n_identity_barcodes + 1):n_cell_barcodes]

  # Generate skewed probabilities using a Gamma distribution
  raw_weights <- rgamma(n_identity_barcodes, shape = skewness, rate = 1)
  probs       <- raw_weights / sum(raw_weights)

  # Sample identity barcodes for all remaining cell barcodes
  remaining_assignments <- sample(identity_barcodes,
                                  size    = length(remaining_cells),
                                  replace = TRUE,
                                  prob    = probs
  )

  # Combine guaranteed and remaining assignments into single lookup table
  cell_identity_barcode_assignment <- dplyr::tibble(
    cell_barcode = c(shuffled_cells[1:n_identity_barcodes], remaining_cells),
    all_BCs      = c(guaranteed_assignments, remaining_assignments)
  )

  return(cell_identity_barcode_assignment)
}


#' @keywords internal
skewed_readname_cell_barcode_assignment <- function(
  readnames,
  cell_barcodes,
  skewness = 0.7
) {
  n_readnames     <- length(readnames)
  n_cell_barcodes <- length(cell_barcodes)

  # Shuffle readnames to ensure randomness
  shuffled_readnames <- sample(readnames)

  guaranteed_assignments <- cell_barcodes
  remaining_readnames <- shuffled_readnames[(n_cell_barcodes + 1):n_readnames]

  # Generate skewed probabilities using a Gamma distribution
  raw_weights <- rgamma(n_cell_barcodes, shape = skewness, rate = 1)
  probs       <- raw_weights / sum(raw_weights)

  # Sample identity barcodes for all remaining cell barcodes
  remaining_assignments <- sample(cell_barcodes,
                                  size    = length(remaining_readnames),
                                  replace = TRUE,
                                  prob    = probs
  )

  # Combine guaranteed and remaining assignments into single lookup table
  readname_cell_barcode_assignment <- dplyr::tibble(
    readname     = c(shuffled_readnames[1:n_cell_barcodes], remaining_readnames),
    cell_barcode = c(guaranteed_assignments, remaining_assignments)
  )

  return(readname_cell_barcode_assignment)
}


#' @keywords internal
generate_baseline_amplicon_sequences <- function(
  design_summary,
  SNP_list,
  ground_truth_table
) {
  n_identity_barcodes = length(unique(ground_truth_table[["all_BCs"]]))

  # Create table with key unique information from constructing baseline sequences
  baseline_amplicon_sequences <- ground_truth_table |>
    dplyr::distinct(all_BCs, SNP_ID, .keep_all = TRUE) |>
    dplyr::select(all_BCs, SNP_ID, dplyr::starts_with("BC"))

  AmpIDs <- design_summary[["AmpID"]]

  for (ampID in AmpIDs) {
    # Get the wildtype amplicon sequence for ampID
    wildtype_seq <- design_summary[["amplicon_seq"]][design_summary[["AmpID"]] == ampID]

    # Edit amplicons
    if (ampID != "AMPL999999") {
      # Initialise all rows with the wildtype amplicon sequence
      amplicon_sequences <- rep(wildtype_seq, times = nrow(baseline_amplicon_sequences))

      # Identify which rows in baseline_amplicon_sequences have SNPs on this AmpID
      SNPs_in_AmpID <- SNP_list[["ID"]][SNP_list[["AmpID"]] == ampID]
      edit_mask <- baseline_amplicon_sequences[["SNP_ID"]] %in% SNPs_in_AmpID

      # For rows matching edit_mask, lookup the matching edited sequence from SNP_list
      target_SNPs <- baseline_amplicon_sequences[["SNP_ID"]][edit_mask]
      SNP_matches <- match(target_SNPs, SNP_list[["ID"]])
      amplicon_sequences[edit_mask] <- SNP_list[["edited_seq"]][SNP_matches]

      # Identity amplicon
    } else {
      # Insert respective barcode combinations into each amplicon_seq
      amplicon_sequences <- insert_identity_barcodes(
        wildtype_seq,
        baseline_amplicon_sequences |> dplyr::select(dplyr::starts_with("BC"))
      )
    }

    # Insert TRUE mutations into baseline sequences at per-base mutation_rate
    amplicon_sequences <- mutate_seqs(amplicon_sequences, mutation_rate = 0.001)

    baseline_amplicon_sequences[[ampID]] <- amplicon_sequences
  }

  return(baseline_amplicon_sequences)
}


#' @keywords internal
generate_soup_reads <- function(
  readnames,
  cell_barcodes,
  design_summary,
  SNP_list,
  identity_barcode_whitelists,
  total_identity_barcodes = 3,
  soup_edit_rate = 0.05
) {
  n_soup_reads <- length(readnames)
  AmpIDs       <- design_summary[["AmpID"]]

  # -- Sample cell_barcodes: quarter from true pool, three quarter generated --

  n_real_cells <- n_soup_reads %/% 4
  n_soup_cells <- n_soup_reads - n_real_cells
  soup_cell_barcodes <- sample(c(
    sample(cell_barcodes, n_real_cells, replace = TRUE),
    generate_cell_barcodes(n_soup_cells, cell_barcode_length = 18)
  ))

  # Assign each soup read to an amplicon at random
  soup_AmpID     <- sample(AmpIDs, n_soup_reads, replace = TRUE)
  soup_sequences <- character(n_soup_reads)

  for (ampID in AmpIDs) {
    rows <- which(soup_AmpID == ampID)
    if (length(rows) == 0) next

    # Pull out the wildtype sequence for ampID
    wildtype_seq <- design_summary[["amplicon_seq"]][design_summary[["AmpID"]] == ampID]

    # Edit amplicons
    if (ampID != "AMPL999999") {
      # Initialise all rows with the wildtype amplicon sequence
      amplicon_sequences <- rep(wildtype_seq, times = length(rows))

      # Identify which rows in baseline_amplicon_sequences have SNPs on this AmpID
      SNPs_in_AmpID <- SNP_list[["ID"]][SNP_list[["AmpID"]] == ampID]

      if (length(SNPs_in_AmpID) > 0) {
        # Sample reads to contain an edit
        has_edit <- runif(length(rows)) < soup_edit_rate

        if (any(has_edit)) {
          insert_SNPs <- sample(SNPs_in_AmpID, sum(has_edit), replace = TRUE)
          amplicon_sequences[has_edit] <- SNP_list[["edited_seq"]][match(insert_SNPs, SNP_list[["ID"]])]
        }
      }

      soup_sequences[rows] <- amplicon_sequences

      # Identity amplicon
    } else {
      # Sample a random barcode combination
      rand_barcodes <- dplyr::tibble(.rows = length(rows))

      for (i in seq_len(total_identity_barcodes)) {
        BC_ID <- paste0("BC", i)
        rand_barcodes[[BC_ID]] <- sample(identity_barcode_whitelists[[BC_ID]],
                                         length(rows), replace = TRUE)
      }

      soup_sequences[rows] <- insert_identity_barcodes(wildtype_seq, rand_barcodes)
    }
  }

  # Insert mutations at the same background rate as baseline sequences
  soup_sequences <- mutate_seqs(soup_sequences, mutation_rate = 0.001)

  soup_reads <- dplyr::tibble(readname     = readnames,
                              cell_barcode = soup_cell_barcodes,
                              seq          = soup_sequences)

  return(soup_reads)
}


#' @keywords internal
generate_FASTQ_reads <- function(
  ground_truth_table,
  baseline_amplicon_sequences,
  soup_reads = NULL,
  amplicon_length = 250,
  paired_read_length = 200
) {
  # Create a data frame for all FASTQ reads from the ground_truth_table
  FASTQ_reads <- ground_truth_table |>
    dplyr::select(readname, cell_barcode, AmpID, all_BCs) |>
    dplyr::left_join(baseline_amplicon_sequences |>
                       tidyr::pivot_longer(
                         cols = dplyr::starts_with("AMPL"),
                         names_to = "AmpID",
                         values_to = "seq"
                       ) |> dplyr::select(all_BCs, AmpID, seq),
                     by = c("all_BCs", "AmpID"))

  # Stack FASTQ_reads with soup_reads
  FASTQ_reads <- dplyr::bind_rows(FASTQ_reads, soup_reads)

  total_reads <- nrow(FASTQ_reads)

  # Draw read lengths from a tight peak + broad tail mixture model
  forward_read_lengths <- mixture_model_read_lengths(total_reads,
                                                     amplicon_length,
                                                     paired_read_length,
                                                     peak_prob = 0.99)

  reverse_read_lengths <- mixture_model_read_lengths(total_reads,
                                                     amplicon_length,
                                                     paired_read_length,
                                                     peak_prob = 0.99)

  # Slice actual sequences into forward and reverse reads based on length
  FASTQ_reads <- FASTQ_reads |>
    dplyr::select(readname, cell_barcode, seq) |>
    dplyr::mutate(fwd_length = forward_read_lengths,
                  rev_length = reverse_read_lengths) |>
    dplyr::mutate(fwd_seq = stringi::stri_sub(seq, from = 1, to = fwd_length),
                  rev_seq = stringi::stri_sub(seq, from = -rev_length, to = -1)) |>
    dplyr::select(-seq)

  # Simulates a per-base sequencing error rate of 0.03%
  FASTQ_reads[["fwd_seq"]] <- insert_sequencing_errors(FASTQ_reads[["fwd_seq"]],
                                                       error_rate = 0.0003)
  FASTQ_reads[["rev_seq"]] <- insert_sequencing_errors(FASTQ_reads[["rev_seq"]],
                                                       error_rate = 0.0003)

  # Modify columns for FASTQ output purposes
  FASTQ_reads <- FASTQ_reads |>
    dplyr::mutate(FASTQ_header = stringi::stri_c("@", readname, " RG:Z:", cell_barcode),
                  separator = rep("+", times = total_reads),
                  fwd_quality = stringi::stri_dup("N", fwd_length),
                  rev_quality = stringi::stri_dup("N", rev_length)) |>
    dplyr::select(FASTQ_header, fwd_seq, rev_seq, separator, fwd_quality, rev_quality)

  return(FASTQ_reads)
}

