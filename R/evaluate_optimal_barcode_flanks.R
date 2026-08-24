#' Evaluate Optimal Barcode Flanks
#'
#' This function takes the extracted identity barcode information as well as the
#' sequenced reads on `identity_AmpID` in order to determine the best choice of
#' flank for each identity barocde.
#'
#' The best choice of flank is determined as the flank that: has 0 reads where
#' the flank has multiple hits and has the most reads with a single hit, with
#' ties broken by selecting the shortest flank size.
#'
#' This function prints the optimal flank for each barcode, but also produces a
#' .csv summary file for each barcode so that users may evaluate each flank
#' possibility for themselves.
#'
#' @inheritParams global_params source_dir sample_name identity_AmpID min_reads output_dir verbose
#'
#' @param min_flank_size An integer that defines the minimum flank size to be
#' evaluated. Cannot be less than 5.
#' @param max_flank_size An integer that defines the maximum flank size to be
#' evaluated.
#'
#' @return The optimal flank (flank size, orientation - forward/reverse read,
#' direction - upstream/downstream) is suggested and added to
#' `{sample_name}_identity_barcode_info.csv`. Summary stats for each possible
#' flank are also saved to `{sample_name}_BC{ID}_flank_stats.csv`.
#'
#' @export
evaluate_optimal_barcode_flanks <- function(
  source_dir,
  sample_name,
  identity_AmpID,
  min_reads = 2,
  min_flank_size = 5,
  max_flank_size = 20,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import datasets --

  identity_amplicon_reads_file_name <- paste0(sample_name, "_", identity_AmpID,
                                              "_reads.tsv.gz")
  identity_amplicon_reads <- import_identity_amplicon_reads(source_dir,
                                                            identity_amplicon_reads_file_name,
                                                            verbose)

  identity_barcode_info_file_name <- paste0(sample_name, "_identity_barcode_info.csv")
  identity_barcode_info <- import_identity_barcode_info(source_dir,
                                                        identity_barcode_info_file_name,
                                                        verbose)

  # If cell_barcode is mis-named as barcode, convert the column name
  if ("barcode" %in% colnames(identity_amplicon_reads)) {
    colnames(identity_amplicon_reads)[colnames(identity_amplicon_reads) == "barcode"] <- "cell_barcode"
  }

  # -- Filter reads according to min_reads threshold --

  invisible(list2env(filter_identity_amplicon_reads(
    identity_amplicon_reads,
    min_reads,
    verbose = FALSE
  ), envir = environment()))

  # -- Set Up Barcode Analysis --

  if (!is.numeric(min_flank_size) | !is.numeric(max_flank_size)) {
    cli::cli_abort("min_flank_size and max_flank_size must be numeric.")
  } else if (min_flank_size >= max_flank_size) {
    cli::cli_abort("min_flank_size must be less than max_flank_size.")
  } else if (min_flank_size < 5) {
    cli::cli_abort("The minimum acceptable barcode flank is 5 bases.")
  }

  flank_sizes <- seq(min_flank_size, max_flank_size)

  # Combine filtered forward and reverse reads
  reads_by_orientation <- list(forward = identity_amplicon_reads_forward,
                               reverse = identity_amplicon_reads_reverse)
  directions <- c("upstream", "downstream")

  # Create a grid of every (barcode x orientation x direction) combination
  plan <- tidyr::expand_grid(
    barcode_row = seq_len(nrow(identity_barcode_info)),
    orientation = names(reads_by_orientation),
    direction   = directions
  )

  # -- Compute hits / multiple-hit stats for every combination --

  all_stats <- purrr::pmap_dfr(plan, function(barcode_row, orientation, direction) {
    barcode_No <- identity_barcode_info$barcode_ID[barcode_row]

    if (verbose) cli::cli_process_start(msg = "BC{barcode_No}: Assessing {orientation} {direction} flanks...",
                                        msg_done = "BC{barcode_No}: {orientation} {direction} flanks assessed.")

    # Test flanks on both reads, in both directions, from min_flank_size to max_flank_size
    stats <- test_flanks(
      reads_by_orientation[[orientation]],
      identity_barcode_info[barcode_row, ],
      direction,
      min_flank_size,
      max_flank_size
    )

    if (verbose) cli::cli_process_done()

    # Record stats in combined dataframe
    tibble::tibble(
      barcode_row   = barcode_row,
      barcode_ID    = barcode_No,
      orientation   = orientation,
      direction     = direction,
      flank_size    = flank_sizes,
      hits          = stats[["hits_vec"]],
      multiple_hits = stats[["multiple_hits_vec"]]
    )
  })

  # -- Optional file saving (only if output_dir is provided) --
  # -- One .CSV per barcode --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # Structure output so each barcode's _hits and _multiple_hits are adjacent
    combo_bases <- paste(
      rep(names(reads_by_orientation), each = length(directions)),
      directions,
      sep = "_"
    )
    ordered_cols <- as.vector(rbind(paste0(combo_bases, "_hits"),
                                    paste0(combo_bases, "_multiple_hits")))


    for (ID in unique(all_stats$barcode_ID)) {
      export_flank_test_stats(ID, all_stats, ordered_cols, sample_name, output_dir)
    }
  }

  # -- Consensus Suggestions --

  total_reads <- nrow(identity_amplicon_reads_forward)

  # Rule: among flank lengths with 0 multiple hits, select the flank with the
  # most single hits. Break ties by selecting the shorter flank.
  consensus <- all_stats |>
    dplyr::filter(multiple_hits == 0) |>
    dplyr::group_by(barcode_ID) |>
    dplyr::slice_max(hits, n = 1, with_ties = TRUE) |>
    dplyr::slice_min(flank_size, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Assess if any barcodes have no flank lengths that prevent multiple hits
  failed_barcodes <- setdiff(identity_barcode_info$barcode_ID, consensus$barcode_ID)
  if (length(failed_barcodes) > 0) {
    cli::cli_alert_warning(
      "No flank size with 0 multiple hits found for BC{failed_barcodes} within size range {min_flank_size}-{max_flank_size}."
    )
  }

  if (verbose) {
    # Print the best barcode flank length and stats for each barcode
    for (row in seq_len(nrow(consensus))) {
      rec <- consensus[row, ]
      cli::cli_alert_success(
        "Optimal BC{rec$barcode_ID} flank: {rec$orientation} reads; {rec$direction}; size = {rec$flank_size}."
      )
      cli::cli_alert_info("Proportion of retained reads with single flank hit: {rec$hits} / {total_reads}.")
    }
  }

  # -- Optional file saving (only if output_dir is specified) --

  if (!is.null(output_dir)) {
    identity_barcode_info <- cbind(identity_barcode_info,
                                   consensus |> dplyr::select(flank_size, orientation, direction) |>
                                     dplyr::rename(optimal_flank_size = flank_size))

    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    readr::write_csv(identity_barcode_info, file = file.path(output_dir, identity_barcode_info_file_name))

    cli::cli_alert_success("Identity barcode information saved to: {normalizePath(output_dir)}.")
  }
}


#' @keywords internal
test_flanks <- function(
  reads,
  identity_barcode_info,
  direction,
  min_size,
  max_size
) {
  # Extract barcode info
  amplicon_seq <- identity_barcode_info$amplicon_seq[1]
  start        <- identity_barcode_info$start[1]
  end          <- identity_barcode_info$end[1]
  flank_sizes  <- min_size:max_size

  # One vectorized str_count() pass per flank size against ALL reads at once,
  # collected into a match-count matrix (n_reads x n_flank_sizes).
  match_counts <- vapply(
    flank_sizes,
    function(i) {
      test_flank <- if (direction == "upstream") {
        stringr::str_sub(amplicon_seq, start - i, start - 1)
      } else if (direction == "downstream") {
        stringr::str_sub(amplicon_seq, end + 1, end + i)
      } else {
        cli::cli_abort("{direction}: Invalid value for direction.")
      }
      stringr::str_count(reads$seq, stringr::fixed(test_flank))
    },
    FUN.VALUE = integer(nrow(reads))
  )

  # Hit / multiple-hit counts come from colSums() on above matrix
  return(list(
    hits_vec          = as.integer(colSums(match_counts >= 1)),
    multiple_hits_vec = as.integer(colSums(match_counts > 1))
  ))
}


#' @keywords internal
export_flank_test_stats <- function(
  ID,
  all_stats,
  ordered_cols,
  sample_name,
  output_dir
) {
  wide_stats <- all_stats |>
    # Reshape flank stats for user-friendly output
    dplyr::filter(ID == barcode_ID) |>
    tidyr::pivot_wider(
      id_cols     = flank_size,
      names_from  = c(orientation, direction),
      names_glue  = "{orientation}_{direction}_{.value}",
      values_from = c(hits, multiple_hits)
    ) |>
    dplyr::select(flank_size, dplyr::all_of(ordered_cols))

  BCX_stats_file <- file.path(output_dir, paste0(sample_name, "_BC", ID, "_flank_stats.csv"))
  readr::write_csv(wide_stats, file = BCX_stats_file)

  cli::cli_alert_success("BC{ID} flank stats exported to: {normalizePath(BCX_stats_file)}")
}



