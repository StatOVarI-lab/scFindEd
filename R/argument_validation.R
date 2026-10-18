# -- Argument Validation Helpers --

#' Used to ensure an argument is an integer
#' @keywords internal
assert_integer <- function(
  argument,
  argument_name,
  min = 1,
  max = Inf
) {
  if (!is.numeric(argument) ||        # Must be a number
      length(argument) != 1 ||        # Must be a single value
      is.na(argument) ||              # Must not be NA
      argument != round(argument)  || # Must be a whole number
      argument < min ||               # Must be greater than or equal to min (1)
      argument > max                  # Must be less    than or equal to max (infinity)
  ) {
    cli::cli_abort("{argument_name} must be an integer between {min} and {max}.")
  }
}


#' Used to ensure an argument is a float
#' @keywords internal
assert_float <- function(
  argument,
  argument_name,
  min = 0,
  max = 1,
  inclusive_max = TRUE
) {
  # Single check for if the value exceeds the max value - inclusive or not
  exceeds_max <- if (inclusive_max) argument > max else argument >= max

  if (!is.numeric(argument) || # Must be a number
      length(argument) != 1 || # Must be a single value
      is.na(argument) ||       # Must not be NA
      argument < min ||        # Must be greater than or equal to min (0)
      exceeds_max              # Must be less    than or equal to max (1)
  ) {
    cli::cli_abort("{argument_name} must be a number between {min} and {max}.")
  }
}


#' Used to ensure an argument is a boolean
#' @keywords internal
assert_boolean <- function(
  argument,
  argument_name
) {
  if (!is.logical(argument) || # Must be a boolean
      length(argument) != 1 || # Must be a single value
      is.na(argument)          # Must not be NA
  ) {
    cli::cli_abort("{argument_name} must be Boolean: TRUE or FALSE.")
  }
}


#' Used to ensure an argument is a string
#' @keywords internal
assert_string <- function(
  argument,
  argument_name
) {
  if (!is.character(argument) || # Must be a character string
      length(argument) != 1   || # Must be a single value
      is.na(argument)            # Must not be NA
  ) {
    cli::cli_abort("{argument_name~ must be a single, valid character string.")
  }
}


# -- Per-Function Argument Validation --

#' @keywords internal
validate_identify_true_mutations_args <- function(
  source_dir,
  sample_name,
  panel_ID,
  min_reads = 10,
  true_mutation_proportion = 0.1,
  barcoded = FALSE,
  min_cells = NULL,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Individual argument checks --
  assert_string(source_dir, "source_dir")
  assert_string(sample_name, "sample_name")
  assert_string(panel_ID, "panel_ID")
  assert_integer(min_reads, "min_reads", min = 3)
  assert_float(true_mutation_proportion, "true_mutation_proportion", inclusive_max = FALSE)
  assert_boolean(barcoded, "barcoded")
  if (!is.null(min_cells)) assert_integer(min_cells, "min_cells")
  if (!is.null(output_dir)) assert_string(output_dir, "output_dir")
  assert_boolean(verbose, "verbose")

  invisible(TRUE)
}


#' @keywords internal
validate_genotype_cells_args<- function(
  source_dir,
  sample_name,
  panel_ID,
  alpha_0,
  ploidy = "haploid",
  credible_interval = 0.95,
  lower_CI_assignment_threshold = NULL,
  output_dir = NULL,
  plots_output_dir = NULL,
  verbose = TRUE
) {
  # -- Individual argument checks --
  assert_string(source_dir, "source_dir")
  assert_string(sample_name, "sample_name")
  assert_string(panel_ID, "panel_ID")
  assert_string(ploidy, "ploidy")
  assert_float(credible_interval, "credible_interval", inclusive_max = FALSE)
  assert_float(lower_CI_assignment_threshold, "lower_CI_assignment_threshold", inclusive_max = FALSE)
  if (!is.null(min_cells)) assert_integer(min_cells, "min_cells")
  if (!is.null(output_dir)) assert_string(output_dir, "output_dir")
  assert_boolean(verbose, "verbose")

  invisible(TRUE)
}


#' @keywords internal
validate_call_identity_barcodes_args <- function(
  source_dir,
  sample_name,
  panel_ID,
  ploidy = "haploid",
  lower_CI_assignment_threshold = 0.7,
  consistency_threshold = 0.5,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Individual argument checks --
  assert_string(source_dir, "source_dir")
  assert_string(sample_name, "sample_name")
  assert_string(panel_ID, "panel_ID")
  assert_string(ploidy, "ploidy")
  assert_float(lower_CI_assignment_threshold, "lower_CI_assignment_threshold", inclusive_max = FALSE)
  assert_float(consistency_threshold, "consistency_threshold", inclusive_max = FALSE)
  if (!is.null(output_dir)) assert_string(output_dir, "output_dir")
  assert_boolean(verbose, "verbose")

  invisible(TRUE)
}


#' @keywords internal
validate_generate_synthetic_experiment_args <- function(
  total_reads,
  reads_per_amplicon = 100000,
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
  # -- Individual argument checks --

  assert_integer(reads_per_amplicon, "reads_per_amplicon", min = 10000)
  assert_integer(total_cell_barcodes, "total_cell_barcodes", min = 1000)
  assert_integer(total_whitelist_barcodes, "total_whitelist_barcodes")
  assert_integer(total_identity_barcode_combinations, "total_identity_barcode_combinations", min = 100)
  assert_float(soup_fraction, "soup_fraction", inclusive_max = FALSE)
  assert_integer(n_amplicons, "n_amplicons", max = 50)
  assert_integer(n_edits, "n_edits", max = 100)
  assert_integer(total_identity_barcodes, "total_identity_barcodes", max = 5)
  assert_integer(amplicon_length, "amplicon_length", min = 225, max = 500)
  assert_integer(paired_read_length, "paired_read_length", min = 120)
  assert_integer(identity_barcode_length, "identity_barcode_length", min = 10, max = 40)
  assert_integer(random_seed, "random_seed", min = -Inf)
  assert_boolean(verbose, "verbose")
  if (!is.null(output_dir)) assert_string(output_dir, "output_dir")

  # -- Combinatorial / Specific argument checks --

  if (total_cell_barcodes <= total_identity_barcode_combinations) {
    cli::cli_abort("total_cell_barcodes must exceed total_identity_barcode_combinations.")
  }
  if (total_reads <= total_cell_barcodes) {
    cli::cli_abort("Total reads: reads_per_amplicon * (n_amplicons + 1), must exceed total_cell_barcodes.")
  }
  if (amplicon_length <= paired_read_length) {
    cli::cli_abort("paired_read_length must exceed amplicon_length.")
  }
  # Ensure specific identity barcodes can fit within the identity amplicon's length
  necessary_length <- 20 + (total_identity_barcodes * (identity_barcode_length + 20))
  if (amplicon_length < necessary_length) {
    cli::cli_abort("{total_identity_barcodes} identity barcodes of length {identity_barcode_length} require an amplicon_length of at least {necessary_length}.")
  }

  invisible(TRUE)
}


