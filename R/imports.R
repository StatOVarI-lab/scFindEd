#' Imports data from a file
#' Expected file name: "ref_amplicons_sorted.fa"
#' @keywords internal
import_reference_amplicons <- function(
  source_dir,
  verbose
) {
  file_name <- "ref_amplicons_sorted.fa"
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load experiment design summary file
  tryCatch(
    if (verbose) {
      data <- Biostrings::readDNAStringSet(file_path)
    } else {
      data <- suppressMessages(Biostrings::readDNAStringSet(file_path))
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "\{sample_name\}_UTR_barcode_info.csv"
#' @keywords internal
import_UTR_barcode_info <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load experiment design summary file
  tryCatch(
    if (verbose) {
      data <- readr::read_csv(file_path, show_col_types = FALSE)
    } else {
      data <- suppressMessages(readr::read_csv(file_path, show_col_types = FALSE))
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "cell_to_extracted_BCs.rds"
#' @keywords internal
import_cell_to_extracted_BCs <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  tryCatch(
    if (verbose) {
      data <- readr::read_csv(file_path, show_col_types = FALSE)
    } else {
      data <- suppressMessages(readr::read_csv(file_path, show_col_types = FALSE))
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )
  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from files
#' Expected file name(s): "\{sample_name\}_\{BC_IDs\}_closest_barcodes.rds"
#' @keywords internal
import_barcode_distance_lists <- function(
  source_dir,
  sample_name,
  BC_IDs,
  verbose
) {
  # Initialise empty list for dynamically storing barcode distance lists
  barcode_distance_lists <- list()

  for (BC_ID in BC_IDs) {
    BC_ID_file_name = paste0(sample_name, "_", BC_ID,
                                           "_closest_barcodes.rds")

    barcode_distance_lists[[BC_ID]] <- import_barcode_distance_list(source_dir,
                                                                    BC_ID_file_name,
                                                                    verbose)
  }

  return(barcode_distance_lists)
}


#' Imports data from a file
#' Expected file name: "\{sample_name\}_\{BC_ID\}_closest_barcodes.rds"
#' @keywords internal
import_barcode_distance_list <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Import list of BC distances
  tryCatch(
    data <- readRDS(file_path),
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "BCX_whitelist_barcodes.csv"
#' @keywords internal
import_barcode_whitelist <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load barcode whitelist from CSV and assign each a unique name
  tryCatch(
    if (verbose) {
      data <- readr::read_csv(file_path, col_names = FALSE, show_col_types = FALSE) |>
        dplyr::pull(X1)
    } else {
      data <- suppressMessages(readr::read_csv(file_path, col_names = FALSE, show_col_types = FALSE)) |>
        dplyr::pull(X1)
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "\{panel_ID\}_design_summary_joined.csv"
#' @keywords internal
import_design_summary_info <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load experiment design summary file
  tryCatch(
    if (verbose) {
      data <- readr::read_csv(file_path, show_col_types = FALSE)
    } else {
      data <- suppressMessages(readr::read_csv(file_path, show_col_types = FALSE))
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "cell_to_UTR_barcode_system.csv"
#' @keywords internal
import_cell_to_UTR_barcode_system <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  tryCatch(
    data <- readr::read_csv(file_path, show_col_types = FALSE),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: e.g. "\{sample_name\}_\{UTR_barcode_AmpID\}_data.tsv.gz"
#' @keywords internal
import_UTR_barcode_reads <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  tryCatch(
    data <- readr::read_tsv(file_path, show_col_types = FALSE),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "\{sample_name\}_\{AmpID\}_fragment_mutations.tsv"
#' @keywords internal
import_fragment_mutations <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load experiment design summary file
  tryCatch(
    if (verbose) {
      data <- readr::read_tsv(file_path, show_col_types = FALSE)
    } else {
      data <- suppressMessages(readr::read_tsv(file_path, show_col_types = FALSE))
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "\{sample_name\}_\{AmpID\}_variant_summary.tsv"
#' @keywords internal
import_variant_summary <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load experiment design summary file
  tryCatch(
    if (verbose) {
      data <- readr::read_tsv(file_path, show_col_types = FALSE)
    } else {
      data <- suppressMessages(readr::read_tsv(file_path, show_col_types = FALSE))
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "\{sample_name\}_\{AmpID\}_true_mutations.tsv"
#' @keywords internal
import_true_mutations <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load experiment design summary file
  tryCatch(
    if (verbose) {
      data <- readr::read_tsv(file_path, show_col_types = FALSE)
    } else {
      data <- suppressMessages(readr::read_tsv(file_path, show_col_types = FALSE))
    },
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}
