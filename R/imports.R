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
    data <- readRDS(file_path),
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )
  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from files
#' Expected file name(s): "list_BC1.rds", "list_BC2_3.rds"
#' @keywords internal
import_barcode_distance_lists <- function(
  source_dir,
  list_BC1_file_name,
  list_BC2_3_file_name,
  verbose
) {
  # -- Importing BC1 list --
  list_BC1 <- import_barcode_distance_list(
    source_dir,
    list_BC1_file_name,
    verbose
  )

  # -- Importing BC2/3 list --
  list_BC2_3 <- import_barcode_distance_list(
    source_dir,
    list_BC2_3_file_name,
    verbose
  )

  return(list(
    list_BC1 = list_BC1,
    list_BC2_3 = list_BC2_3
  ))
}


#' Imports data from a file
#' Expected file name: "list_BCX.rds"
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
#' Expected file name: "barcode_list_name.csv"
#' @keywords internal
import_BC1_whitelist <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load BC1 whitelist from CSV and assign each a unique name
  tryCatch(
    data <- readr::read_csv(file_path, col_names = FALSE,
                            show_col_types = FALSE) |>
      dplyr::pull(X1),
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_process_done()

  return(data)
}


#' Imports data from a file
#' Expected file name: "barcode_list_name.csv"
#' @keywords internal
import_BC2_BC3_whitelist <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_process_start(msg = "Importing {file_path}...",
                                      msg_done = "Successfully imported {file_name}.")

  # Load barcode pair whitelist
  tryCatch(
    data <- readr::read_csv(file_path, show_col_types = FALSE),
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
