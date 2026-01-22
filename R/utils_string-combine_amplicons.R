#' @keywords internal
#' Imports data from a file
#' Expected file name: e.g. "final_summary_minimum_2counts_70_30_70_all_amplicons.csv"
import_combined_summaries_data <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- readr::read_csv(file_path, show_col_types = FALSE),
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
  data <- tibble::as_tibble(data, show_col_types = FALSE)

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}
