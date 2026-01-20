#' @keywords internal
#' Imports data from a file
#' Expected file name: "HDR_file_barcoded.csv"
import_HDR_barcoded_data <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- readr::read_csv(file_path, show_col_types = FALSE),
    error = function(e) {
      cli::cli_abort("Failed to import {file_name} from {source_dir}.")
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}


#' @keywords internal
#' Imports data from a file
#' Expected file name: "6915_design_summary_extended_information.csv"
import_design_summary_info <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

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

  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")

  return(data)
}


#' @keywords internal
#' Takes a list of data frame names with vectors of column names to check for in
#' each, and returns TRUE/FALSE dependent on whether every required column is
#' present.
check_required_columns <- function(
  required_cols
) {
  columns_present = TRUE

  for (df_name in names(required_cols)) {
    check_df <- get(df_name, envir = parent.frame())
    cols_to_check <- required_cols[[df_name]]

    # Finds columns present in required_cols but not in check_df
    missing_cols <- setdiff(cols_to_check, colnames(check_df))

    # Update the user if columns are missing
    if(length(missing_cols) > 0) {
      cli::cli_abort(
        "Missing column(s) in {df_name}: {paste(missing_cols, collapse = ', ')}")
      columns_present = FALSE
    }
  }
  return(columns_present)
}
