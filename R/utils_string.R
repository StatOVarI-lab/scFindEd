#' @keywords internal
#' Imports data from a file
#' Expected file name: "cell_to_extracted_BCs.rds"
import_cell_to_extracted_BCs <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- readRDS(file_path),
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
#' Expected file name: "barcode_list_name.csv"
import_BC1_whitelist <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  # Load BC1 whitelist from CSV and assign each a unique name
  tryCatch(
    data <- readr::read_csv(file_path, col_names = FALSE,
                            show_col_types = FALSE) %>%
      dplyr::pull(X1),
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
#' Expected file name: "barcode_list_name.csv"
import_BC2_BC3_whitelist <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  # Load barcode pair whitelist
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


# Function to compute Hamming distances between a barcode and a whitelist,
# returning the 3 smallest distances
identify_barcode_ID <- function(BC, BC_vec){
  temp <- sapply(BC_vec,
                 function(x) stringdist::stringdist(x, BC, method = "hamming"))
  return(sort(temp)[1:3])
}
