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
#' Expected file name: "list_BCX.rds"
import_barcode_distance_list <- function(
    source_dir,
    file_name,
    verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  # Import list of BC distances
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
#' Corrects a barcode by selecting the closest match, as long as the closest
#' match has its distance below a defined threshold.
correct_BC <- function(distances, thresh){
  # If the distance < the threshold, select the barcode as a match
  if (distances[1] < thresh){
    output <- names(distances)[1]
  }
  # Otherwise, return no output
  else {
    output <- ""
  }
  return(output)
}
