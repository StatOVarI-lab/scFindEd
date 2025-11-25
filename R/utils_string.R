#' @keywords internal
import_cell_to_extracted_BCs <- function(
  source_dir,
  cell_to_extracted_BCs_file_name,
  verbose
) {
  if (verbose) cli::cli_inform("Importing cell to extracted barcodes table...")

  tryCatch(
    cell_to_extracted_BCs <- readRDS(
      file.path(source_dir, cell_to_extracted_BCs_file_name)),
    error = function(e) {
      cli::cli_abort(
        "Failed to import cell to extract BCs data table from {.file {file_path}}."
      )
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported cell to extracted barcodes table.")

  return(cell_to_extracted_BCs)
}


#' @keywords internal
import_barcode_distance_lists <- function(
    source_dir,
    list_BC1_file_name,
    list_BC2_3_file_name,
    verbose
) {
  # -- Importing BC1 list --
  if (verbose) cli::cli_inform("Importing list of BC1 distances...")

  # Import list of BC1 distances
  file_path <- file.path(source_dir, list_BC1_file_name)
  tryCatch(
    list_BC1 <- readRDS(file_path),
    error = function(e) {
      cli::cli_abort(
        "Failed to import list of BC1 distances from {.file {file_path}}."
      )
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported list of BC1 distances.")

  # -- Importing BC2/3 list --
  if (verbose) cli::cli_inform("Importing list of BC2/3 distances...")

  # Import list of BC2/3 distances
  file_path <- file.path(source_dir, list_BC2_3_file_name)
  tryCatch(
    list_BC2_3 <- readRDS(file_path),
    error = function(e) {
      cli::cli_abort(
        "Failed to import list of BC2/3 distances from {.file {file_path}}."
      )
    }
  )

  if (verbose) cli::cli_alert_success(
    "Successfully imported list of BC2/3 distances.\n")

  return(list(
    list_BC1 = list_BC1,
    list_BC2_3 = list_BC2_3
  ))
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
