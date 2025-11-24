#' @keywords internal
# Function to compute Hamming distances between a barcode and a whitelist,
# returning the 3 smallest distances
identify_barcode_ID <- function(BC, BC_vec){
  temp <- sapply(BC_vec,
                 function(x) stringdist::stringdist(x, BC, method = "hamming"))
  return(sort(temp)[1:3])
}


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
    "Successfully imported BC1 whitelist.")

  return(cell_to_extracted_BCs)
}
