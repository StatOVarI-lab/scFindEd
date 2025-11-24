#' Identifying Closest Whitelisted Barcodes
#'
#' This function identifies the three shortest Hamming distances between barcode
#' and whitelisted barcode, for each of BC1 and BC2/BC3.
#' @param source_dir String containing the relative path to the directory that
#' source data should be imported from.
#' @param BC1_whitelist_file_name String containing the file name structure for
#' the BC1 barcode whitelist.
#' @param BC2_BC3_whitelist_file_name String containing the file name structure
#' for the BC2/3 barcode whitelist.
#' @param cell_to_extracted_BCs_file_name A string containing the file name
#' structure for the RDS data frame linking cells to extracted barcodes.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#' @param verbose Boolean deciding whether detailed outputs should be printed.

#' @return A list of barcode matches for BC1 and a list of barcode matches for
#' BC2/BC3.

#' @examples
#' identify_whitelisted_barcodes(
#' source_dir = "data",
#' BC1_whitelist_file_name = "75K_barcode.csv",
#' BC2_BC3_whitelist_file_name = "Wong_barcode_pairs.csv",
#' cell_to_extracted_BCs_file_name = "cell_to_extracted_BCs.rds",
#' output_dir = "data",
#' verbose = TRUE
#' )

#' @export
identify_whitelisted_barcodes <- function(
  source_dir,
  BC1_whitelist_file_name,
  BC2_BC3_whitelist_file_name,
  cell_to_extracted_BCs_file_name,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import whitelists --

  imported_whitelists <- import_whitelisted_barcodes(
    source_dir,
    BC1_whitelist_file_name,
    BC2_BC3_whitelist_file_name,
    verbose
  )

  invisible(list2env(imported_whitelists, envir = environment()))

  # -- Import cell to extracted BCs --

  cell_to_extracted_BCs <- import_cell_to_extracted_BCs(
    source_dir,
    cell_to_extracted_BCs_file_name,
    verbose
  )

  cli::cli_alert_success("Successfully imported datasets.")

  # -- Identify the Three BC1s with the Shortest Hamming Distance

  list_BC1 <- identify_BC1_closest_barcodes(
    cell_to_extracted_BCs,
    BC1_whitelist,
    verbose
  )

  # -- Identify the Three BC2/BC3s with the Shortest Hamming Distance

  list_BC2_3 <- identify_BC2_3_closest_barcodes(
    cell_to_extracted_BCs,
    BC2_BC3_whitelist,
    verbose
  )

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(list_BC1, file = file.path(output_dir, "list_BC1.rds"))

    cli::cli_alert_success(
      "BC1: List of closest barcodes saved to: {normalizePath(output_dir)}.")

    saveRDS(list_BC2_3, file = file.path(output_dir, "list_BC2_3.rds"))

    cli::cli_alert_success(
      "BC2/3: List of closest barcodes saved to: {normalizePath(output_dir)}.")
  }

  invisible(list(list_BC1 = list_BC1, list_BC1 = list_BC2_3))
}


#' @keywords internal
import_whitelisted_barcodes <- function(
  source_dir,
  BC1_whitelist_file_name,
  BC2_BC3_whitelist_file_name,
  verbose
) {
  if (verbose) cli::cli_inform("Importing barcode whitelist files...")

  # Load BC1 whitelist from CSV and assign each a unique name
  file_path <- file.path(source_dir, BC1_whitelist_file_name)
  tryCatch(
    BC1_whitelist <- utils::read.table(file_path,sep = ",")$V1,
    error = function(e) {
      cli::cli_abort("Failed to import BC1 whitelist from {.file {file_path}}.")
    }
  )
  names(BC1_whitelist) <- paste0("BC1_", 1:length(BC1_whitelist))

  if (verbose) cli::cli_alert_success(
    "Successfully imported BC1 whitelist.")

  # Load barcode pair whitelist and extract unique BC2/BC3 sequences
  file_path <- file.path(source_dir, BC2_BC3_whitelist_file_name)
  tryCatch(
    BC2_BC3_whitelist_raw <- utils::read.table(file_path, sep = ",", header = TRUE),
    error = function(e) {
      cli::cli_abort("Failed to import BC2/3 whitelist from {.file {file_path}}.")
    }
  )
  # Extract the same strings for BC2 and BC3
  BC2_BC3_whitelist <- unique(BC2_BC3_whitelist_raw[, 1])
  names(BC2_BC3_whitelist) <- paste0("BC2_3_", 1:length(BC2_BC3_whitelist))

  if (verbose) cli::cli_alert_success(
    "Successfully imported BC2/3 whitelist.")

  return(list(BC1_whitelist = BC1_whitelist,
              BC2_BC3_whitelist = BC2_BC3_whitelist))
}


#' @keywords internal
identify_BC1_closest_barcodes <- function(
  cell_to_extracted_BCs,
  BC1_whitelist,
  verbose
) {
  if (verbose) cli::cli_inform("Processing BC1 whitelist...")

  # Get unique BC1 sequences from the dataset
  cell_BC1_BC2_BC3 <- cell_to_extracted_BCs
  unique_BC_1 <- unique(cell_BC1_BC2_BC3$BC1)

  if (verbose) cli::cli_alert_info("Number of Unique BC1s: {length(unique_BC_1)}")

  # Loop through each unique BC1 sequence and find the 3 closest whitelist
  # barcodes to the current BC1 using Hamming distance
  cli::cli_progress_bar("BC1: Identifying closest barcodes",
                   total = length(unique_BC_1))
  list_BC1 <- list()
  for (j in 1:length(unique_BC_1)){
    list_BC1[[j]] <- identify_barcode_ID(unique_BC_1[j], BC1_whitelist)
    cli::cli_progress_update() # Update progress
  }
  cli::cli_progress_done() # Close the progress bar

  # Assign BC1 sequences as names to the list entries
  names(list_BC1) <- unique_BC_1

  cli::cli_alert_success("Successfully identified closest BC1 barcodes.")

  return(list_BC1)
}


#' @keywords internal
identify_BC2_3_closest_barcodes <- function(
  cell_to_extracted_BCs,
  BC2_BC3_whitelist,
  verbose
) {
  if (verbose) cli::cli_inform("Processing BC2/3 whitelist...")

  # Get unique BC2 and BC3 sequences from the dataset
  cell_BC1_BC2_BC3 <- cell_to_extracted_BCs
  unique_BC_2_3 <- unique(c(cell_BC1_BC2_BC3$BC2, cell_BC1_BC2_BC3$BC3))

  if (verbose) {
    cli::cli_alert_info("Number of Unique BC2_3s: {length(unique_BC_2_3)}")
  }

  # Loop through each unique BC2/BC3 sequence and find the 3 closest whitelist
  # barcodes to the current BC2/BC3 using Hamming distance
  cli::cli_progress_bar("BC2/3: Identifying closest barcodes",
                   total = length(unique_BC_2_3))
  list_BC2_3 <- list()
  for (j in 1:length(unique_BC_2_3)){
    list_BC2_3[[j]] <- identify_barcode_ID(unique_BC_2_3[j], BC2_BC3_whitelist)
    cli::cli_progress_update() # Update progress
  }
  cli::cli_progress_done() # Close the progress bar

  # Assign BC2/BC3 sequences as names to the list entries
  names(list_BC2_3) <- unique_BC_2_3

  cli::cli_alert_success("Successfully identified closest BC2/3 barcodes.")

  return(list_BC2_3)
}
