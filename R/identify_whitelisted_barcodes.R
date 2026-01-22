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

#' @importFrom dplyr %>%

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

  if (verbose) cli::cli_inform("Exporting lists of closest barcodes...")

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
  BC1_whitelist <- import_BC1_whitelist(
    source_dir,
    BC1_whitelist_file_name,
    verbose
  )

  names(BC1_whitelist) <- paste0("BC1_", 1:length(BC1_whitelist))

  BC2_BC3_whitelist_raw <- import_BC2_BC3_whitelist(
    source_dir,
    BC2_BC3_whitelist_file_name,
    verbose
  )

  # Extract the same strings for BC2 and BC3
  BC2_BC3_whitelist <- as.character(unique(BC2_BC3_whitelist_raw[, 1])[[1]])
  names(BC2_BC3_whitelist) <- paste0("BC2_3_", 1:length(BC2_BC3_whitelist))

  return(list(BC1_whitelist = BC1_whitelist,
              BC2_BC3_whitelist = BC2_BC3_whitelist))
}

#' @keywords internal
identify_BC1_closest_barcodes <- function(
  cell_to_extracted_BCs,
  BC1_whitelist,
  verbose
) {
  cli::cli_inform("Processing BC1 whitelist...")

  # Get unique BC1 sequences from the dataset
  unique_BC_1 <- unique(cell_to_extracted_BCs$BC1)

  # Extract the complete list of barcode names
  whitelist_names <- names(BC1_whitelist)

  if (verbose) cli::cli_alert_info("Number of Unique BC1s: {length(unique_BC_1)}")

  cli::cli_process_start("Computing distance matrix...",
                         msg_done = "Matrix computed.")

  # Calculate distance matrix between each extracted barcode and whitelist barcode
  dist_matrix <- stringdist::stringdistmatrix(
    unique_BC_1,
    BC1_whitelist,
    method = "hamming",
    useBytes = TRUE
  )

  cli::cli_process_done()

  cli::cli_process_start("Extracting three closest BC1s...",
                         msg_done = "Successfully extracted three closest BC1s.")

  # Use lapply to extract the top 3 from each row of the matrix
  list_BC1 <- lapply(1:nrow(dist_matrix), extract_closest_barcodes,
                     dist_matrix = dist_matrix,
                     whitelist_names = whitelist_names,
                     n_barcodes = 3)

  # Assign BC1 sequences as names to the list entries
  names(list_BC1) <- unique_BC_1

  cli::cli_process_done()

  return(list_BC1)
}


#' @keywords internal
identify_BC2_3_closest_barcodes <- function(
    cell_to_extracted_BCs,
    BC2_BC3_whitelist,
    verbose
) {
  cli::cli_inform("Processing BC2/3 whitelist...")

  # Get unique BC2 and BC3 sequences from the dataset
  cell_BC1_BC2_BC3 <- cell_to_extracted_BCs
  unique_BC_2_3 <- unique(c(cell_BC1_BC2_BC3$BC2, cell_BC1_BC2_BC3$BC3))

  # Extract the complete list of barcode names
  whitelist_names <- names(BC2_BC3_whitelist)

  if (verbose) {
    cli::cli_alert_info("Number of Unique BC2_3s: {length(unique_BC_2_3)}")
  }

  cli::cli_process_start("Computing distance matrix...",
                         msg_done = "Matrix computed.")

  # Calculate distance matrix between each extracted barcode and whitelist barcode
  dist_matrix <- stringdist::stringdistmatrix(
    unique_BC_2_3,
    BC2_BC3_whitelist,
    method = "hamming",
    useBytes = TRUE
  )

  cli::cli_process_done()

  cli::cli_process_start("Extracting three closest BC2/3s...",
                         msg_done = "Successfully extracted three closest BC2/3s.")

  # Use lapply to extract the top 3 from each row of the matrix
  list_BC2_3 <- lapply(1:nrow(dist_matrix), extract_closest_barcodes,
                       dist_matrix = dist_matrix,
                       whitelist_names = whitelist_names,
                       n_barcodes = 3)

  # Assign BC1 sequences as names to the list entries
  names(list_BC2_3) <- unique_BC_2_3

  cli::cli_process_done()

  return(list_BC2_3)
}
