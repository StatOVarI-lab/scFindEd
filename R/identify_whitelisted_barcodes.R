#' Identifying Closest Whitelisted Barcodes
#'
#' This function identifies the three shortest Hamming distances between barcode
#' and whitelisted barcode, for each of barcode 1 and barcode 2-3.
#'
#' @inheritParams global_params source_dir sample_name output_dir verbose
#'
#' @param BC1_whitelist_file_name String containing the exact file name (including
#' the file type) for the barcode 1 whitelist.
#' @param BC2_BC3_whitelist_file_name String containing the exact file name
#' (including the file type) for the barcode 2-3 whitelist.
#'
#' @return A list of barcode matches for barcode 1 and a list of barcode matches
#' for barcode 2-3.
#'
#' @export
identify_whitelisted_barcodes <- function(
  source_dir,
  sample_name,
  BC1_whitelist_file_name,
  BC2_BC3_whitelist_file_name,
  output_dir = NULL,
  verbose = TRUE
) {
  # -- Import whitelists --

  imported_whitelists <- import_whitelisted_barcodes(source_dir,
                                                     BC1_whitelist_file_name,
                                                     BC2_BC3_whitelist_file_name,
                                                     verbose)

  invisible(list2env(imported_whitelists, envir = environment()))

  # -- Import cell to extracted BCs --

  cell_to_extracted_BCs_file_name <- paste0(sample_name, "_cell_to_extracted_BCs.csv")
  cell_to_extracted_BCs <- import_cell_to_extracted_BCs(source_dir,
                                                        cell_to_extracted_BCs_file_name,
                                                        verbose)

  # -- Identify the Three BC1s with the Shortest Hamming Distance

  list_BC1 <- identify_BC1_closest_barcodes(cell_to_extracted_BCs,
                                            BC1_whitelist,
                                            verbose)

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(list_BC1, file = file.path(output_dir, "list_BC1.rds"))

    cli::cli_alert_success(
      "BC1: List of closest barcodes saved to: {normalizePath(output_dir)}.")
  }

  # -- Identify the Three BC2_3 pairs with the Shortest Hamming Distance

  list_BC2_3 <- identify_BC2_3_closest_barcodes(
    cell_to_extracted_BCs,
    BC2_BC3_whitelist,
    verbose
  )

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    saveRDS(list_BC2_3, file = file.path(output_dir, "list_BC2_3.rds"))

    cli::cli_alert_success(
      "BC2_3: List of closest barcode pairs saved to: {normalizePath(output_dir)}.")
  }
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
  cli::cli_h2("Processing BC1 whitelist")

  # Get unique BC1 sequences from the dataset
  unique_BC_1 <- unique(cell_to_extracted_BCs$BC1)

  # Extract the complete list of barcode names
  whitelist_names <- names(BC1_whitelist)

  if (verbose) cli::cli_alert_info("Number of Unique BC1s: {length(unique_BC_1)}")
  cli::cli_process_start("Computing BC1 distances and extracting three closest...",
                         msg_done = "Three closest BC1s extracted.")
  chunk_size <- 1000L
  n_query <- length(unique_BC_1)
  list_BC1 <- vector("list", n_query)

  for (start in seq(1L, n_query, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n_query)
    idx <- start:end

    dist_block <- stringdist::stringdistmatrix(
      unique_BC_1[idx],
      BC1_whitelist,
      method = "hamming",
      useBytes = TRUE
    )

    # Extract the top 3 from each row of this block
    list_BC1[idx] <- lapply(seq_len(nrow(dist_block)), extract_closest_barcodes,
                            dist_matrix = dist_block,
                            whitelist_names = whitelist_names,
                            n_barcodes = 3)
  }

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
  cli::cli_h2("Processing BC2_3 whitelist")

  # Get unique BC2 and BC3 sequences from the dataset
  cell_BC1_BC2_BC3 <- cell_to_extracted_BCs
  unique_BC_2_3 <- unique(c(cell_BC1_BC2_BC3$BC2, cell_BC1_BC2_BC3$BC3))

  # Extract the complete list of barcode names
  whitelist_names <- names(BC2_BC3_whitelist)

  if (verbose) {
    cli::cli_alert_info("Number of Unique BC2_3 pairs: {length(unique_BC_2_3)}")
  }

  cli::cli_process_start("Extracting three closest BC2_3 pairs...",
                         msg_done = "Three closest BC2_3 pairs extracted.")

  chunk_size <- 1000L
  n_query <- length(unique_BC_2_3)
  list_BC2_3 <- vector("list", n_query)

  for (start in seq(1L, n_query, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n_query)
    idx <- start:end

    dist_block <- stringdist::stringdistmatrix(
      unique_BC_2_3[idx],
      BC2_BC3_whitelist,
      method = "hamming",
      useBytes = TRUE
    )

    # Extract the top 3 from each row of this block
    list_BC2_3[idx] <- lapply(seq_len(nrow(dist_block)), extract_closest_barcodes,
                              dist_matrix = dist_block,
                              whitelist_names = whitelist_names,
                              n_barcodes = 3)
  }

  # Assign BC2_3 sequences as names to the list entries
  names(list_BC2_3) <- unique_BC_2_3

  cli::cli_process_done()

  return(list_BC2_3)
}
