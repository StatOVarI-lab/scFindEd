#' Cell-Barcode Assignment
#'
#' This function begins by correcting barcode assignments by applying Hamming
#' distance-based threshold correction. It then assigns cell-barcode combinations
#' as follows: a barcode combination is assigned to a cell if and only if there
#' are at least `min_reads` reads of that barcode combination in the cell, and
#' they constitute more than `min_reads_proportion` of the barcode reads in the
#' cell.
#'
#' @param cell_to_extracted_BCs Data frame linking cells to extracted barcodes.
#' @param list_BC1 Named list of nearest BC1 barcode candidates.
#' @param list_BC2_3 Named list of nearest BC2/BC3 barcode candidates.
#' @param BC_thresholds List containing the following barcode thresholds:
#' `BC1_threshold`, `BC2_3_threshold`.
#' @param reads_thresholds List containing the following barcode correction
#' thresholds: `min_reads`, `min_reads_proportion`.
#' @param output_dir String containing the relative path to the directory that
#' the output data should be exported to.
#'
#' @return A filtered data frame of barcode assignments that meet the read count
#' and proportion thresholds.

#' @examples
#' correct_barcodes(
#' cell_to_extracted_BCs = cell_to_extracted_BCs,
#' list_BC1 = list_BC1,
#' list_BC2_3 = list_BC2_3,
#' BC_thresholds = list(BC1_threshold = 3, BC2_3_threshold = 5),
#' reads_thresholds = list(min_reads = 5, min_reads_proportion = 0.9),
#' output_dir = "data"
#' )

#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom dplyr %>% filter mutate group_by summarise
#' @importFrom glue glue

#' @export
cell_barcode_assignment <- function(
    cell_to_extracted_BCs,
    list_BC1,
    list_BC2_3,
    BC_thresholds,
    reads_thresholds,
    output_dir
) {
  # -- Barcode Correction --

  message("Performing Barcode Correction...")

  cell_to_BCs_corrected <- correct_barcodes(
    cell_to_extracted_BCs,
    list_BC1,
    list_BC2_3,
    BC_thresholds
  )

  message("Barcodes Corrected.")

  # -- Cell Barcode Assignment --

  message("Performing Cell Barcode Assignment...")

  # Creates column to store the proportion of reads per barcode within each cell
  cell_to_BCs_corrected$proportion_reads_in_cell <- NA

  # Loops through each unique cell barcode  to identify rows corresponding to the
  # current cell barcode, and calculate the proportion of reads for each barcode
  # assignment within the cell
  cli_progress_bar("Assigning Barcodes", total = length(unique(cell_to_BCs_corrected$cell_barcode)))
  for (j in 1:length(unique(cell_to_BCs_corrected$cell_barcode))){
    xx <- cell_to_BCs_corrected$cell_barcode ==
      unique(cell_to_BCs_corrected$cell_barcode)[j]

    cell_to_BCs_corrected$proportion_reads_in_cell[xx] <-
      cell_to_BCs_corrected$n_corrected[xx]/sum(cell_to_BCs_corrected$n_corrected[xx])

    cli_progress_update() # Update progress
  }
  cli_progress_done() # Close the progress bar

  # Filters barcode assignments by read count and read count proportion thresholds
  cell_to_BCs_filtered <- cell_to_BCs_corrected %>%
    dplyr::filter(n_corrected >= reads_thresholds$min_reads &
                    proportion_reads_in_cell >= reads_thresholds$min_reads_proportion)

  # Creates a combined barcode string for each assignment
  cell_to_BCs_filtered <- cell_to_BCs_filtered %>%
    dplyr::mutate(all_BCs = paste(BC1, BC2, BC3, sep = "_"))

  message("Cell Barcodes Assigned.")

  head(cell_to_BCs_filtered)

  # -- Optional file saving (only if output_dir is specified) --
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    file_save_name <- file.path(output_dir, "cell_to_UTR_barcode_system.csv")

    write.table(cell_to_BCs_filtered, file_save_name, sep = ",",
                col.names = TRUE, row.names = FALSE, quote = FALSE)

    message(glue::glue("Cell to UTR barcode system saved to: ",
                       "{normalizePath(output_dir)}."))
  }

  return(cell_to_BCs_filtered)
}


#' Barcode Correction
#'
#' This function corrects barcode assignments by applying Hamming distance-based
#' threshold correction.
#'
#' @param cell_to_extracted_BCs Data frame linking cells to extracted barcodes.
#' @param list_BC1 Named list of nearest BC1 barcode candidates.
#' @param list_BC2_3 Named list of nearest BC2/BC3 barcode candidates.
#' @param BC_thresholds List containing the following barcode thresholds:
#' `BC1_threshold`, `BC2_3_threshold`.
#'
#' @return A data frame of corrected barcode IDs, grouped by cell barcode and
#' barcode combination, and the corrected number of counts.

#' @examples
#' correct_barcodes(
#' cell_to_extracted_BCs = cell_to_extracted_BCs,
#' list_BC1 = list_BC1,
#' list_BC2_3 = list_BC2_3,
#' BC_thresholds = list(BC1_threshold = 3, BC2_3_threshold = 5)
#' )

#' @export
correct_barcodes <- function(
  cell_to_extracted_BCs,
  list_BC1,
  list_BC2_3,
  BC_thresholds
) {
  # Extract correction thresholds
  BC1_threshold <- BC_thresholds$BC1_threshold
  BC2_3_threshold <- BC_thresholds$BC2_3_threshold

  # Correct BC1 barcodes using threshold BC1_threshold
  df_correction_BC1 <- data.frame(
    BC1_extracted = names(list_BC1),
    BC1_corrected_ID = unlist(lapply(list_BC1, correct_BC,
                                     thresh = BC1_threshold))
  )

  # Correct BC2_3 barcodes using threshold BC2_3_threshold
  df_correction_BC2_3 <- data.frame(
    BC2_3_extracted = names(list_BC2_3),
    BC2_3_corrected_ID = unlist(lapply(list_BC2_3, correct_BC,
                                       thresh = BC2_3_threshold))
  )

  # Map corrected barcodes back to each cell by matching IDs
  cell_to_BCs <- data.frame(
    cell_barcode = cell_to_extracted_BCs$cell_barcode,
    BC1 = df_correction_BC1$BC1_corrected_ID[match(
      cell_to_extracted_BCs$BC1,df_correction_BC1$BC1_extracted
    )],
    BC2 = df_correction_BC2_3$BC2_3_corrected_ID[match(
      cell_to_extracted_BCs$BC2,df_correction_BC2_3$BC2_3_extracted
    )],
    BC3 = df_correction_BC2_3$BC2_3_corrected_ID[match(
      cell_to_extracted_BCs$BC3,df_correction_BC2_3$BC2_3_extracted
    )],
    n = cell_to_extracted_BCs$n
  )

  # Filter out rows with missing barcode corrections
  cell_to_BCs_all_BCs <- cell_to_BCs %>%
    dplyr::filter(BC1 != "", BC2 != "", BC3 != "")

  # Group by full barcode combination and sum the counts
  cell_to_BCs_corrected <- cell_to_BCs_all_BCs %>%
    dplyr::group_by(cell_barcode, BC1, BC2, BC3) %>%
    dplyr::summarise(n_corrected = sum(n), .groups = "drop")

  return(cell_to_BCs_corrected)
}
