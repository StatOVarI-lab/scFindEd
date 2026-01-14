#' @keywords internal
#' Imports data from a file
#' Expected file name: "barcode_amplicon_info.txt"
import_barcode_amplicon_info <- function(
    source_dir,
    file_name,
    verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")
  
  # Import barcode amplicon data table
  tryCatch(
    data <- readr::read_table(file_path, col_names = FALSE,
                              show_col_types = FALSE),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )
  
  # Assign column names
  names(data) <- c("readname", "amplicon_name", "start",
                   "quality",  "cigar",         "sequence")
  
  if (verbose) cli::cli_alert_success(
    "Successfully imported {file_name} from {source_dir}.")
  
  return(data)
}


#' @keywords internal
#' Imports data from a file
#' Expected file name: "sample_name_barcodes_fastq.txt"
import_barcode_reads_data <- function(
    source_dir,
    file_name,
    verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")
  
  tryCatch(
    data <- readLines(file_path),
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
#' Expected file name: "sample_name_R2_001_readnames_fastq.txt"
import_readnames_data <- function(
    source_dir,
    file_name,
    verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")
  
  tryCatch(
    data <- data.table::fread(
      file_path, header = FALSE, showProgress = verbose)$V1,
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
# Checks if both seq1 and seq2 are present in each element of vector v
grepl_seq1_2 <- function(v,seq1,seq2){
  temp_1 <- sapply(v,function(x) grepl(seq1,x,fixed = TRUE)) # Check for presence of seq1
  temp_2 <- sapply(v,function(x) grepl(seq2,x,fixed = TRUE)) # Check for presence of seq2
  return(temp_1 & temp_2) # Return TRUE when both temp_1 and temp_2 are TRUE
}

#' @keywords internal
# Extracts up to 20 bases following seq_before_barcode_1
extract_BC1 <- function(x, seq_before_barcode_1){
  # Find seq_before_barcode_1 within x
  barcode_check <- regexpr(seq_before_barcode_1, x, fixed = TRUE)
  
  # Report error if seq_before_barcode_1 is not present inside of x
  if (barcode_check == -1) {
    cli::cli_abort("extract_BC1(): flanking sequence not found.
                   This likely indicates an upstream filtering issue.")
  }
  
  barcode_start <- barcode_check + attr(barcode_check, "match.length")
  
  # Extracts up to 20 bases following seq_before_barcode_1
  substr(x, barcode_start, barcode_start + 19)
}

#' @keywords internal
# Extracts up to 30 bases following seq_before_barcode_2
extract_BC2 <- function(x, seq_before_barcode_2){
  # Find seq_before_barcode_2 within x
  barcode_check <- regexpr(seq_before_barcode_2, x, fixed = TRUE)
  
  # Report error if seq_before_barcode_2 is not present inside of x
  if (barcode_check == -1) {
    cli::cli_abort("extract_BC2(): flanking sequence not found.
                   This likely indicates an upstream filtering issue.")
  }
  
  barcode_start <- barcode_check + attr(barcode_check, "match.length")
  
  # Extracts up to 30 bases following seq_before_barcode_2
  substr(x, barcode_start, barcode_start + 29)
}

#' @keywords internal
# Extracts up to 30 bases preceding seq_after_barcode_3
extract_BC3 <- function(x, seq_after_barcode_3){
  # Find seq_after_barcode_3 within x
  barcode_check <- regexpr(seq_after_barcode_3, x, fixed = TRUE)
  
  # Report error if seq_after_barcode_3 is not present inside of x
  if (barcode_check == -1) {
    cli::cli_abort("extract_BC3(): flanking sequence not found.
                   This likely indicates an upstream filtering issue.")
  }
  
  # Get position of the barcode start from the position of the barcode end
  barcode_end <- barcode_check - 1
  barcode_start <- max(1, barcode_end - 29)
  
  # Extracts up to 30 bases preceding seq_after_barcode_3
  substr(x, barcode_start, barcode_end)
}


#' @keywords internal
#' Imports data from files
#' Expected file name(s): "list_BC1.rds", "list_BC2_3.rds"
import_barcode_distance_lists <- function(
    source_dir,
    list_BC1_file_name,
    list_BC2_3_file_name,
    verbose
) {
  # -- Importing BC1 list --
  list_BC1 <- import_barcode_distance_list(
    source_dir,
    list_BC1_file_name,
    verbose
  )
  
  # -- Importing BC2/3 list --
  list_BC2_3 <- import_barcode_distance_list(
    source_dir,
    list_BC2_3_file_name,
    verbose
  )
  
  return(list(
    list_BC1 = list_BC1,
    list_BC2_3 = list_BC2_3
  ))
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
#' Converts a list of barcode distances into a data frame, better suited for
#' building the downstream QC ggplot histogram
prepare_hamming_plot_dataframe <- function(
    BC_list,
    BC_name,
    verbose
) {
  if (verbose) cli::cli_alert_info("Preparing data frame for {BC_name}...")
  
  # Distance from nearest whitelist barcode
  closest <- as.numeric(unlist(lapply(BC_list, function(x) x[1])))
  
  # Distance from second nearest whitelist barcode
  second_closest <- as.numeric(unlist(lapply(BC_list, function(x) x[2])))
  
  # Combine into a data frame
  BC_df <- data.frame(
    distance = c(closest, second_closest),
    type = factor(c(rep("Closest", length(closest)),
                    rep("Second Closest", length(second_closest))))
  )
  
  if (verbose) {
    cli::cli_alert_success("Successfully created data frame for {BC_name}.")
  }
  
  return(BC_df)
}


#' @keywords internal
#' Produces a histogram of Hamming distances to the nearest whitelist barcode.
barcode_hamming_plot <- function(
    data,
    title
) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = distance, fill = type)) +
    ggplot2::geom_histogram(position = "dodge", binwidth = 0.5, color = "black", alpha = 0.8) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = title,
      subtitle = "Distribution of distances by barcode match type",
      x = "Hamming Distance",
      y = "Number of Cells",
      fill = "Barcode Match"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor.y = ggplot2::element_line(color = "grey95", linewidth = 0.2)
    )
  
  return(plot)
}


