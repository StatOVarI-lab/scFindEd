#' @keywords internal
#' Imports data from a file
#' Expected file name: "reads_classified_{sample_name}HDR_{amplicon_gRNA_number}.csv"
import_amplicon_summary_file <- function(
  source_dir,
  file_name,
  verbose
) {
  file_path <- file.path(source_dir, file_name)
  if (verbose) cli::cli_inform("Importing {file_path}...")

  tryCatch(
    data <- as.data.frame(readr::read_csv(file_path, show_col_types = FALSE)),
    error = function(e) {
      cli::cli_abort(
        "Failed to import {file_name} from {source_dir}.")
    }
  )

  return(data)
}
