#'
#'
#' @export
process_reads <- function(
    source_dir = "",
    sample_name = "sample",
    amplicon_gRNA_number = "n",
    amplicon_data_file_path,
    output_dir = NULL
) {
  message(glue::glue(
    "-- Processing amplicon_{amplicon_gRNA_number} from {sample_name} file --"
  ))

  # Read in the amplicon data
  amplicon_data <- read.table(amplicon_data_file_path, header = TRUE)
  amplicon_data <- amplicon_data %>%
    mutate(amplicon_gRNA = paste0(AMPLICON_NAME, "_", sgRNA_SEQUENCE)) %>%
    filter(!is.na(sgRNA_SEQUENCE))
}
