#' @keywords internal
# Checks if both seq1 and seq2 are present in each element of vector v
grepl_seq1_2 <- function(v,seq1,seq2){
  temp_1 <- sapply(v,function(x) grepl(seq1,x,fixed = TRUE)) # Check for presence of seq1
  temp_2 <- sapply(v,function(x) grepl(seq2,x,fixed = TRUE)) # Check for presence of seq2
  return(temp_1 & temp_2) # Return TRUE when both temp_1 and temp_2 are TRUE
}

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
