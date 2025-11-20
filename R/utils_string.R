#' @keywords internal
# Checks if both seq1 and seq2 are present in each element of vector v
grepl_seq1_2 <- function(v,seq1,seq2){
  temp_1 <- sapply(v,function(x) grepl(seq1,x,fixed = TRUE)) # Check for presence of seq1
  temp_2 <- sapply(v,function(x) grepl(seq2,x,fixed = TRUE)) # Check for presence of seq2
  return(temp_1 & temp_2) # Return TRUE when both temp_1 and temp_2 are TRUE
}

# Extracts 20 bases following seq_before_barcode_1
extract_BC1 <- function(x,seq_before_barcode_1){
  parts <- strsplit(x, seq_before_barcode_1)[[1]]

  # If no split occurred, flank missing
  if (length(parts) < 2) {
    stop("extract_BC1(): flanking sequence not found.
         This likely indicates an upstream filtering issue.")
  }

  temp <- parts[2]
  substr(temp, 1, 20)
}

# Extracts 30 bases following seq_before_barcode_2
extract_BC2 <- function(x,seq_before_barcode_2){
  temp <- strsplit(x,seq_before_barcode_2)[[1]][2]
  return(substr(temp,1,30))
}

# Extracts 30 bases preceding seq_after_barcode_3
extract_BC3 <- function(x,seq_after_barcode_3){
  temp <- strsplit(x,seq_after_barcode_3)[[1]][1]
  nc <- nchar(temp)
  return(substr(temp,nc-29,nc))
}

# Function to compute Hamming distances between a barcode and a whitelist,
# returning the 3 smallest distances
identify_barcode_ID <- function(BC,BC_vec){
  temp <- sapply(BC_vec,function(x) stringdist(x,BC,method="hamming"))
  return(sort(temp)[1:3])
}

