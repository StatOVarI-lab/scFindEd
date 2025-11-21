# Function to compute Hamming distances between a barcode and a whitelist,
# returning the 3 smallest distances
identify_barcode_ID <- function(BC, BC_vec){
  temp <- sapply(BC_vec,
                 function(x) stringdist::stringdist(x, BC, method = "hamming"))
  return(sort(temp)[1:3])
}
