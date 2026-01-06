#' Example barcode sequence
#'
#' A 10-base DNA barcode used for demonstrating flanking sequence extraction,
#' validating the presence of sequences within a vector, and calculating and
#' comparing Hamming distances.
#'
#' @format A character string of length 10.
#' @examples
#' grepl_seq1_2(x_example, barcode_example, barcode_example_2)
#' extract_BC1(x_example, barcode_example)
#' identify_barcode_ID(barcode_example, whitelist_example)
"barcode_example"


#' Secondary example barcode sequence
#'
#' A 10-base DNA barcode used for validating the presence of sequences within a
#' vector.
#'
#' @format A character string of length 10.
#' @examples
#' grepl_seq1_2(x_example, barcode_example, barcode_example_2)
"barcode_example_2"


#' Example vector with sequences containing both barcodes
#'
#' A character vector of DNA sequences, each 30 bases long, that contain both
#' `ATACGCATCA` and `GGTACCTAGC`.
#'
#' @format A character vector of length 5.
#' @examples
#' grepl_seq1_2(v_TRUE_example, "ATACGCATCA", "GGTACCTAGC")
"v_TRUE_example"


#' Example vector with sequences missing one or both barcodes
#'
#' A character vector of DNA sequences, each 30 bases long, that do not all
#' contain both `ATACGCATCA` and `GGTACCTAGC`.
#'
#' @format A character vector of length 5.
#' @examples
#' grepl_seq1_2(v_FALSE_example, "ATACGCATCA", "GGTACCTAGC")
"v_FALSE_example"


#' Example DNA sequence for searching for barcode flanking sequences
#'
#' A single DNA string of 90 bases used to demonstrate the extraction of a
#' barcode's flanking sequences.
#'
#' @format A character string of length 90.
#' @examples
#' extract_BC1(x_example, barcode_example)
"x_example"


#' Example barcode whitelist
#'
#' A character vector of 10 DNA barcodes used for demonstrating Hamming distance
#' calculation and comparison.
#'
#' @format A character vector of length 10.
#' @examples
#' identify_barcode_ID(barcode_example, whitelist_example)
"whitelist_example"





