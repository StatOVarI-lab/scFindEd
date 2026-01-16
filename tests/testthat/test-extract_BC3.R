test_that("extract_BC3 extracts the correct 30 bp sequence when flank is present", {
  seq_after_barcode_3 <- "TTGAC"
  flank_30bp <- paste(rep("A", 30), collapse = "")
  prefix <- paste0("GGGCCC", flank_30bp)
  full_seq <- paste0(prefix, seq_after_barcode_3, "AGTCAGTC")

  result <- extract_BC3(full_seq, seq_after_barcode_3)

  expect_equal(result, flank_30bp)
})

test_that("extract_BC3 throws an error if the flanking sequence is missing", {

  seq_after <- "TTT"
  full_seq <- "AAAA"

  expect_error(
    extract_BC3(full_seq, seq_after),
    "flanking sequence"
  )
})
