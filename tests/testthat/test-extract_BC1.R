test_that("extract_BC1 extracts the correct 20bp barcode when flank is present", {

  seq_before <- "AAA"
  barcode <- paste0(rep("C", 20), collapse = "")
  full_seq <- paste0(seq_before, barcode, "TTTT")  # trailing bases allowed

  result <- extract_BC1(full_seq, seq_before)

  expect_identical(result, barcode)
})

test_that("extract_BC1 throws an error if the flanking sequence is missing", {

  seq_before <- "AAA"
  full_seq <- "GGGCCCCCCCCCCCCCCCCCCCCGGG"

  expect_error(
    extract_BC1(full_seq, seq_before),
    "flanking sequence"
  )
})
