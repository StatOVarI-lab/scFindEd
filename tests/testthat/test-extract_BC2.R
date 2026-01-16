test_that("extract_BC2 extracts the correct 30bp barcode when flank is present", {
  seq_before <- "GGG"
  barcode <- paste0(rep("A", 30), collapse = "")  # 30 As

  full_seq <- paste0(seq_before, barcode, "TTTT")

  result <- extract_BC2(full_seq, seq_before)

  expect_identical(result, barcode)
})

test_that("extract_BC2 throws an error if the flanking sequence is missing", {

  seq_before <- "GGG"
  full_seq <- "TTTAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # no 'GGG' here

  expect_error(
    extract_BC2(full_seq, seq_before),
    "flanking sequence"
  )
})
