test_that("grepl_seq1_2 correctly detects the presence of both motifs", {
  seq1 <- "AAA"
  seq2 <- "TTT"

  input <- c(
    "AAATTT",  # contains both
    "AAAGGG",  # seq1 only
    "GGGTTT",  # seq2 only
    "GGGCCC"   # neither
  )

  result <- grepl_seq1_2(input, seq1, seq2)

  expect_identical(
    unname(result),
    c(TRUE, FALSE, FALSE, FALSE)
  )
})
