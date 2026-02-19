test_that("filter_cell_barcode_reads returns a list with the correct names", {

  # Minimal input data frame with required columns
  mapped_cell_barcode_reads <- data.frame(
    cell_barcode = c("BC01", "BC01"),
    quality      = c(30, 30),
    start        = c(40, 145),
    readname     = c("r1", "r2"),
    stringsAsFactors = FALSE
  )

  result <- filter_cell_barcode_reads(
    mapped_cell_barcode_reads = mapped_cell_barcode_reads,
    min_reads = 0,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("filtered_cell_barcode_reads_forward", "filtered_cell_barcode_reads_reverse")
  )
})

test_that("filter_cell_barcode_reads returns data frames for forward and reverse", {

  mapped_cell_barcode_reads <- data.frame(
    cell_barcode = c("BC01", "BC01"),
    quality      = c(30, 30),
    start        = c(40, 145),
    readname     = c("r1", "r2"),
    stringsAsFactors = FALSE
  )

  result <- filter_cell_barcode_reads(
    mapped_cell_barcode_reads = mapped_cell_barcode_reads,
    min_reads = 0,
    verbose = FALSE
  )

  expect_s3_class(
    result$filtered_cell_barcode_reads_forward,
    "data.frame"
  )

  expect_s3_class(
    result$filtered_cell_barcode_reads_reverse,
    "data.frame"
  )
})


test_that("extract_barcodes returns a data frame", {

  # Minimal forward data
  cell_barcode_reads_forward <- data.frame(
    readname = "r1",
    cell_barcode = "BC01",
    sequence = "AAA",
    stringsAsFactors = FALSE
  )

  # Minimal reverse data
  cell_barcode_reads_reverse <- data.frame(
    readname = "r1",
    sequence = "TTT",
    stringsAsFactors = FALSE
  )

  # Mock helper functions so logic is controlled
  mockery::stub(
    extract_barcodes,
    "grepl_seq1_2",
    function(x, ...) TRUE
  )

  mockery::stub(
    extract_barcodes,
    "extract_BC1",
    function(x, ...) "BC1"
  )

  mockery::stub(
    extract_barcodes,
    "extract_BC2",
    function(x, ...) "BC2"
  )

  mockery::stub(
    extract_barcodes,
    "extract_BC3",
    function(x, ...) "BC3"
  )

  result <- extract_barcodes(
    cell_barcode_reads_forward = cell_barcode_reads_forward,
    cell_barcode_reads_reverse = cell_barcode_reads_reverse,
    seq_before_barcode_1 = "A",
    seq_before_barcode_2 = "A",
    seq_after_barcode_3 = "T",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
})

test_that("extract_barcodes returns required columns", {

  cell_barcode_reads_forward <- data.frame(
    readname = "r1",
    cell_barcode = "BC01",
    sequence = "AAA",
    stringsAsFactors = FALSE
  )

  cell_barcode_reads_reverse <- data.frame(
    readname = "r1",
    sequence = "TTT",
    stringsAsFactors = FALSE
  )

  mockery::stub(
    extract_barcodes,
    "grepl_seq1_2",
    function(x, ...) TRUE
  )

  mockery::stub(
    extract_barcodes,
    "extract_BC1",
    function(x, ...) "BC1"
  )

  mockery::stub(
    extract_barcodes,
    "extract_BC2",
    function(x, ...) "BC2"
  )

  mockery::stub(
    extract_barcodes,
    "extract_BC3",
    function(x, ...) "BC3"
  )

  result <- extract_barcodes(
    cell_barcode_reads_forward = cell_barcode_reads_forward,
    cell_barcode_reads_reverse = cell_barcode_reads_reverse,
    seq_before_barcode_1 = "A",
    seq_before_barcode_2 = "A",
    seq_after_barcode_3 = "T",
    verbose = FALSE
  )

  expect_true(
    all(c("readname", "cell_barcode", "BC1", "BC2", "BC3") %in% colnames(result))
  )
})

test_that("group_barcodes returns a data frame with required columns", {

  readname_cell_table <- data.frame(
    cell_barcode = c("CELL1", "CELL1"),
    BC1 = c(strrep("A", 20), strrep("A", 20)),
    BC2 = c(strrep("C", 30), strrep("C", 30)),
    BC3 = c(strrep("G", 30), strrep("G", 30)),
    stringsAsFactors = FALSE
  )

  result <- group_barcodes(readname_cell_table, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("cell_barcode", "BC1", "BC2", "BC3", "n") %in% names(result)))
})




