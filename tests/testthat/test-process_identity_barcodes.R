test_that("filter_identity_amplicon_reads returns a list with the correct names", {

  # Minimal input data frame with required columns
  identity_amplicon_reads <- data.frame(
    cell_barcode = c("BC01", "BC01"),
    mapq         = c(30, 30),
    pos          = c(40, 145),
    read_name    = c("r1", "r2"),
    stringsAsFactors = FALSE
  )

  result <- filter_identity_amplicon_reads(
    identity_amplicon_reads = identity_amplicon_reads,
    min_reads = 0,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("filtered_identity_amplicon_reads_forward",
      "filtered_identity_amplicon_reads_reverse")
  )
})

test_that("filter_identity_amplicon_reads returns data frames for forward and reverse", {

  identity_amplicon_reads <- data.frame(
    cell_barcode = c("BC01", "BC01"),
    mapq         = c(30, 30),
    pos          = c(40, 145),
    read_name    = c("r1", "r2"),
    stringsAsFactors = FALSE
  )

  result <- filter_identity_amplicon_reads(
    identity_amplicon_reads = identity_amplicon_reads,
    min_reads = 0,
    verbose = FALSE
  )

  expect_s3_class(
    result$filtered_identity_amplicon_reads_forward,
    "data.frame"
  )

  expect_s3_class(
    result$filtered_identity_amplicon_reads_reverse,
    "data.frame"
  )
})


# test_that("extract_UTR_barcodes returns a data frame", {

  # # Minimal forward data
  # UTR_barcode_reads_forward <- data.frame(
    # read_name    = "r1",
    # cell_barcode = "BC01",
    # seq          = "AAA",
    # stringsAsFactors = FALSE
  # )

  # # Minimal reverse data
  # UTR_barcode_reads_reverse <- data.frame(
    # read_name = "r1",
    # seq       = "TTT",
    # stringsAsFactors = FALSE
  # )

  # result <- extract_UTR_barcodes(
    # UTR_barcode_reads_forward = UTR_barcode_reads_forward,
    # UTR_barcode_reads_reverse = UTR_barcode_reads_reverse,
    # UTR_barcode_info = ,
    # verbose = FALSE
  # )

  # expect_s3_class(result, "data.frame")
# })

# test_that("extract_UTR_barcodes returns required columns", {

  # UTR_barcode_reads_forward <- data.frame(
    # read_name    = "r1",
    # cell_barcode = "BC01",
    # seq          = "AAA",
    # stringsAsFactors = FALSE
  # )

  # UTR_barcode_reads_reverse <- data.frame(
    # read_name = "r1",
    # seq       = "TTT",
    # stringsAsFactors = FALSE
  # )

  # result <- extract_UTR_barcodes(
    # UTR_barcode_reads_forward = UTR_barcode_reads_forward,
    # UTR_barcode_reads_reverse = UTR_barcode_reads_reverse,
    # UTR_barcode_info = ,
    # verbose = FALSE
  # )

  # expect_true(
    # all(c("read_name", "cell_barcode", "BC1", "BC2", "BC3") %in% colnames(result))
  # )
# })

test_that("group_identity_barcodes returns a data frame with required columns", {

  readname_barcodes_table <- data.frame(
    cell_barcode = c("CELL1", "CELL1"),
    BC1 = c(strrep("A", 20), strrep("A", 20)),
    BC2 = c(strrep("C", 30), strrep("C", 30)),
    BC3 = c(strrep("G", 30), strrep("G", 30)),
    stringsAsFactors = FALSE
  )

  result <- group_identity_barcodes(readname_barcodes_table, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("cell_barcode", "BC1", "BC2", "BC3", "n") %in% names(result)))
})




