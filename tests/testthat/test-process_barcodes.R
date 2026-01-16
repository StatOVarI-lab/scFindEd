test_that("import_process_barcodes_data returns a list with correct names", {

  # Mock the import functions so the test is isolated
  mockery::stub(
    import_process_barcodes_data,
    "import_barcode_amplicon_info",
    data.frame(readname = character())
  )

  mockery::stub(
    import_process_barcodes_data,
    "import_barcode_reads_data",
    character()
  )

  mockery::stub(
    import_process_barcodes_data,
    "import_readnames_data",
    character()
  )

  result <- import_process_barcodes_data(
    source_dir = "dummy",
    barcode_amplicon_file_name = "a.txt",
    cell_barcodes_file_name = "b.txt",
    readnames_file_name = "c.txt",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("barcode_amplicon_data", "cell_barcodes", "readnames")
  )
})

test_that("cell_barcodes is a named vector with names given by readnames", {

  mockery::stub(
    import_process_barcodes_data,
    "import_barcode_amplicon_info",
    data.frame(readname = c("r1", "r2"))
  )

  mockery::stub(
    import_process_barcodes_data,
    "import_barcode_reads_data",
    c("BC01", "BC02")
  )

  mockery::stub(
    import_process_barcodes_data,
    "import_readnames_data",
    c("r1", "r2")
  )

  result <- import_process_barcodes_data(
    source_dir = "dummy",
    barcode_amplicon_file_name = "a",
    cell_barcodes_file_name = "b",
    readnames_file_name = "c",
    verbose = FALSE
  )

  expect_type(result$cell_barcodes, "character")
  expect_named(result$cell_barcodes, c("r1", "r2"))
})

test_that("barcode_amplicon_data has readname and cell_barcode", {

  mockery::stub(
    import_process_barcodes_data,
    "import_barcode_amplicon_info",
    data.frame(readname = c("r1", "r2"))
  )

  mockery::stub(
    import_process_barcodes_data,
    "import_barcode_reads_data",
    c("BC01", "BC02")
  )

  mockery::stub(
    import_process_barcodes_data,
    "import_readnames_data",
    c("r1", "r2")
  )

  result <- import_process_barcodes_data(
    source_dir = "dummy",
    barcode_amplicon_file_name = "a",
    cell_barcodes_file_name = "b",
    readnames_file_name = "c",
    verbose = FALSE
  )

  expect_type(result$barcode_amplicon_data, "list")
  expect_named(
    result$barcode_amplicon_data,
    c("readname", "cell_barcode")
  )
})


test_that("filter_amplicon_data returns a list with the correct names", {

  # Minimal input data frame with required columns
  barcode_amplicon_data <- data.frame(
    cell_barcode = c("BC01", "BC01"),
    quality      = c(30, 30),
    start        = c(40, 145),
    readname     = c("r1", "r2"),
    stringsAsFactors = FALSE
  )

  result <- filter_amplicon_data(
    barcode_amplicon_data = barcode_amplicon_data,
    min_reads = 0,
    quality_score_threshold = 0,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("barcode_amplicon_data_forward", "barcode_amplicon_data_reverse")
  )
})

test_that("filter_amplicon_data returns data frames for forward and reverse", {

  barcode_amplicon_data <- data.frame(
    cell_barcode = c("BC01", "BC01"),
    quality      = c(30, 30),
    start        = c(40, 145),
    readname     = c("r1", "r2"),
    stringsAsFactors = FALSE
  )

  result <- filter_amplicon_data(
    barcode_amplicon_data = barcode_amplicon_data,
    min_reads = 0,
    quality_score_threshold = 0,
    verbose = FALSE
  )

  expect_s3_class(
    result$barcode_amplicon_data_forward,
    "data.frame"
  )

  expect_s3_class(
    result$barcode_amplicon_data_reverse,
    "data.frame"
  )
})


test_that("extract_barcodes returns a data frame", {

  # Minimal forward data
  barcode_amplicon_data_forward <- data.frame(
    readname = "r1",
    cell_barcode = "BC01",
    sequence = "AAA",
    stringsAsFactors = FALSE
  )

  # Minimal reverse data
  barcode_amplicon_data_reverse <- data.frame(
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
    barcode_amplicon_data_forward = barcode_amplicon_data_forward,
    barcode_amplicon_data_reverse = barcode_amplicon_data_reverse,
    seq_before_barcode_1 = "A",
    seq_before_barcode_2 = "A",
    seq_after_barcode_3 = "T",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
})

test_that("extract_barcodes returns required columns", {

  barcode_amplicon_data_forward <- data.frame(
    readname = "r1",
    cell_barcode = "BC01",
    sequence = "AAA",
    stringsAsFactors = FALSE
  )

  barcode_amplicon_data_reverse <- data.frame(
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
    barcode_amplicon_data_forward,
    barcode_amplicon_data_reverse,
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




