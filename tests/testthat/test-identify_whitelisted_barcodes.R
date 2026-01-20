test_that("import_whitelisted_barcodes returns a list with correct element names", {

  # Mock return for BC1 whitelist
  mockery::stub(
    import_whitelisted_barcodes,
    "import_BC1_whitelist",
    c("A", "B")
  )

  # Mock return for BC2/BC3 whitelist (data frame)
  mockery::stub(
    import_whitelisted_barcodes,
    "import_BC2_BC3_whitelist",
    data.frame(V1 = c("X", "Y"), stringsAsFactors = FALSE)
  )

  result <- import_whitelisted_barcodes(
    source_dir = "dummy",
    BC1_whitelist_file_name = "file1",
    BC2_BC3_whitelist_file_name = "file2",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("BC1_whitelist", "BC2_BC3_whitelist"))
})

test_that("BC1_whitelist gets names BC1_1, BC1_2, etc of correct length", {

  # Mock BC1 whitelist with two items
  mockery::stub(
    import_whitelisted_barcodes,
    "import_BC1_whitelist",
    c("AAA", "BBB")
  )

  # Mock BC2/BC3 whitelist (minimal data frame)
  mockery::stub(
    import_whitelisted_barcodes,
    "import_BC2_BC3_whitelist",
    data.frame(V1 = "XXX", stringsAsFactors = FALSE)
  )

  result <- import_whitelisted_barcodes(
    source_dir = "dummy",
    BC1_whitelist_file_name = "file1",
    BC2_BC3_whitelist_file_name = "file2",
    verbose = FALSE
  )

  expect_equal(
    names(result$BC1_whitelist),
    c("BC1_1", "BC1_2")
  )

  # Values are unchanged
  expect_equal(
    unname(result$BC1_whitelist),
    c("AAA", "BBB")
  )
})

test_that("identify_BC1_closest_barcodes returns a named list keyed by unique BC1s", {

  # Small input with repeated BC1 so we can check `unique()`
  cell_to_extracted_BCs <- data.frame(
    BC1 = c("AAAA", "BBBB", "AAAA"),
    stringsAsFactors = FALSE
  )

  # Whitelist: same length strings (needed for Hamming distance)
  BC1_whitelist <- c("CCCC", "DDDD", "EEEE")
  names(BC1_whitelist) <- c("BC1_1", "BC1_2", "BC1_3")

  # Mock extract_closest_barcodes so we don't test its internals here
  mockery::stub(
    identify_BC1_closest_barcodes,
    "extract_closest_barcodes",
    function(i, dist_matrix, whitelist_names, n_barcodes) {
      # return the first min(n_barcodes, length(whitelist_names)) names
      whitelist_names[seq_len(min(n_barcodes, length(whitelist_names)))]
    }
  )

  result <- identify_BC1_closest_barcodes(
    cell_to_extracted_BCs = cell_to_extracted_BCs,
    BC1_whitelist = BC1_whitelist,
    verbose = FALSE
  )

  unique_BC1 <- unique(cell_to_extracted_BCs$BC1)

  # Type: must be a list
  expect_type(result, "list")

  # Length: one entry per unique BC1
  expect_length(result, length(unique_BC1))

  # Names: exactly the unique BC1 values, in order
  expect_equal(names(result), unique_BC1)
})

test_that("identify_BC1_closest_barcodes uses extract_closest_barcodes output for each list element", {

  # Two unique BC1s so we expect two list elements
  cell_to_extracted_BCs <- data.frame(
    BC1 = c("AAAA", "BBBB"),
    stringsAsFactors = FALSE
  )

  # Whitelist: same length sequences (needed for Hamming)
  BC1_whitelist <- c("CCCC", "DDDD")
  names(BC1_whitelist) <- c("BC1_1", "BC1_2")

  # Stub extract_closest_barcodes so we control its output
  mockery::stub(
    identify_BC1_closest_barcodes,
    "extract_closest_barcodes",
    function(i, dist_matrix, whitelist_names, n_barcodes) {
      # Return something that depends only on i
      paste0("mock_result_", i)
    }
  )

  result <- identify_BC1_closest_barcodes(
    cell_to_extracted_BCs = cell_to_extracted_BCs,
    BC1_whitelist = BC1_whitelist,
    verbose = FALSE
  )

  # We expect one list element per row of the distance matrix,
  # i.e. per unique BC1, with the value from our mock
  expect_identical(
    unname(result),
    list("mock_result_1", "mock_result_2")
  )
})

test_that("duplicate BC1 values result in one list element per unique BC1", {

  # Note: "AAAA" appears twice, "BBBB" once
  cell_to_extracted_BCs <- data.frame(
    BC1 = c("AAAA", "BBBB", "AAAA"),
    stringsAsFactors = FALSE
  )

  # Whitelist: sequences of equal length for Hamming
  BC1_whitelist <- c("CCCC", "DDDD")
  names(BC1_whitelist) <- c("BC1_1", "BC1_2")

  # Stub extract_closest_barcodes: we don't care what it returns here,
  # only that it is called the right number of times
  mockery::stub(
    identify_BC1_closest_barcodes,
    "extract_closest_barcodes",
    function(i, dist_matrix, whitelist_names, n_barcodes) {
      paste0("mock_", i)
    }
  )

  result <- identify_BC1_closest_barcodes(
    cell_to_extracted_BCs = cell_to_extracted_BCs,
    BC1_whitelist = BC1_whitelist,
    verbose = FALSE
  )

  unique_BC1 <- unique(cell_to_extracted_BCs$BC1)

  # One list element per unique BC1
  expect_length(result, length(unique_BC1))

  # Names are exactly the unique BC1 values
  expect_equal(names(result), unique_BC1)
})

test_that("identify_BC2_3_closest_barcodes returns a named list keyed by unique BC2/3s", {

  # Small input with BC2 and BC3 values
  cell_to_extracted_BCs <- data.frame(
    BC2 = c("B2_A", "B2_B"),
    BC3 = c("B3_A", "B3_B"),
    stringsAsFactors = FALSE
  )

  # Whitelist: same-length sequences (needed for Hamming)
  BC2_BC3_whitelist <- c("W_A", "W_B")
  names(BC2_BC3_whitelist) <- c("BC2_3_1", "BC2_3_2")

  # Stub extract_closest_barcodes so we don't depend on its logic here
  mockery::stub(
    identify_BC2_3_closest_barcodes,
    "extract_closest_barcodes",
    function(i, dist_matrix, whitelist_names, n_barcodes) {
      paste0("mock_result_", i)
    }
  )

  result <- identify_BC2_3_closest_barcodes(
    cell_to_extracted_BCs = cell_to_extracted_BCs,
    BC2_BC3_whitelist = BC2_BC3_whitelist,
    verbose = FALSE
  )

  unique_BC_2_3 <- unique(c(cell_to_extracted_BCs$BC2, cell_to_extracted_BCs$BC3))

  # Type: must be a list
  expect_type(result, "list")

  # Length: one entry per unique BC2/3
  expect_length(result, length(unique_BC_2_3))

  # Names: exactly the unique BC2/3 values, in order
  expect_equal(names(result), unique_BC_2_3)
})

test_that("identify_BC2_3_closest_barcodes uses extract_closest_barcodes output for each list element", {

  # Two BC2s and two BC3s → 4 unique BC2/3s in total
  cell_to_extracted_BCs <- data.frame(
    BC2 = c("B2_A", "B2_B"),
    BC3 = c("B3_A", "B3_B"),
    stringsAsFactors = FALSE
  )

  # Whitelist: same-length sequences (needed for Hamming)
  BC2_BC3_whitelist <- c("W_A", "W_B")
  names(BC2_BC3_whitelist) <- c("BC2_3_1", "BC2_3_2")

  # Stub extract_closest_barcodes so we fully control its output
  mockery::stub(
    identify_BC2_3_closest_barcodes,
    "extract_closest_barcodes",
    function(i, dist_matrix, whitelist_names, n_barcodes) {
      paste0("mock_result_", i)
    }
  )

  result <- identify_BC2_3_closest_barcodes(
    cell_to_extracted_BCs = cell_to_extracted_BCs,
    BC2_BC3_whitelist = BC2_BC3_whitelist,
    verbose = FALSE
  )

  # unique_BC_2_3 has length 4 here: B2_A, B2_B, B3_A, B3_B
  unique_BC_2_3 <- unique(c(cell_to_extracted_BCs$BC2, cell_to_extracted_BCs$BC3))
  expect_length(result, length(unique_BC_2_3))

  # The un-named list elements should be exactly what our stub returns
  expect_identical(
    unname(result),
    as.list(paste0("mock_result_", seq_along(unique_BC_2_3)))
  )
})

test_that("identify_BC2_3_closest_barcodes includes BC2 and BC3 and de-duplicates across them", {

  # Here BC2 and BC3 share a value "SHARED"
  cell_to_extracted_BCs <- data.frame(
    BC2 = c("B2_ONLY", "SHARED"),
    BC3 = c("SHARED", "B3_ONLY"),
    stringsAsFactors = FALSE
  )

  # unique(c(BC2, BC3)) should be: "B2_ONLY", "SHARED", "B3_ONLY"
  expected_unique <- unique(c(cell_to_extracted_BCs$BC2, cell_to_extracted_BCs$BC3))

  # Whitelist: same-length sequences (needed for Hamming)
  BC2_BC3_whitelist <- c("W_A", "W_B")
  names(BC2_BC3_whitelist) <- c("BC2_3_1", "BC2_3_2")

  # Stub extract_closest_barcodes; here assume the internal is correct
  mockery::stub(
    identify_BC2_3_closest_barcodes,
    "extract_closest_barcodes",
    function(i, dist_matrix, whitelist_names, n_barcodes) {
      paste0("mock_", i)
    }
  )

  result <- identify_BC2_3_closest_barcodes(
    cell_to_extracted_BCs = cell_to_extracted_BCs,
    BC2_BC3_whitelist = BC2_BC3_whitelist,
    verbose = FALSE
  )

  # One entry per unique BC2/3 value across both columns
  expect_length(result, length(expected_unique))

  # Names must match unique(c(BC2, BC3)) exactly
  expect_equal(names(result), expected_unique)
})


