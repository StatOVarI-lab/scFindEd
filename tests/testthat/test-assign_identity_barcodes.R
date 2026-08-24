test_that("correct_identity_barcodes returns a data frame with expected columns", {

  # Minimal input: one cell with extracted BCs and a count
  cell_to_extracted_BCs <- data.frame(
    cell_barcode = "CELL1",
    BC1 = "BC1_raw",
    BC2 = "BC2_raw",
    BC3 = "BC3_raw",
    n   = 1L,
    stringsAsFactors = FALSE
  )

  # list_BC1, list_BC2 and list_BC3: names must match the extracted BCs
  list_BC1 <- list("dummy")        # value doesn't matter here
  names(list_BC1) <- "BC1_raw"

  list_BC2 <- list("dummy")        # value doesn't matter here
  names(list_BC2) <- "BC2_raw"

  list_BC3 <- list("dummy")        # value doesn't matter here
  names(list_BC3) <- "BC3_raw"

  # Collated barcode distance lists
  barcode_distance_lists <- list(
    BC1 = list_BC1,
    BC2 = list_BC2,
    BC3 = list_BC3
  )

  # Thresholds (values not important for this test)
  BC_thresholds <- list(
    BC1 = 1,
    BC2 = 2,
    BC3 = 3
  )

  # Stub correct_BC so that it always returns a non-empty ID
  mockery::stub(
    correct_identity_barcodes,
    "correct_BC",
    function(x, threshold) "CORRECTED"
  )

  result <- correct_identity_barcodes(
    cell_to_extracted_BCs  = cell_to_extracted_BCs,
    barcode_distance_lists = barcode_distance_lists,
    BC_thresholds          = BC_thresholds,
    verbose                = FALSE
  )

  # It should return a data frame / tibble
  expect_s3_class(result, "data.frame")

  # It should contain at least these columns
  expected_cols <- c("cell_barcode", "BC1", "BC2", "BC3", "n_corrected")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("correct_identity_barcodes maps corrected IDs back to each row via match()", {

  # Two cells with distinct extracted BCs
  cell_to_extracted_BCs <- data.frame(
    cell_barcode = c("CELL1", "CELL2"),
    BC1 = c("B1_raw1", "B1_raw2"),
    BC2 = c("B2_raw1", "B2_raw2"),
    BC3 = c("B3_raw1", "B3_raw2"),
    n   = c(5L, 7L),
    stringsAsFactors = FALSE
  )

  # list_BC1: names are extracted BC1s, values are arbitrary but distinct
  list_BC1 <- list("VAL_B1_raw1", "VAL_B1_raw2")
  names(list_BC1) <- c("B1_raw1", "B1_raw2")

  # list_BC2: names are extracted BC1s, values are arbitrary but distinct
  list_BC2 <- list("VAL_B2_raw1", "VAL_B2_raw2")
  names(list_BC2) <- c("B2_raw1", "B2_raw2")

  # list_BC3: names are extracted BC1s, values are arbitrary but distinct
  list_BC3 <- list("VAL_B3_raw1", "VAL_B3_raw2")
  names(list_BC3) <- c("B3_raw1", "B3_raw2")

  # Collated barcode distance lists
  barcode_distance_lists <- list(
    BC1 = list_BC1,
    BC2 = list_BC2,
    BC3 = list_BC3
  )

  BC_thresholds <- list(
    BC1 = 1,
    BC2 = 2,
    BC3 = 3
  )

  # Stub correct_BC so the "corrected ID" is a deterministic function of the value
  # correct_BC(x, threshold) -> paste0("CORR_", x)
  mockery::stub(
    correct_identity_barcodes,
    "correct_BC",
    function(x, threshold) paste0("CORR_", x)
  )

  result <- correct_identity_barcodes(
    cell_to_extracted_BCs  = cell_to_extracted_BCs,
    barcode_distance_lists = barcode_distance_lists,
    BC_thresholds          = BC_thresholds,
    verbose                = FALSE
  )

  # Build the expected mapping from extracted -> corrected IDs
  BC1_corrected_map <- setNames(
    paste0("CORR_", unlist(list_BC1)),
    names(list_BC1)      # B1_raw1, B1_raw2
  )

  BC2_corrected_map <- setNames(
    paste0("CORR_", unlist(list_BC2)),
    names(list_BC2)      # B2_raw1, B2_raw2
  )

  BC3_corrected_map <- setNames(
    paste0("CORR_", unlist(list_BC3)),
    names(list_BC3)      # B3_raw1, B3_raw2
  )

  # For each row, the final BCs should equal the corrected IDs
  # corresponding to the extracted BCs in that row
  expect_equal(
    as.character(result$BC1),
    as.character(BC1_corrected_map[cell_to_extracted_BCs$BC1])
  )

  expect_equal(
    as.character(result$BC2),
    as.character(BC2_corrected_map[cell_to_extracted_BCs$BC2])
  )

  expect_equal(
    result$BC3,
    as.character(BC3_corrected_map[cell_to_extracted_BCs$BC3])
  )
})

test_that("correct_identity_barcodes aggregates n_corrected as sum(n) per (cell_barcode, BC1, BC2, BC3)", {

  # Two rows that should end up in the SAME corrected barcode combination
  cell_to_extracted_BCs <- data.frame(
    cell_barcode = c("CELL1", "CELL1"),
    BC1 = c("B1_raw1", "B1_raw2"),
    BC2 = c("B2_raw1", "B2_raw2"),
    BC3 = c("B3_raw1", "B3_raw2"),
    n   = c(5L, 7L),
    stringsAsFactors = FALSE
  )

  # Names match extracted barcodes; values are what correct_BC() will transform into IDs
  list_BC1 <- list("ID_B1", "ID_B1")
  names(list_BC1) <- c("B1_raw1", "B1_raw2")

  list_BC2 <- list("ID_B2", "ID_B2")
  names(list_BC2) <- c("B2_raw1", "B2_raw2")

  list_BC3 <- list("ID_B3", "ID_B3")
  names(list_BC3) <- c("B3_raw1", "B3_raw2")

  # Collated barcode distance lists
  barcode_distance_lists <- list(
    BC1 = list_BC1,
    BC2 = list_BC2,
    BC3 = list_BC3
  )

  BC_thresholds <- list(
    BC1 = 1,
    BC2 = 2,
    BC3 = 3
  )

  # Stub correct_BC so corrected ID is based on the list value
  mockery::stub(
    correct_identity_barcodes,
    "correct_BC",
    function(x, threshold) as.character(x)
  )

  result <- correct_identity_barcodes(
    cell_to_extracted_BCs  = cell_to_extracted_BCs,
    barcode_distance_lists = barcode_distance_lists,
    BC_thresholds          = BC_thresholds,
    verbose                = FALSE
  )

  # Because both rows correct to the same (CELL1, ID_B1, ID_B23, ID_B23),
  # the result should have one row with n_corrected = 5 + 7 = 12
  expect_equal(nrow(result), 1L)
  expect_equal(result$cell_barcode, "CELL1")
  expect_equal(result$BC1, "ID_B1")
  expect_equal(result$BC2, "ID_B2")
  expect_equal(result$BC3, "ID_B3")
  expect_equal(result$n_corrected, 12L)
})

test_that("assign_identity_barcodes_to_cells computes proportion_reads_in_cell correctly per cell", {

  cell_to_BCs_corrected <- data.frame(
    cell_barcode = c("CELL1", "CELL1", "CELL2", "CELL2"),
    BC1 = c("A", "B", "C", "D"),
    BC2 = c("A", "B", "C", "D"),
    BC3 = c("A", "B", "C", "D"),
    n_corrected = c(1L, 3L, 2L, 2L),
    stringsAsFactors = FALSE
  )

  result <- assign_identity_barcodes_to_cells(cell_to_BCs_corrected, verbose = FALSE)

  # Column was added
  expect_true("proportion_reads_in_cell" %in% names(result))

  # CELL1: 1 / (1+3) = 0.25, 3 / (1+3) = 0.75
  cell1 <- result[result$cell_barcode == "CELL1", ]
  expect_equal(as.numeric(cell1$proportion_reads_in_cell), c(0.25, 0.75))

  # CELL2: 2 / (2+2) = 0.5, 2 / (2+2) = 0.5
  cell2 <- result[result$cell_barcode == "CELL2", ]
  expect_equal(as.numeric(cell2$proportion_reads_in_cell), c(0.5, 0.5))
})

test_that("assign_identity_barcodes_to_cells proportions sum to 1 within each cell", {

  cell_to_BCs_corrected <- data.frame(
    cell_barcode = c("CELL1", "CELL1", "CELL1", "CELL2", "CELL2"),
    BC1 = c("A", "B", "C", "D", "E"),
    BC2 = c("A", "B", "C", "D", "E"),
    BC3 = c("A", "B", "C", "D", "E"),
    n_corrected = c(2L, 3L, 5L, 1L, 4L),
    stringsAsFactors = FALSE
  )

  result <- assign_identity_barcodes_to_cells(cell_to_BCs_corrected, verbose = FALSE)

  sums <- tapply(
    as.numeric(result$proportion_reads_in_cell),
    result$cell_barcode,
    sum
  )

  # Floating point safety: allow a tiny tolerance
  expect_equal(as.numeric(sums[["CELL1"]]), 1, tolerance = 1e-12)
  expect_equal(as.numeric(sums[["CELL2"]]), 1, tolerance = 1e-12)
})

test_that("filter_assigned_barcodes filters by read count and proportion with inclusive thresholds", {

  cell_to_assigned_BCs <- data.frame(
    cell_barcode = c("CELL1", "CELL1", "CELL1", "CELL2"),
    BC1 = c("A", "B", "C", "D"),
    BC2 = c("A", "B", "C", "D"),
    BC3 = c("A", "B", "C", "D"),
    n_corrected = c(10L, 9L, 10L, 10L),
    proportion_reads_in_cell = c(0.5, 0.5, 0.49, 0.5),
    stringsAsFactors = FALSE
  )

  reads_thresholds <- list(
    min_reads = 10L,
    min_reads_proportion = 0.5
  )

  result <- filter_assigned_barcodes(
    cell_to_assigned_BCs = cell_to_assigned_BCs,
    reads_thresholds = reads_thresholds,
    verbose = FALSE
  )

  # Only rows meeting BOTH conditions (>= thresholds) should remain
  # Row 1: n_corrected == 10, proportion == 0.5  -> kept (boundary)
  # Row 2: n_corrected == 9                     -> dropped
  # Row 3: proportion == 0.49                   -> dropped
  # Row 4: n_corrected == 10, proportion == 0.5 -> kept (boundary)

  expect_equal(nrow(result), 2L)

  expect_equal(
    result$cell_barcode,
    c("CELL1", "CELL2")
  )

  expect_equal(
    result$BC1,
    c("A", "D")
  )
})


