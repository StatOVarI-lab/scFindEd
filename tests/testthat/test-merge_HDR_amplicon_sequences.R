test_that("import_merge_HDR_data returns a list with names HDR_df and design_df", {

  fake_HDR <- data.frame(AMPLICON_NAME = "amp1", stringsAsFactors = FALSE)
  fake_design <- data.frame(AmpID = "amp1", `Amplicon sequence` = "ATCG",
                            stringsAsFactors = FALSE)

  mockery::stub(
    import_merge_HDR_data,
    "import_HDR_barcoded_data",
    fake_HDR
  )

  mockery::stub(
    import_merge_HDR_data,
    "import_design_summary_info",
    fake_design
  )

  # Pretend columns are fine
  mockery::stub(
    import_merge_HDR_data,
    "check_required_columns",
    TRUE
  )

  result <- import_merge_HDR_data(
    source_dir = "dummy",
    HDR_barcoded_file_name = "hdr.csv",
    design_summary_file_name = "design.csv",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("HDR_df", "design_df"))
})

test_that("import_merge_HDR_data aborts when required columns are missing", {

  # Minimal fake data frames (content doesn't matter because we force the check to fail)
  fake_HDR <- data.frame(AMPLICON_NAME = "amp1", stringsAsFactors = FALSE)
  fake_design <- data.frame(AmpID = "amp1", `Amplicon sequence` = "ATCG",
                            stringsAsFactors = FALSE)

  mockery::stub(
    import_merge_HDR_data,
    "import_HDR_barcoded_data",
    fake_HDR
  )

  mockery::stub(
    import_merge_HDR_data,
    "import_design_summary_info",
    fake_design
  )

  # Force the required-column check to fail
  mockery::stub(
    import_merge_HDR_data,
    "check_required_columns",
    FALSE
  )

  expect_error(
    import_merge_HDR_data(
      source_dir = "dummy",
      HDR_barcoded_file_name = "hdr.csv",
      design_summary_file_name = "design.csv",
      verbose = FALSE
    ),
    "Missing columns detected"
  )
})

test_that("amplicon_data_exploration returns unique amplicons from HDR and design data", {

  HDR_df <- data.frame(
    AMPLICON_NAME = c("amp1", "amp2", "amp1"),
    stringsAsFactors = FALSE
  )

  design_df <- data.frame(
    AmpID = c("amp2", "amp3", "amp3"),
    stringsAsFactors = FALSE
  )

  result <- amplicon_data_exploration(
    HDR_df = HDR_df,
    design_df = design_df,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("unique_HDR_amplicons", "unique_design_amplicons"))

  expect_equal(result$unique_HDR_amplicons, unique(HDR_df$AMPLICON_NAME))
  expect_equal(result$unique_design_amplicons, unique(design_df$AmpID))
})

test_that("amplicon_matching_analysis returns intersection of HDR and design amplicons", {

  unique_HDR_amplicons <- c("amp1", "amp2", "amp3")
  unique_design_amplicons <- c("amp2", "amp3", "amp4")

  result <- amplicon_matching_analysis(
    unique_HDR_amplicons = unique_HDR_amplicons,
    unique_design_amplicons = unique_design_amplicons,
    verbose = FALSE
  )

  expect_equal(result, c("amp2", "amp3"))
})

test_that("perform_amplicon_merge left-joins amplicon sequence onto HDR by AMPLICON_NAME", {

  HDR_df <- data.frame(
    AMPLICON_NAME = c("amp1", "amp2", "amp_missing"),
    HDR_value = c(10, 20, 30),
    stringsAsFactors = FALSE
  )

  design_df <- data.frame(
    AmpID = c("amp1", "amp2"),
    Name = c("ignored", "ignored"),
    `Amplicon sequence` = c("SEQ1", "SEQ2"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  merged_df <- perform_amplicon_merge(
    HDR_df = HDR_df,
    design_df = design_df,
    verbose = FALSE
  )

  # Left join keeps all HDR rows
  expect_equal(nrow(merged_df), nrow(HDR_df))

  # Column from design was added
  expect_true("Amplicon sequence" %in% names(merged_df))

  # Matching rows get sequences
  expect_equal(
    merged_df$`Amplicon sequence`[merged_df$AMPLICON_NAME == "amp1"],
    "SEQ1"
  )
  expect_equal(
    merged_df$`Amplicon sequence`[merged_df$AMPLICON_NAME == "amp2"],
    "SEQ2"
  )

  # Non-matching rows get NA
  expect_true(
    is.na(merged_df$`Amplicon sequence`[merged_df$AMPLICON_NAME == "amp_missing"])
  )
})

test_that("export_merged_amplicon_data writes CSV file when output_dir is provided", {

  merged_df <- data.frame(
    AMPLICON_NAME = c("amp1", "amp2"),
    `Amplicon sequence` = c("SEQ1", "SEQ2"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tmp_dir <- tempfile("amplicon_export_")

  export_merged_amplicon_data(
    merged_df = merged_df,
    output_dir = tmp_dir,
    verbose = FALSE
  )

  output_file <- file.path(
    tmp_dir,
    "HDR_file_barcoded_with_amplicon_sequences.csv"
  )

  expect_true(dir.exists(tmp_dir))
  expect_true(file.exists(output_file))
})


