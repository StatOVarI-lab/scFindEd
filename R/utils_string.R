#' @keywords internal
#' Takes a list of data frame names with vectors of column names to check for in
#' each, and returns TRUE/FALSE dependent on whether every required column is
#' present.
check_required_columns <- function(
  required_cols
) {
  columns_present = TRUE

  for (df_name in names(required_cols)) {
    check_df <- get(df_name, envir = parent.frame())
    cols_to_check <- required_cols[[df_name]]

    # Finds columns present in required_cols but not in check_df
    missing_cols <- setdiff(cols_to_check, colnames(check_df))

    # Update the user if columns are missing
    if(length(missing_cols) > 0) {
      message(glue::glue("Missing column(s) in {df_name}: {paste(missing_cols, collapse = ', ')}"))
      columns_present = FALSE
    }
  }
  return(columns_present)
}
