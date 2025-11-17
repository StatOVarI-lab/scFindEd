#' @keywords internal
#' Corrects a barcode by selecting the closest match, as long as the closest
#' match has its distance below a defined threshold.
correct_BC <- function(distances, thresh){
  # If the distance < the threshold, select the barcode as a match
  if (distances[1] < thresh){
    output <- names(distances)[1]
  }
  # Otherwise, return no output
  else {
    output <- ""
  }
  return(output)
}
