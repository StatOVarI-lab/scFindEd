#' Say hello to someone
#'
#' This function greets a person by name.
#' It uses a small helper function to format the greeting.
#'
#' @param name A character string with the person's name.
#' @return A greeting message.
#' @examples
#' say_hello("Rux")
#' @export
say_hello <- function(name) {
  message <- make_greeting(name)  # call helper
  return(message)
}

#' @keywords internal
make_greeting <- function(name) {
  paste("Hello,", name, "!")
}
