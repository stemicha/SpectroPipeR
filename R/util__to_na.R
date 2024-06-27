#' Replace NaN and Inf with NA from https://github.com/STATWORX/helfRlein
#'
#' @param x vector
#'
#' @author Daniel Luettgau
#'
#' @keywords internal
#' @noRd
to_na <- function(x) {
  # check input
  if (!is.vector(x)) {
    stop("input must be a vector")
  }

  # convert input
  if (is.character(x)) {
    return(x)
  } else {
    ifelse(is.infinite(x) | is.nan(x), NA, x)
  }
}
