#' nin function from https://github.com/STATWORX/helfRlein
#'
#' @param x x
#' @param table table
#'
#' @keywords internal
#' @noRd
#'
`%nin%` <- function(x, table) {
  !match(x, table, nomatch = 0) > 0
}
