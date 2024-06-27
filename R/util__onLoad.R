#' SpectroPipeR onload function
#'
#' @param libname libname
#' @param pkgname pkgname
#'
#' @keywords internal
#' @noRd
#'
#'
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.SpectroPipeR <- list(
    ragg.max_dim = 5E10,
    ggrepel.max.overlaps = 100
  )
  toset <- !(names(op.SpectroPipeR) %in% names(op))
  if (any(toset)) options(op.SpectroPipeR[toset])

  #overwrite dplyr.summarise.inform
  options(dplyr.summarise.inform = FALSE)
  #set global variables
  utils::globalVariables(c("i"))
  # hide GGally overwrite method
  Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
  invisible()
}
