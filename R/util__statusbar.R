#' statusprogressbar --> statusbar from https://github.com/STATWORX/helfRlein
#'
#' @param run index of run
#' @param max.run max run
#' @param width width
#' @param info info character
#' @param percent.max max percentage
#'
#' @author Jakob Gepp
#' @keywords internal
#' @noRd
#'
statusprogressbar <- function(run,
                      max.run,
                      width = 20L,
                      info = run,
                      percent.max = width) {
  # check for old parameter
  if ("percent.max" %in% names(match.call())) {
    warning("'percent.max' is deprecated, please use 'width' instead")
    width <- percent.max
  }

  # check run
  if (length(run) > 1) {
    stop("run needs to be of length one!")
  }
  # check max.run
  if (length(max.run) == 0) {
    stop("max.run has length 0")
  }

  if (length(max.run) > 1 || is.character(max.run)) {
    percent <- which(run == max.run) / length(max.run)
  } else {
    percent <- run / max.run
  }

  percent_step <- round(percent * width, 0)
  progress <- paste0("[",
                     paste0(rep("=", percent_step), collapse = ""),
                     paste0(rep(" ", width - percent_step),
                            collapse = ""),
                     "] ",
                     sprintf("%7.2f", percent * 100),
                     "% - ",
                     info)
  last_step <- ifelse(run == max.run[length(max.run)], "\n", "")
  cat("\r", progress, last_step)
  flush.console()
}
