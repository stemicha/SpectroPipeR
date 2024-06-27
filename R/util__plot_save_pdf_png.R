#' save pdf & png
#'
#' @param filename filename
#' @param plot plot
#' @param width width
#' @param height height
#'
#' @return saves plots as png
#' @keywords internal
#' @noRd
#'
plot_save_png <- function(filename = filename,
                           plot = plot,
                           width = width,
                           height = height){
  #png
  png(paste0(filename,".png"),
      width = width,
      height = height,
      units = "in",
      res = 200)
  plot
  dev.off()

}
