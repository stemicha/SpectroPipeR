#' ggsave pdf & png
#'
#' @param filename filename
#' @param plot plot
#' @param limitsize limitsize
#' @param width width
#' @param height height
#' @param dpi dpi
#'
#' @return saves ggplots as pdf and png
#' @keywords internal
#' @noRd
#'
ggsave_pdf_png <- function(filename = filename,
                           plot = plot,
                           limitsize = F,
                           width = width,
                           height = height,
                           dpi = 300){
  #png
  ggsave(filename = paste(filename,".png",sep = ""),
         plot = plot,
         limitsize = limitsize,
         width = width,
         height = height,
         dpi = dpi,
         device="png")
  #pdf
  ggsave(filename = paste(filename,".pdf",sep = ""),
         plot = plot,
         limitsize = limitsize,
         width = width,
         height = height,
         device="pdf")
}
