#' ggsave pdf & png
#'
#' @param filename filename
#' @param plot plot
#' @param limitsize limitsize
#' @param number_of_conditions number of condition add 2 inch /// width +  ((number_of_conditions %/% 15) * 2)
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
                           number_of_conditions = NULL,
                           width = width,
                           height = height,
                           dpi = 300){

  if(is.null(number_of_conditions)){
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
  }else{
    #png
    ggsave(filename = paste(filename,".png",sep = ""),
           plot = plot,
           limitsize = limitsize,
           width = width +  ((number_of_conditions %/% 15) * 2),
           height = height,
           dpi = dpi,
           device="png")
    #pdf
    ggsave(filename = paste(filename,".pdf",sep = ""),
           plot = plot,
           limitsize = limitsize,
           width = width + ((number_of_conditions %/% 15) * 2),
           height = height,
           device="pdf")

  }

}
