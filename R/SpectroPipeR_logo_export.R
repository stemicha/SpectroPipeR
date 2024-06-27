#' SpectroPipeR - logo export
#'
#' @param output_location folder location where the SpectroPipeR logo should be copied
#'
#' @return
#' copies the SpectroPipeR logo to the designated folder
#'
#' @export
#'
#' @examples
#'
#' SpectroPipeR_logo_export(output_location = "../SpectroPipeR_test_folder")
#'
SpectroPipeR_logo_export <- function(output_location = ""){

  #create output dir
  if(!dir.exists(output_location)){
    dir.create(paste0(output_location),recursive = T,showWarnings = F)
  }

  # copy Spectronaut report scheme
  file.copy(from = system.file("extdata", "SpectroPipeR_hexbin_logo.png", package="SpectroPipeR"),
            to = paste0(output_location,"/","SpectroPipeR_hexbin_logo.png"),
            overwrite = T)
  message(paste("The SpectroPipeR_hexbin_logo.png file was sucessfully copied to",output_location))

}
