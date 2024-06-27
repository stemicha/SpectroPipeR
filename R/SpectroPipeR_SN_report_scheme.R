#' SpectroPipeR - Spectronaut report export scheme
#'
#' @param output_location folder location where the SpectroPipeR_report.rs report export scheme should be copied
#'
#' @return
#'
#' Generates a Spectronaut report export template that can be imported into Spectronaut, encompassing all necessary columns SpectroPipeR needs for the Spectronaut export of the analyzed data.
#' @export
#'
#' @examples
#'
#' Spectronaut_export_scheme(output_location = "../SpectroPipeR_test_folder")
#'
Spectronaut_export_scheme <- function(output_location = ""){

  #create output dir
  if(!dir.exists(output_location)){
    dir.create(paste0(output_location),recursive = T,showWarnings = F)
  }

  # copy Spectronaut report scheme
  file.copy(from = system.file("extdata", "SpectroPipeR_report.rs", package="SpectroPipeR"),
            to = paste0(output_location,"/","SpectroPipeR_report.rs"),
            overwrite = T)
  message(paste("The SpectroPipeR_report.rs file was sucessfully copied to",output_location))

}
