#' SpectroPipeR onload function
#'
#' @param libname libname
#' @param pkgname pkgname
#'
#' @keywords internal
#' @noRd
#'
#'
.onAttach <- function(libname, pkgname) {
  # package loading startup message
  #if no quarto CLI are found
  if(Sys.which("quarto")==""){
    packageStartupMessage("
#################   Welcome to   ####################
 _____             _           _____ _         _____
|   __|___ ___ ___| |_ ___ ___|  _  |_|___ ___| __  |
|__   | . | -_|  _|  _|  _| . |   __| | . | -_|    -|
|_____|  _|___|___|_| |_| |___|__|  |_|  _|___|__|__|
      |_|                             |_|
#####################################################")
    packageStartupMessage(paste("version:",packageVersion("SpectroPipeR")))
    packageStartupMessage("please cite: https://doi.org/10.1093/bioinformatics/btaf086")
    packageStartupMessage("!!! no quarto CLI found please install and restart !!!")
    packageStartupMessage("https://quarto.org/docs/get-started/")
    browseURL("https://quarto.org/docs/get-started/")

  }else{

    packageStartupMessage("
#################   Welcome to  #####################
 _____             _           _____ _         _____
|   __|___ ___ ___| |_ ___ ___|  _  |_|___ ___| __  |
|__   | . | -_|  _|  _|  _| . |   __| | . | -_|    -|
|_____|  _|___|___|_| |_| |___|__|  |_|  _|___|__|__|
      |_|                             |_|
#####################################################")
    packageStartupMessage(paste("version:",packageVersion("SpectroPipeR")))
    packageStartupMessage("please cite: https://doi.org/10.1093/bioinformatics/btaf086")
  }


}
