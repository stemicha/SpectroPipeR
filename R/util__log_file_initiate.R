#' initiate log file
#'
#' @param file file path
#'
#' @keywords internal
#' @noRd
#'
init_log_file <- function(file = ""){
  write("
 _____             _           _____ _         _____
|   __|___ ___ ___| |_ ___ ___|  _  |_|___ ___| __  |
|__   | . | -_|  _|  _|  _| . |   __| | . | -_|    -|
|_____|  _|___|___|_| |_| |___|__|  |_|  _|___|__|__|
      |_|                             |_|
##################  log file  #######################
",file = file)
  write(x = paste("version:",packageVersion("SpectroPipeR")),file = file,append = T)
  write(x = " ",file = file,append = T)
}
