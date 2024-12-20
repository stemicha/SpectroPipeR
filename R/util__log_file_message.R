#' message function that also writes log file
#'
#' @param text text
#' @param color color
#' @param log_file_name log file path
#'
#' @return messages and write to log file
#' @keywords internal
#' @noRd
#'
message_function <- function(text = "", color = "blue", log_file_name = ""){

  if(color == "blue"){
    message(crayon::blue(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " MESSAGE: ",text),file = log_file_name,append = T)
  }
  if(color == "red"){
    message(crayon::red(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " ERROR: ",text),file = log_file_name,append = T)
  }
  if(color == "green"){
    message(crayon::green(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " MESSAGE: ",text),file = log_file_name,append = T)
  }
  if(color == "yellow"){
    message(crayon::yellow(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " WARNING: ",text),file = log_file_name,append = T)
  }
  if(color == "magenta"){
    message(crayon::magenta(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " MESSAGE: ",text),file = log_file_name,append = T)
  }
  if(color == "cyan"){
    message(crayon::cyan(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " MESSAGE: ",text),file = log_file_name,append = T)
  }
  if(color == "white"){
    message(crayon::white(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " MESSAGE: ",text),file = log_file_name,append = T)
  }
  if(color == "silver"){
    message(crayon::silver(text))
    write(x = paste0(format(Sys.time(), "%d.%m.%Y %X"), " MESSAGE: ",text),file = log_file_name,append = T)
  }
}
