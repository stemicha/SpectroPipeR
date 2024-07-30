#' figure width estimation on condition and sample length
#'
#' @param sample_length number of samples
#' @param condition_length number of conditons
#'
#' @return messages and write to log file
#' @keywords internal
#' @noRd
#'

fig_width_estimation <- function(sample_length = 717, condition_length = 4){
  if((0.4*sample_length)<10){
    return(15)
    }else{
      if(sample_length/condition_length<5){
        return(4*condition_length)
      }else{
        return(0.4*sample_length)
      }
    }
}
