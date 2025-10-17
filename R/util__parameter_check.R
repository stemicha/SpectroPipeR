#' Parameter list object input cheking function
#'
#' @param parameter provide paramter list object for check parameter inputs
#' @param log_file_name logfile path
#' @keywords internal
#' @noRd
#'
parameter_check <- function(parameter = NULL, log_file_name = NULL){

  message(crayon::blue("checking parameters ..."))

  #ion_q_value_cutoff check
  if(as.character(class(parameter$ion_q_value_cutoff))!="numeric"){
    message_function(text = "ion_q_value_cutoff: value is not numeric",color = "red",log_file_name = log_file_name)
    stop()
  }else{
    if(parameter$ion_q_value_cutoff>0.01){
      message_function(text = "ion_q_value_cutoff bigger than 0.01 !!! filter criteria might be to relaxed !!!",color = "yellow",log_file_name = log_file_name)
    }
  }

  #id_drop_cutoff check
  if(as.character(class(parameter$id_drop_cutoff))!="numeric"){
    message_function(text = "id_drop_cutoff: value is not numeric",color = "red",log_file_name = log_file_name)
    stop()
  }

  #normalization_factor_cutoff_outlier check
  if(as.character(class(parameter$normalization_factor_cutoff_outlier))!="numeric"){
    message_function(text = "normalization_factor_cutoff_outlier: value is not numeric",color = "red",log_file_name = log_file_name)
    stop()
  }

  #filter_oxidized_peptides check
  if(as.character(class(parameter$filter_oxidized_peptides))!="logical"){
    message_function(text = "filter_oxidized_peptides: value is not logical",color = "red",log_file_name = log_file_name)
    stop()
  }

  #stat_test check
  if(sum(parameter$stat_test=="rots" | parameter$stat_test=="modt" | parameter$stat_test=="t")==0){
    message_function(text = "stat_test is not set to: rots or modt or t",color = "red",log_file_name = log_file_name)
    stop()
  }

  #type_slr check
  if(sum(parameter$type_slr=="median" | parameter$type_slr=="tukey")==0){
    message_function(text = "type_slr is not set to: median or tukey",color = "red",log_file_name = log_file_name)
    stop()
  }

  #fold_change check
  if(as.character(class(parameter$fold_change))!="numeric"){
    message_function(text = "fold_change: value is not numeric",color = "red",log_file_name = log_file_name)
    stop()
  }

  #fold_change check
  if(as.character(class(parameter$p_value_cutoff))!="numeric"){
    message_function(text = "p_value_cutoff: value is not numeric",color = "red",log_file_name = log_file_name)
    stop()
  }else{
    if(parameter$p_value_cutoff>0.05){
      message_function(text = "p_value_cutoff bigger than 0.05 !!! filter criteria might be to relaxed !!!",color = "yellow",log_file_name = log_file_name)
      }
  }

  #protein intensity calculation
  if(as.character(class(parameter$protein_intensity_estimation))!="character"){
    message_function(text = "protein_intensity_estimation: value is not character",color = "red",log_file_name = log_file_name)
    stop()
  }else{
    if(sum(is.element(el = c("Hi3","MaxLFQ", "directLFQ"),set = parameter$protein_intensity_estimation))!=1){
      message_function(text = "protein_intensity_estimation is not one of these options !!! character: 'Hi3' = Hi3 protein intensity estimation OR 'MaxLFQ' = MaxLFQ protein intensity estimation OR 'directLFQ' = directLFQ protein intensity estimation !!!",color = "red",log_file_name = log_file_name)
      stop()
    }
  }

}
