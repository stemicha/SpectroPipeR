

#' Check the directlfq environment
#'
#' This function check and install (if absent) a python conda environment contains the directlfq library
#'
#' @param confirm_or_not Whether provide a second comfirm.
#'
#' @return
#' TURE if the environment is fine, and FALSE if the environment is missed.
#' @keywords internal
#' @noRd
check_directlfq_Module <- function(confirm_or_not = T){
  conda_table <- try(reticulate::conda_list(),silent = T)
  if(class(conda_table) == "try-error"){
    if(grepl("Unable to find conda",as.character(conda_table))){
      install_conda()
    }

    install_condares = install_directlfq_module(confirm_or_not)
  }

  conda_table <- try(reticulate::conda_list(),silent = T)
  if(class(conda_table) == "try-error"){
    cat("Conda is still missing in your system. Please confirm the config in the reticulate package.\n",
        "In windows, installation of software reuqired R needs to be executed with administrator privileges.\n",
        "If you have any uncertainties about how to proceed, please consult the author...")
    return(F)
  }

  if("directlfq" %in% conda_table$name){
    module_table = reticulate::py_list_packages(envname = "directlfq")
    if("directlfq" %in% module_table$package){
      return(TRUE)
    }else{
      message("The directlfq in absent in the conda environment. Installing via pip....")
      reticulate::py_install("directlfq",envname = "directlfq", pip = T)
    }
  }else{
    message("The directlfq environment is absent.\nCreacting the conda environment contains python 3.8.")
    reticulate::conda_create(envname = "directlfq", python_version = 3.8)
    message("Installing directlfq via pip")
    reticulate::py_install("directlfq",envname = "directlfq", pip = T)
  }
}



install_conda <- function(){
  confirm <- confirm("The conda can't be found in the path, which is essential for directlfq module.\nWould your like to install miniconda for python? ")

  if(confirm){
    file.remove(list.files(tempdir(),pattern = "^Miniconda.*.exe$",full.names = T))
    installres = try(reticulate::install_miniconda(force = T),silent = F)
    try(reticulate::conda_run2(cmd_line = "conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main"),silent = F)
    try(reticulate::conda_run2(cmd_line = "conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r"),silent = F)
    try(reticulate::conda_run2(cmd_line = "conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2"),silent = F)
    if(class(installres) != "try-error"){
      return(T)
    }
  }else{
    message("Refuse to install")
    return(F)
  }

  message("Installation failed")
  file.remove(list.files(tempdir(),pattern = "^Miniconda.*.exe$",full.names = T))
  return(F)
}

install_directlfq_module <- function(confirm_or_not){
  if(confirm_or_not){
    confirm <- confirm("The directlfq environment is absent.\nCreact the environment now")
  }else{
    confirm = T
  }
  options(timeout = 240)
  reticulate::conda_create(envname = "directlfq", python_version = 3.8)
}

confirm <- function(message) {
  while (TRUE) {
    cat(message)
    response <- readline(paste("(Y/n): "))
    if (toupper(response) == "Y") {
      return(TRUE)
    } else if (toupper(response) == "N") {
      return(FALSE)
    } else if (response == ""){
      return(TRUE)
    }else {
      cat("Invalid response. Please enter Y or N.\n")
    }
  }
}
