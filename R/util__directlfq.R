

#' Run directlfq
#'
#' This function perform directlfq via the [directlfq]("https://github.com/MannLabs/directlfq/") module.
#'
#' @param df A dataframe in the generic input format of directlfq with the following columns: "protein", "ion", runs.
#' @param ncores Use howmany cores, default is 4. If it is 0, will use all available threads.
#' @param tempdir The path of temp files.
#' @param reformatted If TRUE means df is in generic input format, else df is a result table from quantity.
#' We recommend use a formatted input
#' @param ... To be completed in the future. Other parameters passed to directlfq.
#'
#' @details
#' Specifically, this R package installs a conda environment containing the directlfq library using reticulate.
#' Then, directlfq calculations are performed by calling Python scripts within that environment.
#' This function will try to
#' Or, you can set up the required conda environment by running the check_directlfq_Module() function.
#' If you are a Windows user, you may need to run this function with administrator privileges.
#' The directlfq function will also attempt to deploy this conda environment if it is missing.
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr reticulate
#' @examples
#' ## A diann report.tsv
#' data(example_diann_res)
#' pro_int = directlfq(example_diann_res,reformatted = F)
#'
#' ## A generic input table, which has been reformatted.
#' data(example_generic_input)
#' pro_int = directlfq(example_generic_input, reformatted = T)
#'
directlfq <- function(Spectronaut_file, ncores = 4, temp_dir = tempdir(), ...){
  if(check_directlfq_Module()){
    cat("Confirmation passed, the directlfq py module is verified as installed in your environment.\n")

    file.copy(from = Spectronaut_file,to = file.path(temp_dir, "temp_input_directlfq.tsv"),overwrite = T)

    in_file <-  file.path(temp_dir, "temp_input_directlfq.tsv")

    py_path = reticulate::conda_list() %>%
      filter(.data$name == "directlfq") %>% select(.data$python)
    cat(paste("python path is:",py_path,"\n"))

    # replace intable_config.yaml with file provided in SpectroPipeR inst/extdata folder
    # Extract the root directory of the conda environment based on the OS
    if ( Sys.info()["sysname"] == "Windows") {
      # Windows uses backslashes, so we'll handle paths with backslashes
      directLFQ_root_path <- sub("\\\\bin\\\\python$", "", py_path)  # Double backslashes for escaping
      directLFQ_root_path <- gsub("\\\\", "/", py_path)  # Convert to forward slashes for consistency
    } else {
      # For Linux and macOS (which both use forward slashes)
      directLFQ_root_path <- sub("/bin/python$", "", py_path)
    }


    directLFQ_path_config_yaml <- file.path(directLFQ_root_path,
                                            "lib",
                                            "python3.8",
                                            "site-packages",
                                            "directlfq",
                                            "configs",
                                            "intable_config.yaml")

    # overwrite config / directLFQ shoud use column names provided by SpectroPipeR report scheme
    try(file.copy(from = system.file("extdata",
                                     "run_directLFQ_script__intable_config.yaml",
                                     package="SpectroPipeR"),
                  to = directLFQ_path_config_yaml,
                  overwrite = T),silent = T)



    script_path <-  system.file("extdata","run_directLFQ_script.py", package = "SpectroPipeR")
    script_path <-  paste0('"',script_path,'"')
    command <-  paste(py_path, script_path, in_file, ncores)
    res = system(command)

    if(res == 0){
      output_directLFQ_calc <-  list(
        directLFQ_protein =  read_tsv(paste0(in_file,".protein_intensities.tsv"),show_col_types = FALSE),
        directLFQ_norm_ion =  read_tsv(paste0(in_file,".ion_intensities.tsv"),show_col_types = FALSE)
      )
      message("Finnished!")
      return(output_directLFQ_calc)
    }else{
      message("Encountered error when running directlfq. Please check your input.")
    }
  }else{
    message("Interrupted running because directlfq is not correctly installed. Run check_directlfq_Module to install the module.")
  }

  return(NULL)
}



