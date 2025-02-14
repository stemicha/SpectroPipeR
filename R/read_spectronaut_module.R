#' SpectroPipeR: loading Spectronaut data module
#' @description
#' Function for loading Spectronaut data and performing identification (ID) analysis, which is an essential first step in the SpectroPipeR workflow.
#'
#' @param file location (path) of Spectronaut output report;
#'   you should use the `Spectronaut_export_scheme()` function for getting a SpectroPipeR report scheme encompassing all mandatory columns
#' @param ID_condition_filtering TRUE or FALSE if a condition-wise filtering should be performed
#' @param ID_condition_filtering_percent (numerical value ranging from 0 - 1, default = 0.5) define the proportion for the condition-wise ID filtering
#' @param max_chars_file_name_capping integer, (default = 35) number of max characters used for raw file name presentation; must be adjusted if function
#' @param print.plot if TRUE --> printing ID plot on ion level coloring corresponds ID outlier estimate by "id_drop_cutoff" variable
#' @param report_copy if TRUE --> copy Spectronaut input report to SpectroPipeR project folder 01_input_data
#' @param parameter __mandatory parameter list element__
#'
#'  _table of list elements:_
#'
#' | <u> __parameter__ </u> | <u> __description__ </u>                |
#' |:-----------------------|:--------------------------------------|
#' | output_folder          | **mandatory !!!** - _character_ - output folder path (abs.) |
#' | ion_q_value_cutoff     | **default = 0.01** - _numeric_ - Q-value used in Spectronaut analysis: Biognosys |
#' |                        | default is 0.01 = 1% error rate |
#' | id_drop_cutoff         | **default = 0.3** - _numeric_ - value between 0-1 (1 = 100%); xx percent lower |
#' |                        | than median of ion ID rate => outlier |
#' | normalization_method   | **default = "median"** - _character_ - "median" or Spectronaut - auto-detection |
#' |                        | is per default ON, meaning if normalization was performed in Spectronaut |
#' |                        | this will be detected and preferred over parameter setting here;|
#' |                        | median normalization is the fallback option|
#' | normalization_factor_cutoff_outlier | **default = 4** - _numeric_ - median off from global median |
#' |                        | (4 means abs. 4fold off) |
#' | filter_oxidized_peptides | **default = TRUE** _logical_ - if oxidized peptides should be removed before |
#' |                          |peptide quantification |
#' | protein_intensity_estimation | **default = "MaxLFQ"** - _character_ - Hi3 = Hi3 protein intensity estimation, |
#' |                              |  MaxLFQ = MaxLFQ protein intensity estimation |
#' | stat_test              | **default = "rots"** - _character_ - choose statistical test: "rots" = reproducibility |
#' |                        | optimized test statistics, "modt" = moderate t-test (lmfit, eBayes),|
#' |                        | "t" = t-test |
#' | type_slr               | **default = "median"** - _character_ - choose ratio aggregation method: |
#' |                        | "median" or "tukey" is used when calculating protein values |
#' | fold_change            | **default = 1.5** - _numeric_ - fold-change used as cutoff e.g. 1.5 |
#' | p_value_cutoff         | **default = 0.05** - _numeric_ - p-value used as cutoff e.g. 0.05 |
#' | paired                 | **default = FALSE** - _logical_ - Should paired statistics be applied? |
#'
#'
#' | <u>example parameters list (default)</u>:                         |
#' |---------------------------------------------------------|
#' | params <- list(output_folder = "../Spectronaut_example",|
#' |               ion_q_value_cutoff = 0.01,|
#' |               id_drop_cutoff = 0.3,|
#' |               normalization_method = "median",|
#' |               normalization_factor_cutoff_outlier = 4,|
#' |               filter_oxidized_peptides = T,|
#' |               protein_intensity_estimation = "MaxLFQ",|
#' |               stat_test = "rots",|
#' |               type_slr = "median",|
#' |               fold_change = 1.5,|
#' |               p_value_cutoff = 0.05,|
#' |               paired = FALSE|
#' |              )|
#'
#'
#' @returns SpectroPipeR_data list object with the loaded raw data and processed data tables, in addition to the automatically saved tables and plots.
#'  For the description of the generated figures and tables please read the manual & vignettes
#'
#' | <u> __list element__ </u> | <u> __description__ </u>              |
#' |:--------------------------|:--------------------------------------|
#' | spectronaut_output        | *tibble:* Spectronaut report tibble provided for the analysis |
#' | SDRF_file                 | *tibble:* intermediate SDRF table of the analysis |
#' | summary_distinct          | *tibble:* distinct ion, modified peptide, stripped peptides and |
#' |                           | protein group count per file filtered by provided Q-value |
#' | raw_file_names            | *tibble:* R.FileNames capped and uncapped version together with |
#' |                           | R.Condition and R.Replicate |
#' | ion_id_median             | *numerical value:* median of ion intensity |
#' | ion_id_cutoff             | *numerical value:* ion ID count threshold to classify sample as outlier |
#' | PG_2_peptides_ID_raw      | *tibble:* with protein groups with at least 2 peptides with peptide
#' |                           | and replicate count |
#' | summary_distinct_outlier  | *tibble:* if outlier are detected they are listed in this tibble |
#' | ID_rate_plot              | *ggplot2 plot:* ID rate plot |
#' | ID_rate_plot_filter       | *ggplot2 plot:* ion ID rate plot with ion ID cutoff line |
#' | sample_length             | *numberical value:* number of samples in the provided Spectronaut report |
#' | parameter                 | *list:* parameters provided by the user |
#' | time_stamp_log_file       | *string:* time stamp of the log file (format: %Y_%m_%d__%H_%M) |
#' | log_file_name             | *string:* analysis log file name |
#'
#'
#' @md
#' @import dplyr
#' @import ggplot2
#' @import tidyverse
#' @import readr
#' @import tidyr
#' @import grDevices
#' @import utils
#' @import openxlsx
#' @import forcats
#' @import patchwork
#' @importFrom stats median sd lm as.formula quantile
#' @export
#'
#' @examples
#' \donttest{
#'#load library
#'library(SpectroPipeR)
#'
#'# use default parameters list
#'params <- list(output_folder = "../SpectroPipeR_test_folder")
#'
#'# example input file
#'example_file_path <- system.file("extdata",
#'                                 "SN_test_HYE_mix_file.tsv",
#'                                 package="SpectroPipeR")
#'
#'# step 1: load Spectronaut data module
#'SpectroPipeR_data <- read_spectronaut_module(file = example_file_path,
#'                                             parameter = params,
#'                                             print.plot = FALSE)
#'}

read_spectronaut_module <- function(file = "",
                                    ID_condition_filtering = FALSE,
                                    ID_condition_filtering_percent = 0.5,
                                    parameter = list(),
                                    max_chars_file_name_capping = 35,
                                    print.plot = FALSE,
                                    report_copy = F){

# setup default parameters ------------------------------------------------
  parameter_user_input <- parameter

  if(is.null(parameter_user_input$ion_q_value_cutoff)){
    parameter_user_input$ion_q_value_cutoff <- 0.01 # default = 0.01
  }
  if(is.null(parameter_user_input$id_drop_cutoff)){
    parameter_user_input$id_drop_cutoff <- 0.3 # default = 0.4
  }
  if(is.null(parameter_user_input$normalization_method)){
    parameter_user_input$normalization_method <- "median" # default = "median"
  }
  if(is.null(parameter_user_input$normalization_factor_cutoff_outlier)){
    parameter_user_input$normalization_factor_cutoff_outlier <- 4 # default = 4
  }
  if(is.null(parameter_user_input$filter_oxidized_peptides)){
    parameter_user_input$filter_oxidized_peptides <- TRUE # default = TRUE
  }
  if(is.null(parameter_user_input$protein_intensity_estimation)){
    parameter_user_input$protein_intensity_estimation <- "MaxLFQ" # default = "MaxLFQ"
  }
  if(is.null(parameter_user_input$stat_test)){
    parameter_user_input$stat_test <- "rots" # default = "rots"
  }
  if(is.null(parameter_user_input$type_slr)){
    parameter_user_input$type_slr <- "median" # default = "median"
  }
  if(is.null(parameter_user_input$fold_change)){
    parameter_user_input$fold_change <- 1.5 # default = 1.5
  }
  if(is.null(parameter_user_input$p_value_cutoff)){
    parameter_user_input$p_value_cutoff <- 0.05 # default = 0.05
  }
  if(is.null(parameter_user_input$paired)){
    parameter_user_input$paired <- FALSE # default = FALSE
  }

  # get output folder main
  out_folder <- parameter_user_input$output_folder

  #define parameter_input
  parameter_input <- parameter_user_input

  #add input file path
  parameter_input$Spectronaut_report_file <- file

  #add package version to parameters
  parameter_input <- c(list(SpectroPipeR_version = packageVersion("SpectroPipeR")),
                       parameter_input)

  #check if output-folder exist
  if (dir.exists(parameter_user_input$output_folder)==FALSE) {
    dir.create(parameter_user_input$output_folder,recursive = T)
    crayon::magenta("generate output folder")
  }else{
    crayon::magenta("provided output folder exists / files and folders will be overwritten")
  }

  # generate log file
  time_stamp_log_file <- format(Sys.time(), "%Y_%m_%d__%H_%M")
  log_file_name <- paste0(out_folder,"/",time_stamp_log_file, "_SpectroPipeR_analysis.log")
  # init log file
  init_log_file(file = log_file_name)
  parameter_input$log_file_name <- log_file_name

  # read input
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "# READ SPECTRONAUT MODULE",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "",
                   color = "white",
                   log_file_name = log_file_name)

  #check parameters
  parameter_check(parameter_input, log_file_name = log_file_name)


  # create folders
  dir.create(paste0(out_folder,"/","01_input_data"),recursive = T,showWarnings = F)
  dir.create(paste0(out_folder,"/","02_ID_rate"),recursive = T,showWarnings = F)

  message_function(text = "loading data ...",
                   color = "blue",
                   log_file_name = log_file_name)

  # fast data input with readr ====
  tmp_data_input <- suppressWarnings(readr::read_delim(file = file,
                               delim = "\t",
                               col_names = T,
                               show_col_types = FALSE,
                               col_types = readr::cols(EG.Qvalue = readr::col_character(),
                                                       PG.IBAQ = readr::col_character())))


  # write input data to project folder ====
  if(report_copy == T){
    # write input file to input data folder
    message_function(text = "write input data to output folder ...",
                     color = "blue",
                     log_file_name = log_file_name)

    readr::write_delim(x = tmp_data_input,
                       file = paste0(out_folder,"/","01_input_data","/",rev(unlist(strsplit(file,split = "/")))[1]),delim = "\t")

  }else{
    # DO not write input file to input data folder
    message_function(text = "write input data to output folder was omitted...",
                     color = "blue",
                     log_file_name = log_file_name)
  }

  # change condition from numeric to character if needed
  if(as.character(class(tmp_data_input$R.Condition))=="numeric"){
    message_function(text = "condition input is NUMERIC --> converting to character",
                     color = "blue",
                     log_file_name = log_file_name)
    tmp_data_input$R.Condition <- as.character(tmp_data_input$R.Condition)
  }

  # perform R.FileName capping ====
  message_function(text = "R.FileName capping ...",
                   color = "blue",
                   log_file_name = log_file_name)
  max_chars = max_chars_file_name_capping
  tmp_data_input$R.FileName_raw <- tmp_data_input$R.FileName

  sample_data_tmp_data_input <- tmp_data_input %>%
    distinct(.data$R.FileName,
             .data$R.FileName_raw,
             .data$R.Condition,
             .data$R.Replicate) %>%
    arrange(.data$R.Condition)

  # capping R.FileName
  sample_data_tmp_data_input <- sample_data_tmp_data_input %>%
    dplyr::rowwise() %>%
    dplyr::mutate(R.FileName = ifelse(nchar(.data$R.FileName)<=max_chars,
                               yes = .data$R.FileName,
                               no =   paste0(substr(x = .data$R.FileName,start = 1,stop = max_chars/2),
                                             "...",
                                             substr(x = .data$R.FileName,
                                                    start = nchar(.data$R.FileName)-max_chars/2,
                                                    stop = nchar(.data$R.FileName))
                               ))
    )

  if(length(unique(sample_data_tmp_data_input$R.FileName))!=length(unique(sample_data_tmp_data_input$R.FileName_raw))){
    # stop if name capping didnt work
    stop_text <- paste("increase character number for R.FileName capping (max_chars_file_name_capping)")
    message_function(text = stop_text,color = "red",log_file_name = log_file_name)
    stop()
  }

  # merge capped R.FileName to data
  tmp_data_input <- dplyr::left_join(tmp_data_input,
                              sample_data_tmp_data_input %>%
                                dplyr::select(.data$R.FileName,
                                       .data$R.FileName_raw) %>%
                                dplyr::rename(R.FileName_capped = .data$R.FileName),
                              by = c("R.FileName" = "R.FileName_raw")) %>%
    dplyr::select(-.data$R.FileName) %>%
    dplyr::rename(R.FileName = .data$R.FileName_capped)




  # add warning if one condition has only one replicate ====
  condition_sample_count<- tmp_data_input %>%
    dplyr::distinct(.data$R.FileName,
                    .data$R.Condition) %>%
    dplyr::group_by(.data$R.Condition) %>%
    dplyr::summarise(sample_count = dplyr::n_distinct(.data$R.FileName)) %>%
    dplyr::ungroup()

  condition_sample_count_n1 <- condition_sample_count %>%
    dplyr::filter(.data$sample_count==1)

  if(dim(condition_sample_count_n1)[1]!=0){
    message_function(text = "input data contains conditions with only 1 replicate !!!",
                     color = "red",
                     log_file_name = log_file_name)
  }

  #number of samples detected
  sample_length <- length(unique(tmp_data_input$R.FileName))

  #create file number specific folder
  if(file.exists(paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis"))){
    message_function(text = paste0("02_ID_rate/",sample_length,"_sample_analysis"," ATTENTION !!! ---> folder already exists - files will be replaced !!!"),
                     color = "magenta",
                     log_file_name = log_file_name)
  }else{
    message_function(text = paste0("create folder: 02_ID_rate/",sample_length,"_sample_analysis"),
                     color = "magenta",
                     log_file_name = log_file_name)
    dir.create(paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis"),recursive = T,showWarnings = F)
  }


  #input cols for check
  data_input_cols <- c("R.FileName",
                       "R.Condition",
                       "R.Replicate",
                       "R.Instrument Name",
                       "R.Raw File Name",
                       "R.MS1 Mass Analyzer",
                       "R.MS2 Mass Analyzer",
                       "R.Run Date",
                       "PG.ProteinGroups",
                       "PG.Organisms",
                       "PG.IBAQ",
                       "PEP.StrippedSequence",
                       "EG.ModifiedPeptide",
                       "PEP.NrOfMissedCleavages",
                       "EG.UserGroup",
                       "EG.Qvalue",
                       "EG.PEP",
                       "EG.Cscore",
                       "EG.NormalizationFactor",
                       "EG.TotalQuantity (Settings)",
                       "EG.SignalToNoise",
                       "EG.Identified",
                       "EG.ApexRT",
                       "EG.IntCorrScore",
                       "EG.DatapointsPerPeak",
                       "EG.DatapointsPerPeak (MS1)",
                       "FG.Charge",
                       "FG.Id",
                       "FG.XICDBID",
                       "FG.LabeledSequence",
                       "FG.ShapeQualityScore",
                       "FG.MS1Quantity",
                       "FG.MS2Quantity",
                       "FG.MS1RawQuantity",
                       "FG.MS2RawQuantity")

  # test if params is a list
  if(as.character(class(parameter)) != "list"){
    stop_text <- paste("input parameter is not a list")
    message_function(text = stop_text,color = "red",log_file_name = log_file_name)
    stop()
  }

  # test if all needed columns are in data_input
  if(sum(colnames(tmp_data_input)%in%data_input_cols)!=length(data_input_cols)){
    stop_text <- paste(paste("essentiell columns are missing in input data:",paste(data_input_cols[data_input_cols%nin%colnames(tmp_data_input)],collapse = "; ")))
    message_function(text = stop_text,color = "red",log_file_name = log_file_name)
    stop()
  }

  # stop if there are missing values
  if(sum(is.na(tmp_data_input$"EG.TotalQuantity (Settings)"))!=0){
    stop_text <- paste("the provided file contains missing values in EG.TotalQuantity; please use a quantification workflow in Spectronaut which either impute or parses the data to omit missing values")
    message_function(text = stop_text,color = "red",log_file_name = log_file_name)
    stop()
    }

  #add ion data
  tmp_data_input <- tmp_data_input %>%
    unite(.data$EG.ModifiedPeptide, .data$FG.Charge,
          sep = "_", remove = F,col = "ion")

  #reformat Q-value
  tmp_data_input <- tmp_data_input %>%
    mutate(EG.Qvalue.formatted = ifelse(test = .data$EG.Qvalue=="Profiled",
                                        yes = 1,
                                        no = .data$EG.Qvalue))

  tmp_data_input$EG.Qvalue.formatted <- to_na(tmp_data_input$EG.Qvalue.formatted)
  tmp_data_input$EG.Qvalue.formatted<-as.numeric(tmp_data_input$EG.Qvalue.formatted)

  # sort according to conditions ====
  tmp_data_input <- tmp_data_input %>%
    mutate(R.FileName = factor(.data$R.FileName,
                               levels = unique(sample_data_tmp_data_input$R.FileName)))


  # intermediate SDRF file generation and export ====
  message_function(text = "generating intermediate SDRF file", color = "blue",log_file_name = log_file_name)

  # extract raw file function
  extract_raw_file<- function(x) {
    vendor_formats <- c("\\.d","\\.raw","\\.wiff","\\.wiff2")
    tmp_x<- unlist(strsplit(x = x,split = "\\\\"))

    tmp_raw_out <- tmp_x[unlist(sapply(X = vendor_formats,
                                       FUN = function(y){
                                         grep(pattern = y,x = tmp_x)
                                       })
    )
    ]
    return(tmp_raw_out)
  }

  # get unique run meta information
  SDRF_file<- tmp_data_input %>%
                  distinct(.data$R.Condition,
                           .data$R.Replicate,
                           .data$`R.Instrument Name`,
                           .data$`R.Raw File Name`,
                           .data$`R.MS1 Mass Analyzer`,
                           .data$`R.MS2 Mass Analyzer`)

  # extract raw file names from path
  SDRF_file <- SDRF_file %>%
    dplyr::mutate(`R.Raw File Name` = extract_raw_file(x = .data$`R.Raw File Name`))

  # add source column --> condition + replicate
  SDRF_file <- SDRF_file %>%
    tidyr::unite(.data$R.Condition,.data$R.Replicate,remove = F,sep = "__",col = "source name")


  # rename columns
  SDRF_file <- SDRF_file %>%
    dplyr::rename(`comment[data file]` = .data$`R.Raw File Name`,
                  `comment[instrument]` = .data$`R.Instrument Name`,
                  `comment[MS1 analyzer type]` = .data$`R.MS1 Mass Analyzer`,
                  `comment[MS2 analyzer type]` = .data$`R.MS2 Mass Analyzer`,
                  `comment[biological replicate]` = .data$R.Replicate) %>%
    dplyr::select(-.data$R.Condition)

  SDRF_file <- dplyr::bind_cols(SDRF_file,
                         tibble::tibble("technology type" = "proteomic profiling by mass spectrometry",
                                        "comment[label]" = "label free sample",
                                        "comment[fraction identifier]" = 1))


  message_function(text = "saving intermediate SDRF file ... please use https://lessdrf.streamlit.app for finalizing/editing file before submitting to a public repository", color = "blue",log_file_name = log_file_name)
  #write SDRF file output
  write_delim(x = SDRF_file,file = paste0(out_folder,"/","SDRF_intermediate_file.tsv"),delim = "\t")


  # filter data ID_condition_filtering ====
  if(ID_condition_filtering==T){
    message_function(text = "perform ion ID condition filtering", color = "green",log_file_name = log_file_name)

    #add to global parameter for reporting
    parameter_input$ID_condition_filtering <- ID_condition_filtering
    parameter_input$ID_condition_filtering_percent <- ID_condition_filtering_percent
    # sample inputs
    total_sample_count<- tmp_data_input %>%
      dplyr::group_by(.data$R.Condition) %>%
      dplyr::summarise(total_sample_count = n_distinct(.data$R.FileName))%>%
      dplyr::ungroup()
    # filter for ions with a Q-value below threshold choosen in parameters
    filtered_sample_count <- tmp_data_input %>%
      dplyr::select( .data$ion,
                     .data$EG.Qvalue.formatted,
                     .data$R.FileName,
                     .data$R.Condition) %>%
      dplyr::filter(.data$EG.Qvalue.formatted<=parameter_input$ion_q_value_cutoff) %>%
      dplyr::group_by(.data$R.Condition,
                      .data$ion) %>%
      dplyr::summarise(sample_count = n_distinct(.data$R.FileName))%>%
      dplyr::ungroup()

    # add total sample count of input
    filtered_sample_count <- dplyr::left_join(filtered_sample_count,
                                       total_sample_count,
                                       by = "R.Condition")
    # calculate percent
    filtered_sample_count <- filtered_sample_count %>%
      mutate(percent = .data$sample_count / .data$total_sample_count)
    selected_filtered_sample_count <- filtered_sample_count %>%
      dplyr::filter(.data$percent >= ID_condition_filtering_percent)

    #filter input data for condition wise percentage ion identification percentage
    tmp_data_input <- tmp_data_input %>%
      dplyr::filter(.data$ion %in% unique(selected_filtered_sample_count$ion))


    message_function(text = "perform ion ID condition filtering: write output in 01_input_data and go on with analysis", color = "green",log_file_name = log_file_name)
    write_csv(x = tmp_data_input,
              file = paste0(out_folder,"/","01_input_data/input_data__filtered_ion_identifications_condition_wise__Qvalue_cutoff_",parameter_input$ion_q_value_cutoff,"__percentage_cutoff_",ID_condition_filtering_percent,".csv"))

  }


  #do short summary after loading

  message_function(text = paste("_________ data set loaded with ... _________"),
                   color = "green",
                   log_file_name = log_file_name)
  message_function(text = paste("number of raw files =", length(unique(tmp_data_input$R.FileName))),
                   color = "green",
                   log_file_name = log_file_name)
  message_function(text = paste("number of conditions =", length(unique(tmp_data_input$R.Condition))),
                   color = "green",
                   log_file_name = log_file_name)
  message_function(text = paste("number of ions without filtering =", length(unique(tmp_data_input$ion))),
                   color = "green",
                   log_file_name = log_file_name)
  message_function(text = paste("number of peptides without filtering =", length(unique(tmp_data_input$PEP.StrippedSequence))),
                   color = "green",
                   log_file_name = log_file_name)
  message_function(text = paste("number of Protein groups without filtering =", length(unique(tmp_data_input$PG.ProteinGroups))),
                   color = "green",
                   log_file_name = log_file_name)


  message_function(text = "count profiled values ...",
                   color = "blue",
                   log_file_name = log_file_name)

  if(sum(tmp_data_input$EG.Qvalue=="Profiled")==0){

    message_function(text = "no profiled ions in data, profiling was omitted in Spectronaut  ...", color = "green",log_file_name = log_file_name)
    tmp_summary_distinct_ions <- tmp_data_input %>%
      group_by(.data$R.FileName,
               .data$R.Condition,
               .data$R.Replicate) %>%
      summarise(distinct_ions_q_value_filtered = n_distinct(.data$ion[.data$EG.Qvalue.formatted < parameter_input$ion_q_value_cutoff]),
                distinct_ions_q_value_larger = n_distinct(.data$ion[.data$EG.Qvalue.formatted > parameter_input$ion_q_value_cutoff | is.na(.data$EG.Qvalue.formatted)])
                ) %>%
      ungroup()
    colnames(tmp_summary_distinct_ions)[4:5] <- c(paste("<",parameter_input$ion_q_value_cutoff,sep=""),
                                                  paste(">",parameter_input$ion_q_value_cutoff,sep=""))

  }else{
    tmp_summary_distinct_ions <- tmp_data_input %>%
      group_by(.data$R.FileName,
               .data$R.Condition,
               .data$R.Replicate) %>%
      summarise(distinct_ions_q_value_filtered = n_distinct(.data$ion[.data$EG.Qvalue.formatted < parameter_input$ion_q_value_cutoff]),
                distinct_ions_q_value_larger = n_distinct(.data$ion[(.data$EG.Qvalue.formatted > parameter_input$ion_q_value_cutoff | is.na(.data$EG.Qvalue.formatted)) & .data$EG.Qvalue != "Profiled"]),
                distinct_ions_profiled = n_distinct(.data$ion[.data$EG.Qvalue == "Profiled"])) %>%
      ungroup()

    colnames(tmp_summary_distinct_ions)[4:6] <- c(paste("<",parameter_input$ion_q_value_cutoff,sep=""),
                                                  paste(">",parameter_input$ion_q_value_cutoff,sep=""),
                                                  "profiled")

  }


  # performing counting of min. sample percentage where an ion was detected
  message_function(text = "performing counting of min. sample percentage where an ion was detected...",
                   color = "blue",
                   log_file_name = log_file_name)
  message_function(text = paste("ion Q-value cutoff <",parameter_input$ion_q_value_cutoff),
                   color = "blue",
                   log_file_name = log_file_name)

  sample_count_detect_tmp<- min(unlist(tmp_data_input %>%
    dplyr::filter(.data$EG.Qvalue.formatted <= parameter_input$ion_q_value_cutoff) %>%
    dplyr::group_by(.data$ion) %>%
    dplyr::summarise(sample_count = dplyr::n_distinct(.data$R.FileName))%>%
    dplyr::ungroup() %>%
    dplyr::select(.data$sample_count)))

  message_function(text = paste(round(sample_count_detect_tmp/sample_length*100,digits = 2),"%",
                                paste("(",sample_count_detect_tmp,"/",sample_length,")",sep=""),
                                "is the min. sample percentage where an ion was detected"),
                   color = "green",
                   log_file_name = log_file_name)


  message_function(text = "performing ID rate filtering ...",
                   color = "blue",
                   log_file_name = log_file_name)
  #ID rate
  tmp_summary_distinct <- tmp_data_input %>%
    dplyr::filter(.data$EG.Qvalue.formatted <= parameter_input$ion_q_value_cutoff) %>%
    dplyr::group_by(.data$R.FileName,
                    .data$R.Condition,
                    .data$R.Replicate) %>%
    dplyr::summarise(distinct_ions = dplyr::n_distinct(.data$ion),
                     distinct_modified_peptides = dplyr::n_distinct(.data$EG.ModifiedPeptide),
                     distinct_peptides = dplyr::n_distinct(.data$PEP.StrippedSequence),
                     distinct_proteins = dplyr::n_distinct(.data$PG.ProteinGroups))%>%
    dplyr::ungroup()

  message_function(text = "performing protein count over replicates (more than or equal 2 peptides) ...",
                   color = "blue",
                   log_file_name = log_file_name)

  replicate_count_ID<- dplyr::left_join(tmp_data_input %>%
                                   dplyr::filter(.data$EG.Qvalue.formatted <= parameter_input$ion_q_value_cutoff) %>%
                                   dplyr::group_by(.data$R.FileName,
                                                   .data$R.Condition,
                                                   .data$R.Replicate,
                                                   .data$PG.ProteinGroups) %>%
                                   dplyr::filter(dplyr::n_distinct(.data$PEP.StrippedSequence)>=2) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::group_by(.data$R.Condition,
                                                   .data$PG.ProteinGroups) %>%
                                   dplyr::summarise(replicate_count = dplyr::n_distinct(.data$R.Replicate))%>%
                                   dplyr::ungroup(),
                                tmp_data_input %>%
                                  dplyr::group_by(.data$R.Condition) %>%
                                  dplyr::summarise(replicate_count_total = dplyr::n_distinct(.data$R.Replicate))%>%
                                  dplyr::ungroup(),
                                by= "R.Condition")

  #calculate percentage
  replicate_count_ID <- replicate_count_ID %>%
    dplyr::mutate(replicate_percentage = .data$replicate_count / .data$replicate_count_total*100)

  #generate wide format
  replicate_count_ID_percentage_wide <- replicate_count_ID %>%
    dplyr::select(.data$R.Condition,
                  .data$PG.ProteinGroups,
                  .data$replicate_percentage) %>%
    pivot_wider(names_from = .data$R.Condition,
                values_from = .data$replicate_percentage)
  #replace NA (not there) with 0
  replicate_count_ID_percentage_wide[is.na(replicate_count_ID_percentage_wide)] <- 0

  message_function(text = "performing protein count (<2 and more than or equal 2 peptides) ...",
                   color = "blue",
                   log_file_name = log_file_name)

  #ID rate proteins with 1 or ≥2 peptides
  protein_count_ID <- dplyr::bind_rows( tmp_data_input %>%
                                   dplyr::filter(.data$EG.Qvalue.formatted <= parameter_input$ion_q_value_cutoff) %>%
                                   dplyr::group_by(.data$R.FileName,
                                                   .data$R.Condition,
                                                   .data$R.Replicate,
                                                   .data$PG.ProteinGroups) %>%
                                   dplyr::filter(dplyr::n_distinct(.data$PEP.StrippedSequence)<2) %>%
                                   dplyr::ungroup() %>%
                                   dplyr::group_by(.data$R.FileName,
                                                   .data$R.Condition,
                                                   .data$R.Replicate) %>%
                                   dplyr::summarise(protein_count = dplyr::n_distinct(.data$PG.ProteinGroups))%>%
                                   dplyr::ungroup() %>%
                                   dplyr::mutate(number_of_peptides = "< 2"),
  tmp_data_input %>%
    dplyr::filter(.data$EG.Qvalue.formatted <= parameter_input$ion_q_value_cutoff) %>%
    dplyr::group_by(.data$R.FileName,
                    .data$R.Condition,
                    .data$R.Replicate,
                    .data$PG.ProteinGroups) %>%
    dplyr::filter(dplyr::n_distinct(.data$PEP.StrippedSequence)>=2) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$R.FileName,
                    .data$R.Condition,
                    .data$R.Replicate) %>%
    dplyr::summarise(protein_count = dplyr::n_distinct(.data$PG.ProteinGroups))%>%
    dplyr::ungroup() %>%
    dplyr::mutate(number_of_peptides = ">= 2")
  )


  protein_count_ID_wo_filtering <- dplyr::bind_rows( tmp_data_input %>%
                                               dplyr::group_by(.data$R.FileName,
                                                               .data$R.Condition,
                                                               .data$R.Replicate,
                                                               .data$PG.ProteinGroups) %>%
                                               dplyr::filter(dplyr::n_distinct(.data$PEP.StrippedSequence)<2) %>%
                                               dplyr::ungroup() %>%
                                               dplyr::group_by(.data$R.FileName,
                                                               .data$R.Condition,
                                                               .data$R.Replicate) %>%
                                               dplyr::summarise(protein_count = dplyr::n_distinct(.data$PG.ProteinGroups)) %>%
                                               dplyr::ungroup() %>%
                                               dplyr::mutate(number_of_peptides = "< 2"),
                                tmp_data_input %>%
                                  dplyr::group_by(.data$R.FileName,
                                                  .data$R.Condition,
                                                  .data$R.Replicate,
                                                  .data$PG.ProteinGroups) %>%
                                  dplyr::filter(dplyr::n_distinct(.data$PEP.StrippedSequence)>=2) %>%
                                  dplyr::ungroup() %>%
                                  dplyr::group_by(.data$R.FileName,
                                                  .data$R.Condition,
                                                  .data$R.Replicate) %>%
                                  dplyr::summarise(protein_count = dplyr::n_distinct(.data$PG.ProteinGroups)) %>%
                                  dplyr::ungroup() %>%
                                  dplyr::mutate(number_of_peptides = ">= 2")
  )


  #ion ID median
  ion_id_median <- stats::median(tmp_summary_distinct$distinct_ions)
  #ion ID cutoff
  ion_id_cutoff <- stats::median(tmp_summary_distinct$distinct_ions)-stats::median(tmp_summary_distinct$distinct_ions)*parameter_input$id_drop_cutoff

  #id_outliers
  tmp_summary_distinct_outlier <- tmp_summary_distinct[base::which(tmp_summary_distinct$distinct_ions<ion_id_cutoff),]


  #add outlier indication
  tmp_summary_distinct <- tmp_summary_distinct %>%
    dplyr::mutate(ion_ID_outlier=ifelse(test = .data$R.FileName %in% tmp_summary_distinct_outlier$R.FileName,
                                        yes = "yes",
                                        no = "no"))

  #write ID percentage over replicates ≥2 peptides
  readr::write_csv(x = replicate_count_ID,
              file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count_over_replicates_min_2_pep_Qvalue_cutoff_",parameter_input$ion_q_value_cutoff,".csv"))
  readr::write_csv(x = replicate_count_ID_percentage_wide,
            file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count_over_replicates_min_2_pep_percentage_WIDE_Qvalue_cutoff_",parameter_input$ion_q_value_cutoff,".csv"))

  #write ID counts proteins <2 and ≥2 peptides
  readr::write_csv(x = protein_count_ID %>%
                mutate(Qvalue_below = parameter_input$ion_q_value_cutoff),
              file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count__strippedPEP.csv"))


  #write ID counts proteins <2 and ≥2 peptides
  readr::write_csv(x = protein_count_ID_wo_filtering,
              file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count__strippedPEP_without_Qvalue_cut.csv"))

  #write samples included in analysis
  readr::write_csv(x = sample_data_tmp_data_input,
              file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/file_list.csv"))

  #write summary ID rate output
  readr::write_csv(x = tmp_summary_distinct,
              file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/ID_counts.csv"))

  #write summary ID rate output
  readr::write_csv(x = tmp_summary_distinct_ions,
              file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/ion_ID_counts_fractions.csv"))


  tmp_summary_distinct_gather <- tmp_summary_distinct %>%
    tidyr::pivot_longer(cols = c("distinct_ions",
                          "distinct_modified_peptides",
                          "distinct_peptides",
                          "distinct_proteins"),
                  names_to = "fraction",
                  values_to = "count") %>%
    dplyr::mutate(ion_ID_outlier = ifelse(test = .data$R.FileName %in% tmp_summary_distinct_outlier$R.FileName,
                                   yes = "yes",
                                   no = "no"))


  #do plotting ====
  ion_ID_rate_plot_colors <- c("dodgerblue4","orangered4","darkgrey")
  names(ion_ID_rate_plot_colors) <- c(paste("<",parameter_input$ion_q_value_cutoff,sep=""),"profiled",paste(">",parameter_input$ion_q_value_cutoff,sep=""))

  #protein ID rate plot (1 and 2 peptides)
  max_value_ID <- max(unlist(protein_count_ID_wo_filtering %>%
                               group_by(.data$R.FileName) %>%
                               summarise(sum = sum(.data$protein_count)) %>%
                               ungroup() %>%
                               select(sum)))

  protein_2_peptide_ID_rate_plot <- ggplot2::ggplot(data = protein_count_ID,
                                                    mapping = aes(x = .data$R.FileName,
                                                                  y = .data$protein_count,
                                                                  fill=.data$number_of_peptides))+
    geom_bar(stat = "identity")+
    scale_fill_brewer(palette = "Set1")+
    theme_minimal(base_size = 16)+
    scale_y_continuous(sec.axis = sec_axis(~(./max_value_ID)*100, name = "%"))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=9),
          strip.text.x = element_text(size=16),
          axis.ticks.y = element_line(color="black"),
          panel.grid = element_blank())+
    labs(title="protein ID rate",
         fill="peptide count",
         subtitle = paste("Q-value filter <",parameter_input$ion_q_value_cutoff),
         x="raw file name",
         caption = "100% >> protein IDs whole dataset")


  #protein ID rate plot (1 and 2 peptides) without Q-value filtering


  protein_2_peptide_ID_rate_plot_wo <- ggplot2::ggplot(data = protein_count_ID_wo_filtering,
                                                        mapping = aes(x = .data$R.FileName,
                                                                      y = .data$protein_count,
                                                                      fill = .data$number_of_peptides))+
    geom_bar(stat = "identity")+
    scale_fill_brewer(palette = "Set1")+
    theme_minimal(base_size = 16)+
    scale_y_continuous(sec.axis = sec_axis(~(./max_value_ID)*100, name = "%"))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=9),
          strip.text.x = element_text(size=16),
          axis.ticks.y = element_line(color="black"),
          panel.grid = element_blank())+
    labs(title="protein ID rate",
         fill="peptide count",
         subtitle = paste("without Q-value filter"),
         x="raw file name",
         caption = "100% >> protein IDs whole dataset")



  #ion ID rate plot
  ion_ID_rate_plot <- tmp_summary_distinct_ions %>%
                      pivot_longer(colnames(tmp_summary_distinct_ions)[4:ncol(tmp_summary_distinct_ions)],
                                   names_to = "fraction",
                                   values_to = "count") %>%
    ggplot2::ggplot(mapping = aes(x = .data$R.FileName,
                                  y = .data$count,
                                  fill = forcats::fct_rev(.data$fraction))) +
                      geom_bar(stat="identity")+
                      theme_minimal(base_size = 16)+
                      scale_fill_manual(values = ion_ID_rate_plot_colors)+
                      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=9),
                            strip.text.x = element_text(size=16),
                            axis.ticks.y = element_line(color="black"),
                            panel.grid = element_blank())+
                      labs(title = "ion ID rate",
                           fill = "filtering/\nprofiled",
                           subtitle = paste("Q-value filter <",parameter_input$ion_q_value_cutoff),x="raw file name")+
                      scale_y_continuous(sec.axis = sec_axis(~(./length(unique(tmp_data_input$ion)))*100, name = "% of total"))


  #overall plot
  ID_rate_plot <- ggplot2::ggplot(data = tmp_summary_distinct_gather,
                                 mapping = aes(x = .data$R.FileName,
                                               y = .data$count,
                                               fill = .data$ion_ID_outlier))+
    geom_bar(stat="identity")+
    facet_wrap(~.data$fraction,scales = "free",ncol=1)+
    theme_minimal(base_size = 16)+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=9),
          strip.text.x = element_text(size=16))+
    scale_fill_manual(values = c("yes"="orangered","no"="darkgrey"))+
    labs(title="ID rate",
         subtitle = paste("Q-value <",parameter_input$ion_q_value_cutoff),
         x="raw file name")

  #ion ID rate cutoff plot
  ID_rate_plot_filter <- ggplot2::ggplot(data = tmp_summary_distinct_gather %>%
                                           dplyr::filter(.data$fraction == "distinct_ions"),
                                mapping = aes(x = .data$R.FileName,
                                              y = .data$count,
                                              fill = .data$ion_ID_outlier))+
    geom_bar(stat="identity")+
    theme_minimal(base_size = 16)+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=9))+
    scale_fill_manual(values = c("yes"="orangered","no"="darkgrey"))+
    geom_hline(yintercept = ion_id_cutoff)+
    geom_hline(yintercept = ion_id_median,
               linetype="dotted")+
    labs(title="ion ID rate",
         subtitle = paste("Q-value <",parameter_input$ion_q_value_cutoff),x="raw file name",
         caption=c("dotted line = median / solid line = cutoff"))

  #writeOutputPlots in folder
  number_of_conditions <- length(unique(tmp_summary_distinct_gather$R.Condition))

  ggsave_pdf_png(filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count__strippedPEP_without_Qvalue_cut",sep=""),
         plot = if(number_of_conditions>1){
           protein_2_peptide_ID_rate_plot_wo+
             facet_wrap(~.data$R.Condition, nrow = 1, scales = "free_x")+
             theme(strip.background = element_rect(fill = "grey90", color = NA),
                   strip.text = element_text(color = "grey10"))
         }else{
           protein_2_peptide_ID_rate_plot_wo
         },
         height = 10,
         width = fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions),
         limitsize = FALSE)

  ggsave_pdf_png(filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count__strippedPEP"),
         plot = if(number_of_conditions>1){
           protein_2_peptide_ID_rate_plot+
             facet_wrap(~.data$R.Condition, nrow = 1, scales = "free_x")+
             theme(strip.background = element_rect(fill = "grey90", color = NA),
                   strip.text = element_text(color = "grey10"))
         }else{
           protein_2_peptide_ID_rate_plot
         },
         height = 10,
         width = fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions),
         limitsize = FALSE)

  ggsave_pdf_png(filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/ID_ion_counts_plot"),
         plot = if(number_of_conditions>1){
           ion_ID_rate_plot+
             facet_wrap(~.data$R.Condition, nrow = 1, scales = "free_x")+
             theme(strip.background = element_rect(fill = "grey90", color = NA),
                   strip.text = element_text(color = "grey10"))
         }else{
           ion_ID_rate_plot
         },
         height = 10,
         width = fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions),
         limitsize = FALSE)

  ggsave_pdf_png(filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/ID_counts_plot"),
         plot = if(number_of_conditions>1){
           ID_rate_plot+
             facet_grid(.data$fraction~.data$R.Condition,scales = "free")+
             theme(strip.background = element_rect(fill = "grey90", color = NA),
                   strip.text = element_text(color = "grey10"))
         }else{
           ID_rate_plot
         },
         height = 20,
         # cap to 200 width with high number of samples project
         width = ifelse(fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions)>200,
                        200,
                        fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions)),
         limitsize = FALSE)

  ggsave_pdf_png(filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/ID_counts_plot_ion_filter"),
         plot = if(number_of_conditions>1){
           ID_rate_plot_filter+
             facet_wrap(~.data$R.Condition, nrow = 1, scales = "free_x")+
             theme(strip.background = element_rect(fill = "grey90", color = NA),
                   strip.text = element_text(color = "grey10"))
         }else{
           ID_rate_plot_filter
         },
         height = 10,
         width = fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions),
         limitsize = FALSE)

  #plot if print out == T
  if(print.plot==T){
    print(ID_rate_plot_filter)
  }



# add ON/OFF analysis per condition -----------------------------------------------------

  message_function(text = "_________ ON/OFF analysis: _________", color = "blue",log_file_name = log_file_name)
  message_function(text = paste0("... filter with Q-value ",parameter_input$ion_q_value_cutoff," ..."), color = "blue",log_file_name = log_file_name)

  condition_replicate_count<- tmp_data_input %>%
    group_by(.data$R.Condition) %>%
    summarise(replicate_total_count = n_distinct(.data$R.Replicate)) %>%
    ungroup()

  #PG with 2 peptides
  PG_2_peptides_ID <- tmp_data_input %>%
    dplyr::filter(.data$EG.Qvalue.formatted <= parameter_input$ion_q_value_cutoff) %>%
    dplyr::group_by(.data$R.FileName,
             .data$R.Condition,
             .data$PG.ProteinGroups) %>%
    dplyr::filter(n_distinct(.data$PEP.StrippedSequence)>=2) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$R.Condition,
             .data$PG.ProteinGroups) %>%
    dplyr::summarise(peptide_count = n_distinct(.data$PEP.StrippedSequence),
              present_replicate_count = n_distinct(.data$R.Replicate))
  PG_2_peptides_ID <- dplyr::left_join(PG_2_peptides_ID,condition_replicate_count,by = "R.Condition") %>%
    dplyr::mutate(present_replicate_percentage = .data$present_replicate_count / .data$replicate_total_count * 100)

  message_function(text = paste0("... filter for 2 peptides and min. present in 50% of replicates ..."), color = "blue",log_file_name = log_file_name)

  #placeholder for raw data before 50% filtering (must be present in 50% of replicates)
  PG_2_peptides_ID_raw <- PG_2_peptides_ID

  #50% filtering over replicates
  PG_2_peptides_ID <- PG_2_peptides_ID %>%
    dplyr::filter(.data$present_replicate_percentage>=50)


  PG_2_peptides_ID_wide<- PG_2_peptides_ID %>%
    dplyr::select(-.data$present_replicate_count,
                  -.data$replicate_total_count,
                  -.data$present_replicate_percentage) %>%
    pivot_wider(names_from = "R.Condition",values_from = "peptide_count")

  #replace NA with 0 >> 0 peptides found
  PG_2_peptides_ID_wide[is.na(PG_2_peptides_ID_wide)] <- 0

  #binary coded table
  PG_2_peptides_ID_wide_binary<- PG_2_peptides_ID %>%
    dplyr::mutate(binary_peptide_count = if_else(condition = .data$peptide_count>0,
                                                 true = 1,
                                                 false = 0)) %>%
    dplyr::select(-.data$peptide_count,
                  -.data$present_replicate_count,
                  -.data$replicate_total_count,
                  -.data$present_replicate_percentage) %>%
    pivot_wider(names_from = "R.Condition",
                values_from = "binary_peptide_count")

  #NA to zero
  PG_2_peptides_ID_wide_binary[is.na(PG_2_peptides_ID_wide_binary)] <- 0

  condition_length <- length(unique(PG_2_peptides_ID$R.Condition))

  message_function(text = "... ON/OFF analysis write outputs ...", color = "blue",log_file_name = log_file_name)

  #save upsetR plot
  png(filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/Detected_ProteinGroups__UpSetR__plot.png"),
      width = 6+condition_length,height = 8,res = 600,units="in")
  if(condition_length==1){
    ggplot()+
      theme_void()+
      labs(title = "no Upset plot possible >> only one condition")

  }else{
    #plot upset plot
    print(UpSetR::upset(as.data.frame(PG_2_peptides_ID_wide_binary),
                        nintersects = NA,
                        nsets = condition_length,
                        order.by = c("freq")))
  }

  dev.off()

  #write PG_2_peptides_ID_wide_binary
  readr::write_csv(x = PG_2_peptides_ID_wide_binary,
            file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/Detected_ProteinGroups__UpSetR__plot__binary_coded.csv"))

  #write PG_2_peptides_ID_wide
  readr::write_csv(x = PG_2_peptides_ID_wide,
            file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/Detected_ProteinGroups__UpSetR__plot__stripped_peptide_count_per_condition.csv"))

  #write excel file
  #create workbook
  excel_output <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = excel_output, sheetName = "PG_with_2_peptides_binary")
  openxlsx::addWorksheet(wb = excel_output, sheetName = "PG_with_2_peptides_pep_count")

  #write data into workbook
  openxlsx::writeDataTable(wb = excel_output,sheet =  1,
                           x = PG_2_peptides_ID_wide_binary,
                           startRow = 1,
                           startCol = 1,
                           tableStyle = "TableStyleMedium1")
  openxlsx::writeDataTable(wb = excel_output,
                           sheet =  2,
                           x = PG_2_peptides_ID_wide,
                           startRow = 1,
                           startCol = 1,
                           tableStyle = "TableStyleMedium1")

  #conditional formating
  openxlsx::conditionalFormatting(wb = excel_output,
                                  sheet = 1,
                                  cols = 2:ncol(PG_2_peptides_ID_wide_binary),
                                  rows = 1:nrow(PG_2_peptides_ID_wide_binary)+1,
                                  style = c("#B71C1C", "#1B5E20"),
                                  rule = c(1, 0),
                                  type = "colourScale")

  openxlsx::conditionalFormatting(wb = excel_output,
                                  sheet = 2,
                                  cols = 2:ncol(PG_2_peptides_ID_wide),
                                  rows = 1:nrow(PG_2_peptides_ID_wide)+1,
                                  style = c("#B71C1C", "#1B5E20"),
                                  rule = c(1, 0),
                                  type = "colourScale")

  openxlsx::saveWorkbook(wb = excel_output,
                         overwrite = T,
                         file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/Detected_ProteinGroups__stripped_peptide_count_per_condition_raw_data.xlsx"))


  #write PG_2_peptides_ID_raw_data
  readr::write_csv(x = PG_2_peptides_ID_raw,
            file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/Detected_ProteinGroups__stripped_peptide_count_per_condition_raw_data.csv"))

  # adding MC plots for ion data --------------------------------------------
  message_function(text = "count missed cleavages ...",
                   color = "blue",
                   log_file_name = log_file_name)
  #calc. missed cleavage count per sample with q-value cut
  mc_count_data <- tmp_data_input %>%
    dplyr::filter(.data$EG.Qvalue.formatted<parameter_input$ion_q_value_cutoff) %>%
    dplyr::distinct(.data$R.FileName,
                    .data$R.Condition,
                    .data$PEP.StrippedSequence,
                    .data$PEP.NrOfMissedCleavages) %>%
    dplyr::group_by(.data$R.FileName,
                    .data$R.Condition,
                    .data$PEP.NrOfMissedCleavages) %>%
    dplyr::summarise(MC_count = n_distinct(.data$PEP.StrippedSequence)) %>%
    ungroup()

  # wide table missed cleavages
  mc_count_data_wide <- mc_count_data %>%
    pivot_wider(names_from = "PEP.NrOfMissedCleavages",
                values_from = "MC_count",
                names_prefix = "MC_")

  # vector for MC
  mc_vec <- colnames(mc_count_data_wide)[-c(1,2)]
  # adding MC percentage
  mc_count_data_wide <- mc_count_data_wide %>%
    group_by(.data$R.FileName,.data$R.Condition) %>%
    mutate(total_peptide_count = sum(c_across(mc_vec[1:length(mc_vec)]),na.rm=T),
           MC_peptide_count = sum(c_across(mc_vec[2:length(mc_vec)]),na.rm=T)) %>%
    mutate(MC_percentage = .data$MC_peptide_count/.data$total_peptide_count*100)

  # write missed cleavage out data
  write_csv(x = mc_count_data_wide,
            file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/missed_cleavages_sample_wise_PERCENTAGE.csv"))

  # missed cleavage count of whole project
  mc_count_data_global <- tmp_data_input %>%
    dplyr::filter(.data$EG.Qvalue.formatted<parameter_input$ion_q_value_cutoff) %>%
    dplyr::distinct(.data$PEP.StrippedSequence,
                    .data$PEP.NrOfMissedCleavages) %>%
    dplyr::group_by(.data$PEP.NrOfMissedCleavages) %>%
    dplyr::summarise(MC_count = n_distinct(.data$PEP.StrippedSequence))

  # total number of peptides in project
  n_total_peptides <- sum(mc_count_data_global$MC_count)

  message_function(text = "save missed cleavage plots...",
                   color = "blue",
                   log_file_name = log_file_name)
  # missed cleavage plot
  global_missed_cleavage_plot <- ggplot(data = mc_count_data_global,
                                        mapping = aes(x = .data$PEP.NrOfMissedCleavages,
                                                      y = .data$MC_count,
                                                      fill = as.factor(.data$PEP.NrOfMissedCleavages)))+
    geom_bar(stat = "identity")+
    scale_fill_brewer(palette = "Set1")+
    theme_light(base_size = 18)+
    scale_y_continuous(sec.axis = sec_axis(~(./n_total_peptides)*100, name = "%"))+
    labs(title =" missed cleavages",
         y = "peptide count",
         x = "missed cleavage",
         fill = "MC",
         caption = "percentage to total peptide count of project")

  # save global missed cleavage plot
  ggsave_pdf_png(plot = global_missed_cleavage_plot,
                 filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/missed_cleavages_global"),
                 width = 8,
                 height = 6)

  missed_cleavage_plot<- ggplot(data = mc_count_data,
                                mapping = aes(x = .data$R.FileName,
                                              y = .data$MC_count,
                                              fill = as.factor(.data$PEP.NrOfMissedCleavages)))+
    geom_bar(stat = "identity", position = "dodge")+
    scale_fill_brewer(palette = "Set1")+
    theme_light(base_size = 16)+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.1))+
    scale_y_continuous(sec.axis = sec_axis(~(./n_total_peptides)*100, name = "%"))+
    labs(title ="count of missed cleavages per sample",
         subtitle = paste("Q-value <",parameter_input$ion_q_value_cutoff),
         y = "peptide count",
         x = "raw file name",
         fill = "MC",
         caption = "percentage to total peptide count of project")

  # save sample wise missed cleavage plot
  ggsave_pdf_png(plot = if(number_of_conditions>1){
                        missed_cleavage_plot+
                          facet_wrap(~.data$R.Condition, nrow = 1, scales = "free_x")+
                          theme(strip.background = element_rect(fill = "grey90", color = NA),
                                strip.text = element_text(color = "grey10"))
                      }else{
                        missed_cleavage_plot
                      },
                 filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/missed_cleavages_sample_wise"),
                 # cap to 200 width with high number of samples project
                 width = fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions),
                 height = 10)

  # missed cleavage plot percentage
  missed_cleavage_plot_percentage<- ggplot(data = mc_count_data_wide,
                                           mapping = aes(x = .data$R.FileName,
                                                         y = .data$MC_percentage))+
    geom_bar(stat = "identity", position = "dodge")+
    scale_fill_brewer(palette = "Set1")+
    theme_light(base_size = 16)+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.1))+
    labs(title ="percentage of missed\ncleavages per sample",
         subtitle = paste("Q-value <",parameter_input$ion_q_value_cutoff),
         y = "MC [%]", x = "raw file name", caption = "")+
    geom_text(mapping = aes(label = paste(round(.data$MC_percentage,digits = 1),"%",sep = "")),
              angle = 90,
              hjust=-0.1,
              color = "black", size = 3)+
    coord_cartesian(ylim = c(0,max(mc_count_data_wide$MC_percentage,na.rm=T)*1.5))

  # save sample wise missed cleavage plot
  ggsave_pdf_png(plot = if(number_of_conditions>1){
                        missed_cleavage_plot_percentage+
                          facet_wrap(~.data$R.Condition, nrow = 1, scales = "free_x")+
                          theme(strip.background = element_rect(fill = "grey90", color = NA),
                                strip.text = element_text(color = "grey10"))
                      }else{
                        missed_cleavage_plot_percentage
                      },
                 filename = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/missed_cleavages_sample_wise_PERCENTAGE"),
                 width = fig_width_estimation(sample_length = sample_length,condition_length = number_of_conditions),
                 height = 10)


  if(dim(tmp_summary_distinct_outlier)[1]==0){
    message_function(text = "_________ ID rate per sample: _________", color = "blue", log_file_name = log_file_name)
    message_function(text = paste("ions: ","median = ",ion_id_median,"; min. = ",min(tmp_summary_distinct$distinct_ions),"; max. = ",max(tmp_summary_distinct$distinct_ions),sep=""), color = "blue", log_file_name = log_file_name)
    message_function(text = paste("modified peptides: ","median = ",median(tmp_summary_distinct$distinct_modified_peptides),"; min. = ",min(tmp_summary_distinct$distinct_modified_peptides),"; max. = ",max(tmp_summary_distinct$distinct_modified_peptides),sep=""), color = "blue", log_file_name = log_file_name)
    message_function(text = paste("stripped peptides: ","median = ",median(tmp_summary_distinct$distinct_peptides),"; min. = ",min(tmp_summary_distinct$distinct_peptides),"; max. = ",max(tmp_summary_distinct$distinct_peptides),sep=""), color = "blue", log_file_name = log_file_name)
    message_function(text = paste("protein groups: ","median = ",median(tmp_summary_distinct$distinct_proteins),"; min. = ",min(tmp_summary_distinct$distinct_proteins),"; max. = ",max(tmp_summary_distinct$distinct_proteins),sep=""), color = "blue", log_file_name = log_file_name)

    output_list <- list(spectronaut_output = tmp_data_input,
                        SDRF_file = SDRF_file,
                        summary_distinct = tmp_summary_distinct,
                        raw_file_names = sample_data_tmp_data_input,
                        ion_id_median = ion_id_median,
                        ion_id_cutoff = ion_id_cutoff,
                        PG_2_peptides_ID_raw = PG_2_peptides_ID_raw,
                        summary_distinct_outlier = tmp_summary_distinct_outlier,
                        ID_rate_plot=ID_rate_plot,
                        ID_rate_plot_filter=ID_rate_plot_filter,
                        sample_length = sample_length,
                        parameter = parameter_input,
                        time_stamp_log_file = time_stamp_log_file,
                        log_file_name = log_file_name)

    message_function(text = paste0("read module done --> please check outputs in folder: ", out_folder,"/","02_ID_rate/"),
                     color = "blue",
                     log_file_name = log_file_name)

    class(output_list) <- "SpectroPipeR_data"

    return(output_list)
  }else{
    message_function(text = "_________ ID rate per sample: _________",color = "blue",log_file_name = log_file_name)
    message_function(text = paste("ions: ","median=",ion_id_median,"; min. = ",min(tmp_summary_distinct$distinct_ions),"; max. = ",max(tmp_summary_distinct$distinct_ions),sep=""),color = "blue",log_file_name = log_file_name)
    message_function(text = paste("modified peptides: ","median = ",median(tmp_summary_distinct$distinct_modified_peptides),"; min. = ",min(tmp_summary_distinct$distinct_modified_peptides),"; max. = ",max(tmp_summary_distinct$distinct_modified_peptides),sep=""),color = "blue",log_file_name = log_file_name)
    message_function(text = paste("stripped peptides: ","median = ",median(tmp_summary_distinct$distinct_peptides),"; min. = ",min(tmp_summary_distinct$distinct_peptides),"; max. = ",max(tmp_summary_distinct$distinct_peptides),sep=""),color = "blue",log_file_name = log_file_name)
    message_function(text = paste("protein groups: ","median=",median(tmp_summary_distinct$distinct_proteins),"; min. = ",min(tmp_summary_distinct$distinct_proteins),"; max. = ",max(tmp_summary_distinct$distinct_proteins),sep=""),color = "blue",log_file_name = log_file_name)
    message_function(text = paste("Attention ---> OUTLIER detected !!! >",parameter_input$id_drop_cutoff*100,"% lower than median ion ID rate"),color = "yellow",log_file_name = log_file_name)
    message_function(text = paste("samples detected with ID rates below cutoff of ",ion_id_cutoff,":"),color = "blue",log_file_name = log_file_name)

    message_function(text = paste(tmp_summary_distinct_outlier$R.FileName,collapse = "\n"),color = "yellow",log_file_name = log_file_name)
    #writing outlier ID rate
    write_csv(x = tmp_summary_distinct_outlier,
                file = paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/ID_counts__detected_outliers.csv"))

    output_list <- list(spectronaut_output = tmp_data_input,
                        SDRF_file = SDRF_file,
                        summary_distinct = tmp_summary_distinct,
                        raw_file_names = sample_data_tmp_data_input,
                        ion_id_median = ion_id_median,
                        ion_id_cutoff = ion_id_cutoff,
                        PG_2_peptides_ID_raw = PG_2_peptides_ID_raw,
                        summary_distinct_outlier = tmp_summary_distinct_outlier,
                        ID_rate_plot = ID_rate_plot,
                        ID_rate_plot_filter = ID_rate_plot_filter,
                        sample_length = sample_length,
                        parameter = parameter_input,
                        time_stamp_log_file = time_stamp_log_file,
                        log_file_name = log_file_name)

    message_function(text = paste0("read module done --> please check outputs in folder: ", out_folder,"/","02_ID_rate/"),
                     color = "blue",
                     log_file_name = log_file_name)

    class(output_list) <- "SpectroPipeR_data"
    return(output_list)

  }

}

