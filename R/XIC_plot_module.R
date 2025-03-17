#' Extraction and plotting of XICs from Spectronaut's exported SQLite Databases (Spectronaut version >  18.7.240506.55695;  XIC DB export generates one SQLite db file per raw file)
#' @description
#' Exports user-selected protein group XIC (MS1/MS2) data in both PDF and CSV formats for all detected ions.
#'
#' @param Spectronaut_report_path file; Path pointing to Spectronaut analysis report using SpectroPipeR report export scheme.
#' @param Spectronaut_xicDB_path folder; Path pointing to Spectronaut all XIC SQLite database export folder from that analysis (use pipeline mode in Spectronaut).
#' @param protein_groups Vector containing user-selected protein group IDs for generating plots. e.g. c("P38720","P29311")
#' @param output_path Define the output directory for the generated ion-specific plots (PDF) and data tables (CSV) for user-selected protein groups.
#' @param export_csv_files Export individual CSV files per ion (boolean, default: FALSE).
#' @param selected_conditions Vector of user-selected conditions for plotting (serves for both selection and ordering).
#' @param number_of_cores Define the number of CPU cores to utilize for parallel XIC extraction and plotting.
#' @param run_specific_y_axis if the XIC plot should be globally scaled y-axis or run specific (default: FALSE)
#' @param Spectronaut_colors TRUE or FALSE if TRUE Spectronaut colors are selected, FALSE (default) SpectroPipeR colors are selected
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom scales squish
#' @import patchwork
#' @import base64enc
#' @import tidyverse
#' @import RSQLite
#' @md
#' @return Exports user-selected protein group XIC data (MS1/MS2) in both PDF and CSV formats for all detected ions of a protein.
#' @details
#' Description of the Spectronaut scores shown in the output plots:
#' | <u> __parameter__ </u> | <u> __description__ </u>|
#' |:---------|:--------------------------------------|
#' | EG.Qvalue  | The Q-value is the multiple test corrected p-value. The Q-value can be used|
#' |            | as an estimate for the FDR at a given Cscore cutoff.|
#' | EG.PEP     | The posterior error probability for a peptide precursor observation. The PEP|
#' |            | is calculated as the decoy probability function for a specficic Cscore divided|
#' |            | by the probability function for decoy plus target. PEP = d/(t+d) where d is the |
#' |            | value of the decoy probability function for a given Cscore and t is the value for|
#' |            | target probability function for a given Cscore. |
#' | FG.ShapeQualityScore | The average of all available peak quality measurements |
#' |                      |(scale: 1.0 `best` to -1.0 `worst`) |
#' | EG.IntCorrScore | The correlation score between the relative fragment intensities as provided|
#' |                 | by the spectral librarxy and the empirical observed fragtment intensitiers |
#' |                 |at XIC peak apex. |
#' | EG.SignalToNoise | The signal to noise ratio of the LC-MS peak for a given peptide precursor.|
#' |                  | The signal is calc. as the maximum intensity of the fragment sum XIC within|
#' |                  | the peak boundaries. The noise is calculated as the average fragment sum XIC|
#' |                  | intensity outside the peak boundaries. |
#' | EG.Cscore  | The discriminant score used for identification. This score is a linear combination|
#' |            | of all applicable scores for a certain workflow. |
#' | EG.Identified | Returns weather a precursor is considered identified or not. In order to be |
#' |               | identified it has to match both the precursor and the protein FDR cutoff. |
#' | dppp (MS1) | Returns the number of data points (scans) that make up the MS1 XIC of the selected|
#' |            |  LC-MS peak for this precursor. Only data points that lie within the start and the|
#' |            | end RT are considered. |
#' | dppp (MS2) | Returns the number of data points (scans) that make up the MS2 XIC of the selected|
#' |            | LC-MS peak for this precursor. Only data points that lie within the start and the|
#' |            | end RT are considered. |
#'
#' The vertical dashed line in the plots indicate the individual ApexRT, which is calculated based on all fragment XIC for a specific peptide precursor ion.
#'
#' @export
#'
#' @examples
#' \donttest{
#' #setup example input paths / protein selection for plotting
#' Spectronaut_report_path <- system.file("extdata/HYE_demo_data",
#'                                        "HYE_demo_data_Report_SpectroPipeR.tsv",
#'                                         package="SpectroPipeR")
#' Spectronaut_xicDB_path <- system.file("extdata/HYE_demo_data/XIC_DBs",package="SpectroPipeR")
#' protein_groups <- c("P29311","P38720")
#' output_path <- "../SpectroPipeR_test_folder/single_XIC_plots"
#'
#' # extracting and plotting of XIC
#' XIC_plot_module(Spectronaut_report_path = Spectronaut_report_path,
#'                 Spectronaut_xicDB_path = Spectronaut_xicDB_path,
#'                 protein_groups = protein_groups,
#'                 output_path = output_path,
#'                 number_of_cores = 2
#' )
#'}

XIC_plot_module <- function(Spectronaut_report_path = NULL,
                            Spectronaut_xicDB_path = NULL,
                            protein_groups = NULL,
                            output_path = NULL,
                            export_csv_files = FALSE,
                            run_specific_y_axis = FALSE,
                            selected_conditions = NULL,
                            Spectronaut_colors = FALSE,
                            number_of_cores = 2){

  #scientific format function of axis
  scientific_10 = function(x) {
    ifelse(
      x==0, "0",
      parse(text = sub("e[+]?", " %*% 10^", scales::scientific(x = x,digits = 1)))
    )
  }

  # XIC colors
  MS1_color<- colorRampPalette(c("#B62F1A","#244B52","#586634","#4499AD","#B58542","#666666"))
  MS2_color<- colorRampPalette(c("#5B478B","#356792","#3D8458","#59A4A4","#893B6C","#80AD56","#E3AF3D","#D4812F","#BD5845","#666666"))


# Spectronaut colors ------------------------------------------------------

  Spectronaut_color_fn <- function(n) {
    # Parameters
    step <- 300.0 / n
    s <- 1
    b <- 0.97

    # Generate colors
    colors <- sapply(0:(n-1), function(i) {
      h <- (i * step) %% 360 # Ensure the hue stays within 0-360 degrees
      hsv(h/360, s, b) # Convert HSV to an RGB color string
    })

    return(colors)
  }

# create output dir ------------------------------------------------------
  if(dir.exists(output_path)){
    message("The output directory already exists, and there is a possibility that files within it may be overwritten.")

  }else{
    message("create output dir ...")
    dir.create(output_path,recursive = T,showWarnings = F)
  }

# load SN report data -----------------------------------------------------
  message("loading Spectronaut report data ...")
  tmp_data_input <- suppressWarnings(readr::read_delim(file = Spectronaut_report_path,
                                        delim = "\t",
                                        col_names = T,
                                        show_col_types = FALSE,
                                        col_types = readr::cols(EG.Qvalue = readr::col_character(),
                                                                PG.IBAQ = readr::col_character())))

# check colnames of SN report ---------------------------------------------

  #input cols for check
  data_input_cols <- c("R.Condition",
                       "R.Replicate",
                       "R.FileName",
                       "PG.ProteinGroups",
                       "FG.IntMID",
                       "FG.Id",
                       "FG.XICDBID",
                       "EG.Qvalue",
                       "EG.PEP",
                       "EG.PeakWidth",
                       "EG.DatapointsPerPeak",
                       "EG.DatapointsPerPeak (MS1)",
                       "EG.Cscore",
                       "EG.SignalToNoise",
                       "EG.Identified",
                       "EG.ApexRT",
                       "EG.IntCorrScore",
                       "FG.ShapeQualityScore")

  # test if all needed columns are in data_input
  if(sum(colnames(tmp_data_input)%in%data_input_cols)!=length(data_input_cols)){
    stop_text <- paste(paste("essentiell columns are missing in input data:",
                             paste(data_input_cols[data_input_cols%nin%colnames(tmp_data_input)],
                                   collapse = "; ")))
    message(stop_text)
    stop()
  }

# get ion meta data ------------------------------------------------------
  message("extract ion meta data from report ...")
  meta_data_ions <- tmp_data_input %>%
    distinct(.data$R.Condition,
             .data$R.Replicate,
             .data$R.FileName,
             .data$PG.ProteinGroups,
             .data$EG.Qvalue,
             .data$EG.PEP,
             .data$EG.DatapointsPerPeak,
             .data$`EG.DatapointsPerPeak (MS1)`,
             .data$EG.Cscore,
             .data$EG.SignalToNoise,
             .data$EG.Identified,
             .data$EG.ApexRT,
             .data$EG.IntCorrScore,
             .data$FG.ShapeQualityScore,
             .data$FG.Id,
             .data$FG.XICDBID) %>%
    dplyr::mutate(EG.Qvalue_numeric = suppressWarnings(as.numeric(.data$EG.Qvalue))) %>%
    dplyr::mutate(EG.Qvalue_formatted = scales::scientific(x = .data$EG.Qvalue_numeric,
                                                           digits = 3)) %>%
    dplyr::mutate(EG.Qvalue_formatted = ifelse(is.na(.data$EG.Qvalue_numeric),
                                               .data$EG.Qvalue,
                                               .data$EG.Qvalue_formatted)
    ) %>%
    dplyr::mutate(EG.PEP_numeric = suppressWarnings(as.numeric(.data$EG.PEP))) %>%
    dplyr::mutate(EG.PEP_formatted = scales::scientific(x = .data$EG.PEP_numeric,
                                                           digits = 3)) %>%
    dplyr::mutate(EG.PEP_formatted = ifelse(is.na(.data$EG.PEP_numeric),
                                               .data$EG.PEP,
                                               .data$EG.PEP_formatted)
    )


  #check selected_conditions
  if(!is.null(selected_conditions)){
    if(sum(unique(meta_data_ions$R.Condition) %in%
           selected_conditions) != length(selected_conditions)){
      stop("The conditions provided by via selected_conditions parameter does not overlap completely the Spectronaut report R.Condition entries.")
    }
  }

# Precursor DB ID of selection --------------------------------------------
  message("extract ion IDs of selected proteins groups ...")
  fgid_selection <- unlist(meta_data_ions %>%
    filter(.data$PG.ProteinGroups %in% protein_groups) %>%
    distinct(.data$FG.XICDBID)
    )

  fgid_selection_PG <- meta_data_ions %>%
    filter(.data$PG.ProteinGroups %in% protein_groups) %>%
    distinct(.data$PG.ProteinGroups)

# Precursor DB ID of selection error check --------------------------------
  if(nrow(fgid_selection_PG)!=length(protein_groups)){
    if(length(protein_groups[protein_groups %nin% meta_data_ions$PG.ProteinGroups])==1){
      stop(paste0("The selected protein group ",
                  paste(protein_groups[protein_groups %nin% meta_data_ions$PG.ProteinGroups],collapse = ","),
                  " is not present in your analysis report. Use only valid protein group entries!")
      )
    }else{
      stop(paste0("The selected protein groups ",
                  paste(protein_groups[protein_groups %nin% meta_data_ions$PG.ProteinGroups],collapse = ","),
                  " are not present in your analysis report. Use only valid protein group entries!")
      )
    }
  }

# extact ions from database -----------------------------------------------
  # single xicDB file Spectronaut version < 18.7.240506.55695
  #
  #message("connect Spectronaut XIC database ...")
  #db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = Spectronaut_xicDB_path)

  #message("extracting ion XICs from the database...")
  #result <- as_tibble(RSQLite::dbGetQuery(conn = db,
  #                              statement = "SELECT * FROM IonTraces INNER JOIN RTAxis ON Iontraces.RTAxis_ID=RTAxis.ID WHERE FGID = ?",
  #                              params = list(fgid_selection)))
  # extract runs from database
  #result_runTable <- as_tibble(RSQLite::dbGetQuery(conn = db,
  #                                                statement = "SELECT * FROM Run"))


  # one SQLite file per raw file extraction Spectronaut version > 18.7.240506.55695
  xicDB_files <- list.files(path = Spectronaut_xicDB_path,pattern = "xic.db",full.names = T)

  #create output tables
  result <- tibble()
  result_runTable <- tibble()

  # multicore extraction of XIC
  message("extracting XIC (this might take a while) ...")
  cl <- parallel::makeCluster(number_of_cores, setup_strategy = "sequential") # number of cores.
  doSNOW::registerDoSNOW(cl) #register cores
  #setup the progressbar
  progress <- function(n) setTxtProgressBar(pb = txtProgressBar(max = length(xicDB_files),
                                                                width=60L,
                                                                style = 3),n)
  opts <- list(progress = progress)
  # multicore plotting
  result  <- foreach::foreach(i = 1:length(xicDB_files),
                   .combine = rbind,
                   .packages = c("tidyverse","RSQLite"),
                   .options.snow = opts,
                   .multicombine = T) %dopar% {
                     progress(n = i) #update progressbar

                     #connect SQLite db
                     db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = xicDB_files[i])
                     #extract results
                     result_tmp <- as_tibble(RSQLite::dbGetQuery(conn = db,
                                                                 statement = "SELECT * FROM IonTraces INNER JOIN RTAxis ON Iontraces.RTAxis_ID=RTAxis.ID WHERE FGID = ?",
                                                                 params = list(fgid_selection)))
                     dbDisconnect(db)
                     result_tmp

                   }

  parallel::stopCluster(cl)#closing cluster for parallel computing
  message("")

  message("extracting XIC runs (this might take a while) ...")
  cl <- parallel::makeCluster(number_of_cores, setup_strategy = "sequential") # number of cores.
  doSNOW::registerDoSNOW(cl) #register cores
  #setup the progressbar
  progress <- function(n) setTxtProgressBar(pb = txtProgressBar(max = length(xicDB_files),
                                                              width=60L,
                                                              style = 3),n)
  opts <- list(progress = progress)

  # multicore extraction of runs
  result_runTable  <- foreach::foreach(i = 1:length(xicDB_files),
                            .combine = rbind,
                            .packages = c("tidyverse","RSQLite"),
                            .options.snow = opts,
                            .multicombine = T) %dopar% {
                              progress(n = i) #update progressbar

                              #connect SQLite db
                              db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = xicDB_files[i])
                              #extract results
                              result_runTable_tmp <- as_tibble(RSQLite::dbGetQuery(conn = db,
                                                                                   statement = "SELECT * FROM Run"))
                              dbDisconnect(db)
                              result_runTable_tmp

                            }

  parallel::stopCluster(cl)#closing cluster for parallel computing
  message("")

  # extraction error check
  if(dim(result)[1]==0){
    stop(paste0("No data could be extracted from ",Spectronaut_xicDB_path,". Have you supplied a legitimate SQLite database that corresponds to the same analysis as the one detailed in the submitted Spectronaut report?"))
  }




  #nest extraction XIC
  message("reformat extracted XIC(s) ...")
  result_data_nested <- result %>%
    dplyr::rowwise() %>%
    #nest extraction XIC
    dplyr::mutate(RTAxis_Length = ifelse(test = .data$RTAxis_XOffset==0,
                                         yes = .data$RTAxis_Length+1,
                                         no = .data$RTAxis_Length)) %>%
    dplyr::mutate(decoded_data_XIC = list(base64enc::base64decode(.data$IntensityDataBytes))) %>%
    dplyr::mutate(num_floats_XIC = length(.data$decoded_data_XIC) / 4) %>%
    dplyr::mutate(XIC = list(readBin(.data$decoded_data_XIC,
                                     what = "numeric",
                                     n = .data$num_floats_XIC,
                                     size = 4,
                                     endian = "little"))) %>%
    #nest extraction RT
    dplyr::mutate(decoded_data_rt = list(base64enc::base64decode(.data$RTDataBytes))) %>%
    dplyr::mutate(num_floats_rt = length(.data$decoded_data_rt) / 4) %>%
    dplyr::mutate(RT = list(readBin(.data$decoded_data_rt,
                                    what = "numeric",
                                    n = .data$num_floats_rt,
                                    size = 4,
                                    endian = "little"))) %>%
      # adjust RT offset
    mutate(RT = list(.data$RT[c(.data$RTAxis_XOffset):(c(.data$RTAxis_XOffset + .data$RTAxis_Length)-1)]))






  # unnest data
  merged_extracted_data_unnested<- result_data_nested %>%
    dplyr::select(-.data$IntensityDataBytes,
                  -.data$num_floats_XIC,
                  -.data$decoded_data_XIC,
                  -.data$RTDataBytes,
                  -.data$decoded_data_rt,
                  -.data$num_floats_rt,
                  -.data$RTAxis_XOffset,
                  -.data$RTAxis_Length,
                  -.data$RTAxis_ID,
                  -.data$ID
                  ) %>%
    tidyr::unnest(cols = c(.data$XIC, .data$RT))


# add ion meta data to extracted XIC --------------------------------------
  message("add ion meta data to extracted XICs ...")
  # add raw file name
  merged_extracted_data_unnested_tmp1 <- dplyr::left_join(merged_extracted_data_unnested,
                                                           result_runTable %>%
                                                            rowwise() %>%
                                                             mutate(RawFileName = unlist(stringr::str_split(.data$RawFileName,pattern = "\\.",n = 2))[1]) %>%
                                                            ungroup(),
                                                           by = c("Run_ID" = "ID"))


  # add meta ion data
  merged_extracted_data_unnested_tmp2 <- left_join(merged_extracted_data_unnested_tmp1 %>%
                                                      mutate(FGID = as.numeric(.data$FGID)),
                                                    meta_data_ions %>%
                                                      distinct(.data$PG.ProteinGroups,
                                                               .data$FG.Id,
                                                               .data$FG.XICDBID),
                                                    by = c("FGID" =  "FG.XICDBID"))

  # add meta run data
  merged_extracted_data_unnested_tmp3 <- left_join(merged_extracted_data_unnested_tmp2 %>%
                                                      mutate(FGID = as.numeric(.data$FGID)),
                                                    meta_data_ions %>%
                                                      dplyr::filter(.data$FG.XICDBID%in%fgid_selection) %>%
                                                      dplyr::distinct(.data$R.Condition,
                                                                      .data$R.Replicate,
                                                                      .data$R.FileName,
                                                                      .data$FG.Id),                                                     ,
                                                    by = c("RawFileName" =  "R.FileName",
                                                           "FG.Id")
                                                    )




  merged_extracted_data_unnested_final <- merged_extracted_data_unnested_tmp3

  # filter meta ion data for conditions -------------------------------------
  if(!is.null(selected_conditions)){
    message("filter ion meta data for user selected conditions...")
    meta_data_ions <- meta_data_ions %>%
      dplyr::filter(.data$R.Condition %in% selected_conditions)
  }

  # filter meta ion data for conditions -------------------------------------
  if(!is.null(selected_conditions)){
    message("filter XIC user selected conditions...")
    merged_extracted_data_unnested_final <- merged_extracted_data_unnested_final %>%
      dplyr::filter(.data$R.Condition %in% selected_conditions)
  }

  #extract total sample length
  sample_length <- length(unique(meta_data_ions$R.FileName))

# plot ion XIC ------------------------------------------------------------

  #get ions
  ions <- unique(merged_extracted_data_unnested_final$FG.Id)

    message("register processor cores")
    cl <- parallel::makeCluster(number_of_cores, setup_strategy = "sequential") # number of cores.
    doSNOW::registerDoSNOW(cl) #register cores

  #setup the progressbar
    progress <- function(n) setTxtProgressBar(pb = txtProgressBar(max = length(ions),width=60L, style = 3), n)
    opts <- list(progress = progress)

message("plotting (this might take a while) ...")
# multicore plotting
foreach::foreach(i = 1:length(ions),
                .combine=rbind,
                .packages = c("tidyverse","scales", "patchwork"),
                .options.snow = opts,
                .multicombine = T) %dopar% {
    progress(n = i) #update progressbar

    ion_lookup <- ions[i]
    #filter data for ion
    plot_data <- merged_extracted_data_unnested_final %>%
      dplyr::filter(.data$FG.Id == ion_lookup) %>%
      rowwise() %>%
      tidyr::unite(.data$R.Condition,.data$R.Replicate,
                   col = "sample",sep = ".",remove = F) %>%
      ungroup()

    #get Q-value and other scores
    plot_data_ion_meta <- meta_data_ions %>%
      dplyr::filter(.data$FG.Id == ion_lookup) %>%
      dplyr::rename(RawFileName = .data$R.FileName) %>%
      tidyr::unite(.data$R.Condition,.data$R.Replicate,
                   col = "sample",sep = ".",remove = F)

    #user condition ordering
    if(!is.null(selected_conditions)){
      plot_data <- plot_data %>%
        dplyr::mutate(R.Condition = factor(.data$R.Condition,
                                           levels = selected_conditions))
      plot_data_ion_meta <- plot_data_ion_meta %>%
        dplyr::mutate(R.Condition = factor(.data$R.Condition,
                                           levels = selected_conditions))
    }


    #get MS level plot data
    plot_data_MS1 <- plot_data %>%
      dplyr::filter(.data$MSLevel==1)
    plot_data_MS2 <- plot_data %>%
      dplyr::filter(.data$MSLevel==2)


    # Q-value plot
    Qvalue_plot<- ggplot(data = plot_data_ion_meta,
                         mapping = aes(x = 1,
                                       y = 1,
                                       fill= -log10(as.numeric(.data$EG.Qvalue))
                         )
    )+
      geom_tile()+
      theme_void()+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,nrow = 1)+
      theme(strip.text = element_blank(),
            axis.text.y = element_text(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0),units = "cm"),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(size = 9,hjust = 0.5, face = "bold"),
            legend.text = element_text(angle = 45,
                                       hjust = 1,
                                       vjust = 1.3),
            legend.key.width = unit(0.6, "cm"),
            legend.key.height = unit(0.25, "cm"))+
      labs(fill = "EG.Qvalue")+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("EG.Qvalue"))+
      geom_text(mapping = aes(1,1,label = .data$EG.Qvalue_formatted))+
      scale_fill_gradientn(colours = c("#E13E2E","#D97B2E","#AFC297","#7CB49B","#6BAF9F","#5EA7A5","#56A0AA","#5597AE"),
                           breaks = -log10(c(0.05,0.01,1E-3,1E-3,1E-4,1E-5,1E-6,1E-7)),
                           labels = c(0.05,0.01,1E-3,1E-3,1E-4,1E-5,1E-6,1E-7),
                          limits = c(1.3,7),
                          oob = squish,
                          na.value = "white")

    # EG.PEP plot
    PEP_plot<- ggplot(data = plot_data_ion_meta,
                         mapping = aes(x = 1,
                                       y = 1,
                                       fill= -log10(as.numeric(.data$EG.PEP))
                         )
    )+
      geom_tile()+
      theme_void()+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,nrow = 1)+
      theme(strip.text = element_blank(),
            axis.text.y = element_text(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0),units = "cm"),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(size = 9,hjust = 0.5, face = "bold"),
            legend.text = element_text(angle = 45,
                                       hjust = 1,
                                       vjust = 1.3),
            legend.key.width = unit(0.6, "cm"),
            legend.key.height = unit(0.25, "cm"))+
      labs(fill = "EG.PEP")+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("EG.PEP"))+
      geom_text(mapping = aes(1,1,label = .data$EG.PEP_formatted))+
      scale_fill_gradientn(colours = c("brown","gold3","#AFC297","#7CB49B","#6BAF9F","#5EA7A5","#56A0AA","seagreen"),
                           breaks = -log10(c(0.05,0.01,1E-3,1E-3,1E-4,1E-5,1E-6,1E-7)),
                           labels = c(0.05,0.01,1E-3,1E-3,1E-4,1E-5,1E-6,1E-7),
                           limits = c(1.3,7),
                           oob = squish,
                           na.value = "white")

    # EG.SignalToNoise plot
    EG.SignalToNoise_plot<- ggplot(data = plot_data_ion_meta,
                         mapping = aes(x = 1,
                                       y = 1,
                                       fill= log2(as.numeric(.data$EG.SignalToNoise))
                         )
    )+
      geom_tile()+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 nrow = 1)+
      theme_void()+
      labs(fill = expression(bold(log[2]~"EG.SignalToNoise")))+
      theme(strip.text = element_blank(),
            axis.text.y = element_text(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0),units = "cm"),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(size = 9,hjust = 0.5, face = "bold"),
            legend.key.width = unit(0.6, "cm"),
            legend.key.height = unit(0.25, "cm"))+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("EG.SignalToNoise"))+
      geom_text(mapping = aes(1,1,label = scales::scientific(x = as.numeric(.data$EG.SignalToNoise),
                                                             digits = 3)))+
      scale_fill_distiller(palette = "RdYlGn",
                           direction = 1,
                           breaks = seq(0,7,by = 1),
                           labels = seq(0,7,by = 1),
                           limits = c(log2(1),log2(128)),
                           oob = squish,
                           na.value = "white")

    # Shape-value plot
    ShapeValue_plot<- ggplot(data = plot_data_ion_meta,
                         mapping = aes(x = 1,
                                       y = 1,
                                       fill = .data$FG.ShapeQualityScore)
                         )+
      geom_tile()+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 nrow = 1)+
      theme_void()+
      theme(strip.text = element_blank(),
            axis.text.y = element_text(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0),units = "cm"),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(size = 9,hjust = 0.5, face = "bold"),
            legend.key.width = unit(0.6, "cm"),
            legend.key.height = unit(0.25, "cm"))+
      labs(fill = "FG.ShapeQualityScore")+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("FG.ShapeQualityScore"))+
      geom_text(mapping = aes(1,1,label = round(x = .data$FG.ShapeQualityScore,digits = 3)))+
      scale_fill_distiller(palette = "BrBG",
                           direction = 1,
                           limits = c(-1,1),
                           oob = squish,na.value = "white")
    # EG.IntCorrScore plot
    EG.IntCorrScore_plot<- ggplot(data = plot_data_ion_meta,
                             mapping = aes(x = 1,
                                           y = 1,
                                           fill= .data$EG.IntCorrScore)
    )+
      geom_tile()+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 nrow = 1)+
      theme_void()+
      theme(strip.text = element_blank(),
            axis.text.y = element_text(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0),units = "cm"),
            legend.position = "bottom",
            legend.title.position = "top",
            legend.title = element_text(size = 9,hjust = 0.5, face = "bold"),
            legend.text = element_text(angle = 45,
                                       hjust = 1,
                                       vjust = 1.3),
            legend.key.width = unit(0.6, "cm"),
            legend.key.height = unit(0.25, "cm"))+
      labs(fill = "EG.IntCorrScore")+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("EG.IntCorrScore"))+
      geom_text(mapping = aes(1,1,label = round(x = .data$EG.IntCorrScore,digits = 3)))+
      scale_fill_gradientn(colours = c("#92442B","#DA8040","#CCCECC","#96B3CE","#7298C0","#4D719E","#436692","#385B86"),
                           breaks = seq(0.3,1.0,by = 0.1),
                           labels = scales::number_format(accuracy = 0.1),
                           limits = c(0.3,1.0),
                           oob = squish,
                           na.value = "white")
    # Cscore plot
    Cscore_plot<- ggplot(data = plot_data_ion_meta,
                             mapping = aes(x = 1,
                                           y = 1,
                                           label = round(x = .data$EG.Cscore,
                                                         digits = 3))
    )+
      geom_tile(fill = "white")+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 nrow = 1)+
      theme_void()+
      theme(strip.text = element_blank(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0),units = "cm"),
            axis.text.y = element_text())+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("EG.Cscore"))+
      geom_text(fontface = "bold")

    # EG.Identified plot
    EG.Identified_plot<- ggplot(data = plot_data_ion_meta %>%
                                  dplyr::mutate(EG.Identified_fontface = ifelse(test = .data$EG.Identified==TRUE,
                                                                                yes = "bold",
                                                                                no = "plain")),
                         mapping = aes(x = 1,y = 1,label = .data$EG.Identified)
    )+
      geom_tile(fill = "white")+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 nrow = 1)+
      theme_void()+
      theme(strip.text = element_blank(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0),units = "cm"),
            axis.text.y = element_text())+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("EG.Identified"))+
      geom_text(mapping = aes(color = .data$EG.Identified,
                              fontface = .data$EG.Identified_fontface))+
      guides(color = "none", )+
      scale_color_manual(values = c("TRUE" = "dodgerblue4","FALSE" = "firebrick4"))

    # dppp MS1
    dpppMS1_plot<- ggplot(data = plot_data_ion_meta %>%
                            dplyr::mutate(dpppMS1_fontface = ifelse(test = .data$`EG.DatapointsPerPeak (MS1)`>=5,
                                                             yes = "bold",
                                                             no = "plain")),
                         mapping = aes(x = 1,y = 1,label = .data$`EG.DatapointsPerPeak (MS1)`)

    )+
      geom_tile(fill = "white")+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 nrow = 1)+
      theme_void()+
      theme(strip.text = element_blank(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                        yes = 1.8,
                                        no = 0),units = "cm"),
            axis.text.y = element_text())+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("dppp (MS1)"))+
      geom_text(mapping = aes(fontface = .data$dpppMS1_fontface))

    # dppp MS2
    dpppMS2_plot<- ggplot(data = plot_data_ion_meta %>%
                            dplyr::mutate(dpppMS2_fontface = ifelse(test = .data$EG.DatapointsPerPeak>=5,
                                                             yes = "bold",
                                                             no = "plain")),
                          mapping = aes(x = 1,y = 1, label = .data$EG.DatapointsPerPeak)

    )+
      geom_tile(fill = "white")+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 nrow = 1)+
      theme_void()+
      theme(strip.text = element_blank(),
            panel.spacing = unit(x = ifelse(test = run_specific_y_axis == TRUE,
                                            yes = 1.8,
                                            no = 0.02),units = "cm"),
            axis.text.y = element_text())+
      scale_y_continuous(position = "right",breaks = c(1),labels = c("dppp (MS2)"))+
      geom_text(mapping = aes(fontface = .data$dpppMS2_fontface))

    # MS1 plot
    tmp_ms1_plot<- ggplot(data = plot_data_MS1,
                          mapping = aes(x = .data$RT,y = .data$XIC, color = .data$IonLabel))+
      geom_vline(data = plot_data_ion_meta,
                 mapping = aes(xintercept = .data$EG.ApexRT),linetype = "dashed", color = "grey")+
      geom_line()+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 scales = ifelse(test = run_specific_y_axis == TRUE,
                                 yes = "free_y",
                                 no = "fixed"),
                 nrow = 1)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
            plot.title = element_text(face = "italic",hjust = 0.5),
            legend.title = element_text(size = 9,hjust = 0.5, face = "bold"))+
      labs(title = "XIC MS1 level", color = "XIC MS1", y = "MS1 intensity")+
      guides(color=guide_legend(ncol = 1,
                                override.aes = list(linewidth = 3)))+
      scale_y_continuous(labels = scientific_10)

    if(Spectronaut_colors==TRUE){
      tmp_ms1_plot <- tmp_ms1_plot+
        scale_color_manual(values = Spectronaut_color_fn(length(unique(plot_data_MS1$IonLabel))))
    }else{
      tmp_ms1_plot <- tmp_ms1_plot+
        scale_color_manual(values = MS1_color(length(unique(plot_data_MS1$IonLabel))))
    }


    # MS2 plot
    tmp_ms2_plot<- ggplot(data = plot_data_MS2,
                          mapping = aes(x = .data$RT, y = .data$XIC, color = .data$IonLabel))+
      geom_vline(data = plot_data_ion_meta,
                 mapping = aes(xintercept = .data$EG.ApexRT),linetype = "dashed", color = "grey")+
      geom_line()+
      facet_wrap(~ .data$R.Condition + .data$R.Replicate,
                 scales = ifelse(test = run_specific_y_axis == TRUE,
                                 yes = "free_y",
                                 no = "fixed"),
                 nrow = 1)+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),,
            plot.title = element_text(face = "italic",hjust = 0.5),
            legend.title = element_text(size = 9,hjust = 0.5, face = "bold"))+
      labs(title = "XIC MS2 level",
           color = "XIC MS2",
           y = "MS2 intensity")+
      guides(color=guide_legend(ncol = 2,
                                override.aes = list(linewidth = 3)
                                )
             )+
      scale_y_continuous(labels = scientific_10)

    if(Spectronaut_colors==TRUE){
      tmp_ms2_plot <- tmp_ms2_plot+
        scale_color_manual(values = Spectronaut_color_fn(length(unique(plot_data_MS2$IonLabel))))
    }else{
      tmp_ms2_plot <- tmp_ms2_plot+
        scale_color_manual(values = MS2_color(length(unique(plot_data_MS2$IonLabel))))
    }

    # merge to final plot
    tmp_ion_plot_out<- Qvalue_plot/
                        PEP_plot/
                        ShapeValue_plot/
                        EG.IntCorrScore_plot/
                        EG.SignalToNoise_plot/
                        Cscore_plot/
                        EG.Identified_plot/
                        dpppMS1_plot/
                        dpppMS2_plot/
                        tmp_ms1_plot/
                        tmp_ms2_plot+
      plot_layout(heights = c(0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,1,1),guides = 'collect')+
      plot_annotation(title = str_wrap(paste(unique(plot_data$PG.ProteinGroups),"-",ion_lookup),width = 60),
                      theme = theme(legend.position = 'bottom',
                                    plot.title = element_text(size = 25,face = "bold")))


    #save plot
    ggsave(plot = tmp_ion_plot_out,filename = str_replace_all(string = paste0(unique(plot_data$PG.ProteinGroups),
                                                                              "_",
                                                                              ion_lookup,
                                                                              ".pdf"),
                                              pattern = ";",
                                              replacement = "_"
                                              ),
           path = output_path,
           device = "pdf",
           limitsize = FALSE,
           width = if((1.5*sample_length)<11){
             13}else{
               ifelse(test = run_specific_y_axis == TRUE,
                      yes = 1.8,
                      no = 1.5)*sample_length},
           height = 9.2)


    #save ion data
    if(export_csv_files==TRUE){
      write_csv(x = plot_data,
                file = paste0(output_path,
                              "/",
                              str_replace_all(string = paste0(unique(plot_data$PG.ProteinGroups),"_",
                                                              ion_lookup,".csv"),
                                              pattern = ";",
                                              replacement = "_"
                              )
                )
      )
    }




    }# end loop
  parallel::stopCluster(cl)#closing cluster for parallel computing
  message("")
  message("DONE ...")

  # close connection
  rm(tmp_data_input);rm(meta_data_ions)
}









