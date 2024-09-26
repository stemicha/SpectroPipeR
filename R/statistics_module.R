#' SpectroPipeR: statistics module
#' @description
#' Function for performing statistical analysis on Spectronaut output reports, representing the fourth step in the pipeline and is dependent on the normalization and quantification module output from step 2.
#' @param SpectroPipeR_data_quant it is the SpectroPipeR_data_quant list object from norm_quant_module() object e.g. `SpectroPipeR_data_quant` see example below
#' @param condition_comparisons condition comparisons for pairwise- comparison; e.g. condition_comparisons <- cbind(c("condition1","control"),c("condition3","control") )
#' @param number_of_cores number of processor cores to be used for the calculations default = 2;
#' `parallel::detectCores()-2` for faster processing (will detect the number of cores in the system and use nearly all cores)
#' @returns SpectroPipeR_statistics list element containing the statistics analysis results in addition to the automatically generated plots and tables in output folder
#'  For the description of the generated figures and tables please read the manual & vignettes
#'
#' | <u> __list element__ </u> | <u> __description__ </u>              |
#' |:--------------------------|:--------------------------------------|
#' | stat_results              | *tibble:* statistical analysis results table |
#' | stat_column_description   | *tibble:* statistical analysis results table column description |
#' | stats_results_iBAQ_quantiles | *tibble:* statistical analysis results table containing|
#' |                              | the iBAQ quantilies (Q1-Q10) of the protein per group for a better |
#' |                              | ratio judgement |
#' | stat_results_filtered      | *tibble:* filtered (user defined FC and adj. p-value) statistical |
#' |                            |analysis results table |
#'
#' @md
#' @export
#'
#' @import ggplot2
#' @import tidyverse
#' @importFrom dplyr filter
#' @import readr
#' @import tidyr
#' @import PECA
#' @importFrom stats median sd lm as.formula quantile
#' @import ggrepel
#' @import openxlsx
#' @importFrom methods is
#' @importFrom scales squish
#' @import effsize
#' @import readxl
#' @import graphics
#' @import foreach
#' @import parallel
#' @import doSNOW
#' @import tibble
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
#'
#'# step 2: normalize & quantification module
#'SpectroPipeR_data_quant <- norm_quant_module(SpectroPipeR_data = SpectroPipeR_data)
#'
#'# step 3: MVA module
#'SpectroPipeR_data_MVA <- MVA_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
#'           HCPC_analysis = FALSE)
#'
#'# step 4: statistics module
#'SpectroPipeR_data_stats <- statistics_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
#'                                             condition_comparisons = cbind(c("HYE mix A",
#'                                                                             "HYE mix B")
#'                                                                           )
#'                                             )
#'}

statistics_module <- function(SpectroPipeR_data_quant = NULL,
                              condition_comparisons = NULL,
                              number_of_cores = 2){

  # check if SpectroPipeR_norm_quant object was provided
  if(!is(SpectroPipeR_data_quant, "SpectroPipeR_norm_quant")){
    stop_text <- paste("no valid SpectroPipeR_norm_quant object was provided: please process data with norm_quant_module() to generate one!")
    message_function(text = stop_text,color = "red",log_file_name = log_file_name)
    stop()
  }

  #output folder
  out_folder <- SpectroPipeR_data_quant$parameter$output_folder

  #log file name
  log_file_name <- SpectroPipeR_data_quant$parameter$log_file_name

  message_function(text = "",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "# STATISTICS MODULE",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "",
                   color = "white",
                   log_file_name = log_file_name)

  # if comparison is NULL
  if(is.null(condition_comparisons)){
    message_function(text = "no condition comparisons provided -- skipping statistics",
                     color = "blue",
                     log_file_name = log_file_name)

    output_list <-  NULL
    return(output_list)

  }else{
    #parameter file
    parameter <- SpectroPipeR_data_quant$parameter

    # test if all conditions are in input
    if(sum(unique(SpectroPipeR_data_quant$data_input_normalized$R.Condition)%in%
           as.character(condition_comparisons))!=length(unique(as.character(condition_comparisons)))
    ){
      stop_text <- paste("STOP: not all conditions in condition_comparisons are covered by Spectronaut report")
      message_function(text = stop_text,color = "red",log_file_name = log_file_name)
      stop()
    }

    # test if duplicates exist in condition file
    user_condition_comparison_input <- apply(condition_comparisons,2,function(x) paste(x,collapse = "_vs_"))
    if(sum(duplicated(user_condition_comparison_input))>0){
      stop_text <- paste("There are duplicated entrie(s) in the condition_comparison input!: number",which(duplicated(user_condition_comparison_input)), "is duplicated")
      message_function(text = stop_text,color = "red",log_file_name = log_file_name)
      stop()
    }

    #test if condition comparison possess only 1 replicate >> should result in skipping

    replicate_condition_count_n1<- SpectroPipeR_data_quant$data_input_normalized %>%
                                    dplyr::distinct(.data$R.Condition,.data$R.Replicate) %>%
                                    dplyr::group_by(.data$R.Condition) %>%
                                    dplyr::summarise(count = n_distinct(.data$R.Replicate)) %>%
                                    dplyr::filter(count == 1)

    # if condition has only 1 replicate !!! skip comparison in statistics
    if(nrow(replicate_condition_count_n1)>=1){
      message_function(text = paste0("condition with only one replicate detected !!!: ", paste0(replicate_condition_count_n1$R.Condition,collapse = "; ")),
                       color = "yellow",
                       log_file_name = log_file_name)

      condition_comparisons_input <- condition_comparisons
      #merge comparisions
      combined_condition_comp<- apply(condition_comparisons,2,function(x){ paste0(x[1],"/",x[2]) })

      # condition comparison removed (N = 1)
      cond_comp_rm_pos<- grepl(pattern = replicate_condition_count_n1$R.Condition,
               x = combined_condition_comp)
      cond_comp_rm <- combined_condition_comp[cond_comp_rm_pos]

      message_function(text = paste0("removing these comparisons because one condition contains only 1 replicate !!!: ", paste0(cond_comp_rm,collapse = "; ")),
                       color = "yellow",
                       log_file_name = log_file_name)

      # filter condition comparison
      condition_comparisons<- condition_comparisons_input[, which(cond_comp_rm_pos==FALSE),drop=FALSE]
      if(ncol(condition_comparisons)==0){
        condition_comparisons <- NULL
      }
    }

    #if no comparison is possible due to replicate issue
 if(is.null(condition_comparisons)){
   message_function(text = "no condition comparisons can be calculated due to low number of replicates !!!",
                    color = "yellow",
                    log_file_name = log_file_name)

   output_list <-  NULL
   return(output_list)
 }else{
    #number of samples detected  ====
    sample_length <- length(unique(SpectroPipeR_data_quant$data_input_normalized$R.FileName))

    #create file number specific folder ====
    if(dir.exists(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"))){
      message_function(text = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"," ATTENTION !!! ---> folder already exists - files will be replaced !!!"),
                       color = "magenta",
                       log_file_name = log_file_name)
    }else{
      dir.create(paste0(out_folder,"/","06_statistics/"),recursive = T,showWarnings = F)
      dir.create(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"),recursive = T,showWarnings = F)
    }

    #reformat data for PECA ====
    message_function(text = "reformatting data ...",
                     color = "blue",
                     log_file_name = log_file_name)

    data_ROPECA_in <- SpectroPipeR_data_quant$peptide_intensity %>%
      unite(.data$R.Condition,
            .data$R.Replicate,
            col = "R.Condition__R.Replicate",
            sep = "__",
            remove = F) %>%
      dplyr::select(.data$PG.ProteinGroups,
                    .data$EG.ModifiedPeptide,
                    .data$PEP.StrippedSequence,
                    .data$R.Condition__R.Replicate,
                    .data$peptide_intensity) %>%
      group_by(.data$PG.ProteinGroups,
               .data$PEP.StrippedSequence,
               .data$R.Condition__R.Replicate) %>%
      # sum up intensities of one peptide over different modifications version of this peptide
      dplyr::summarise(peptide_intensity = sum(.data$peptide_intensity)) %>%
      ungroup() %>%
      pivot_wider(names_from = .data$R.Condition__R.Replicate,
                  values_from = .data$peptide_intensity)


    #habilliage ====
    habil <- SpectroPipeR_data_quant$peptide_intensity %>%
      unite(.data$R.Condition,
            .data$R.Replicate,
            col = "R.Condition__R.Replicate",
            sep = "__",
            remove = F) %>%
      dplyr::distinct(.data$R.Condition,
                      .data$R.Condition__R.Replicate,
                      .data$R.Replicate) %>%
      dplyr::arrange(.data$R.Condition)

    #Start a cluster ====
    message_function(text = "register processor cores ...",
                     color = "blue",
                     log_file_name = log_file_name)
    cl <- parallel::makeCluster(number_of_cores, setup_strategy = "sequential") # number of cores.
    doSNOW::registerDoSNOW(cl) #register cores

    message_function(text = "performing statistical analysis (this might take a while) ...",
                     color = "blue",
                     log_file_name = log_file_name)
    start.time <- Sys.time()
    #setup the progressbar  ====
    pb <- txtProgressBar(max = ncol(condition_comparisons), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    #do calculation in parallel ====

    #_________________________________________
    # Run PECA
    #_________________________________________

    results_peca <- foreach(i=1:ncol(condition_comparisons), .combine=rbind, .options.snow = opts,.multicombine = T) %dopar% { #use the foreach function in parallel processing mode (dopar)
      setTxtProgressBar(pb, i) #update progressbar
      #lookup group1
      group1 <- as.character(habil$R.Condition__R.Replicate[habil$R.Condition %in% condition_comparisons[1,i]])
      #lookup group2
      group2 <- as.character(habil$R.Condition__R.Replicate[habil$R.Condition %in% condition_comparisons[2,i]])
      #calculate ROPECA
      if(length(group1)>=2 & length(group2)>=2){ # at least 2 replicates each !!! checking !!!
        peca.out <- PECA::PECA_df(df = as.data.frame(data_ROPECA_in),#input data.frame
                                  id =  "PG.ProteinGroups",#set the protein ID column
                                  samplenames1 =  group1,#string input for the group1
                                  samplenames2 =  group2,#string input for the group2
                                  test = parameter$stat_test,#use the ROTS test
                                  type = parameter$type_slr,#do character string indicating whether ("median") or ("tukey") is used when calculating gene/protein values
                                  normalize = FALSE,#no normalization since normalization should be done upfront (e.g. median norm.)
                                  paired = parameter$paired,
                                  progress = F)
        peca.out$PG.ProteinGroups <- rownames(peca.out) #add the rownames as column
        peca.out$group1 <- condition_comparisons[1,i]#add group1 from the pairwise comparison
        peca.out$group2 <- condition_comparisons[2,i]#add group2 from the pairwise comparison
        peca.out$slr_ratio_meta <- paste(condition_comparisons[1,i],"/",condition_comparisons[2,i],sep="")
        peca.out$test <- parameter$stat_test
        peca.out$type <- parameter$type_slr
        peca.out# output for the for loop
      }

    }
    message("")#status bar breakline mimic
    message_function(text = "close cores ...",
                     color = "blue",
                     log_file_name = log_file_name)
    close(pb)#closing progressbar
    parallel::stopCluster(cl)#closing cluster for parallel computing
    end.time <- Sys.time()

    message_function(text = paste("start to end time comparison for stat. analysis:",difftime(time1 = end.time,time2 = start.time,units = "hours"),"hours"),
                     color = "blue",
                     log_file_name = log_file_name)
    #perfoming effect size calculation on log data ====
    # using cohens.d
    message_function(text = "estimating effect sizes ...",
                     color = "blue",
                     log_file_name = log_file_name)
    effect_size_data_in<- SpectroPipeR_data_quant$peptide_intensity %>%
      unite(.data$R.Condition,
            .data$R.Replicate,
            col = "R.Condition__R.Replicate",
            sep = "__",
            remove = F) %>%
      dplyr::select(.data$PG.ProteinGroups,
                    .data$R.Condition,
                    .data$EG.ModifiedPeptide,
                    .data$PEP.StrippedSequence,
                    .data$R.Condition__R.Replicate,
                    .data$peptide_intensity) %>%
      group_by(.data$PG.ProteinGroups,
               .data$PEP.StrippedSequence,
               .data$R.Condition,
               .data$R.Condition__R.Replicate) %>%
      # sum up intensities of one peptide over different modifications version of this peptide
      dplyr::summarise(peptide_intensity = sum(.data$peptide_intensity))

    # calculating effect sizes
    effect_sizes<- effect_size_estimation(data = effect_size_data_in,
                                          condition_comparisons = condition_comparisons)


    #convert df to tibble ====
    message("") #status bar breakline mimic
    message_function(text = "join and tidy tables ...",
                     color = "blue",
                     log_file_name = log_file_name)

    results_peca <- tibble::as_tibble(results_peca) %>%
      dplyr::mutate(significant_changed = ifelse(test = .data$slr>=log2(parameter$fold_change) &
                                                   .data$p.fdr<=parameter$p_value_cutoff,
                                                 yes = "up",
                                                 no = ifelse(test = .data$slr <= log2(1/parameter$fold_change) &
                                                               .data$p.fdr<=parameter$p_value_cutoff,
                                                             yes = "down",
                                                             no = "none")),
                    significant_changed_raw_p = ifelse(test = .data$slr >= log2(parameter$fold_change) &
                                                         .data$p<=parameter$p_value_cutoff,
                                                       yes = "up",
                                                       no = ifelse(test = .data$slr<=log2(1/parameter$fold_change) &
                                                                     .data$p<=parameter$p_value_cutoff,
                                                                   yes = "down",
                                                                   no = "none")),
                    significant_changed_fc = parameter$fold_change,
                    significant_changed_p_value = parameter$p_value_cutoff) %>%
      dplyr::mutate(fold_change_absolute = 2^abs(.data$slr)) %>%
      dplyr::mutate(fold_change_direction = ifelse(test = .data$slr==0,
                                                   yes = "none",
                                                   no = ifelse(test = .data$slr>0,
                                                               yes = "up",
                                                               no = "down"))) %>%
      dplyr::mutate(fold_change = ifelse(test = .data$slr<0,
                                         yes = .data$fold_change_absolute*-1,
                                         no = .data$fold_change_absolute))

    #merge effect size tibble to peca tibble ====

    results_peca <- left_join(results_peca,
                              effect_sizes,
                              by = c("PG.ProteinGroups","slr_ratio_meta"))

    #filtering stats ====
    message_function(text = "filtering statistical table using supplied cutoffs ...",
                     color = "blue",
                     log_file_name = log_file_name)
    results_peca_filtered <- results_peca %>%
      dplyr::filter(abs(.data$slr)>=log2(parameter$fold_change) &
                      .data$p.fdr<=parameter$p_value_cutoff)
    results_peca_filtered_raw_p <- results_peca %>%
      dplyr::filter(abs(.data$slr)>=log2(parameter$fold_change) &
                      .data$p<=parameter$p_value_cutoff)

    #writing outputs ====
    message_function(text = "writing output files ...",
                     color = "blue",
                     log_file_name = log_file_name)
    write_csv(x = results_peca,
              file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis.csv"))
    write_csv(x = results_peca_filtered,
              file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_filtered__fc_",parameter$fold_change,"__adjusted_p_value_",parameter$p_value_cutoff,".csv"))
    write_csv(x = results_peca_filtered %>%
                dplyr::filter(.data$n>=2),
              file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_filtered__fc_",parameter$fold_change,"__adjusted_p_value_",parameter$p_value_cutoff,"_2_more_peptides_per_protein.csv"))
    write_csv(x = results_peca_filtered_raw_p,
              file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_filtered__fc_",parameter$fold_change,"__raw_p_value_",parameter$p_value_cutoff,".csv"))
    write_csv(x = results_peca_filtered_raw_p %>%
                dplyr::filter(.data$n>=2),
              file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_filtered__fc_",parameter$fold_change,"__raw_p_value_",parameter$p_value_cutoff,"_2_more_peptides_per_protein.csv"))

    #add iBAQ quantiles to outputs ====
    message_function(text = "adding iBAQ quantiles to statistics table ...",
                     color = "blue",
                     log_file_name = log_file_name)
    iBAQ_intensities_summary <- SpectroPipeR_data_quant$iBAQ_intensities_summary

    # add iBAQ quantiles
    results_peca_iBAQ <- left_join(results_peca,
                                   iBAQ_intensities_summary %>%
                                     dplyr::select(.data$R.Condition,
                                                   .data$PG.ProteinGroups,
                                                   .data$mean_iBAQ_intensities,
                                                   .data$iBAQ_quantiles) %>%
                                     dplyr::rename(group1 = .data$R.Condition,
                                                   group1__mean_iBAQ = .data$mean_iBAQ_intensities,
                                                   group_1__iBAQ_quantiles = .data$iBAQ_quantiles),
                                   by = c("PG.ProteinGroups", "group1")
    )

    results_peca_iBAQ <- left_join(results_peca_iBAQ,
                                   iBAQ_intensities_summary %>%
                                     dplyr::select(.data$R.Condition,
                                                   .data$PG.ProteinGroups,
                                                   .data$mean_iBAQ_intensities,
                                                   .data$iBAQ_quantiles) %>%
                                     dplyr::rename(group2 = .data$R.Condition,
                                                   group2__mean_iBAQ = .data$mean_iBAQ_intensities,
                                                   group_2__iBAQ_quantiles = .data$iBAQ_quantiles),
                                   by = c("PG.ProteinGroups", "group2")
    )
    #generate iBAQ quantile comparison column
    results_peca_iBAQ <- results_peca_iBAQ %>%
      rowwise() %>%
      dplyr::mutate(iBAQ_quantile_comp = paste0(.data$group_1__iBAQ_quantiles,
                                                "/",
                                                .data$group_2__iBAQ_quantiles)) %>%
      dplyr::mutate(iBAQ_quantile_comp = ifelse(test = is.na(.data$group_1__iBAQ_quantiles) |
                                                  is.na(.data$group_2__iBAQ_quantiles),
                                                yes = NA,
                                                no = .data$iBAQ_quantile_comp)
      )

    write_csv(x = results_peca_iBAQ,
              file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_iBAQ_added.csv"))


    #generating wide Excel outputs ====
    message_function(text = "generating Excel outputs ...",
                     color = "blue",
                     log_file_name = log_file_name)
    #load quantile color gradient
    #quantile_color_gradient <- read_excel("single_module_R_scripts/iBAQ_dynamic_range_2D-Color_gradient.xlsx")
    quantile_color_gradient <- readxl::read_excel(path = system.file("extdata", "iBAQ_dynamic_range_2D_Color_gradient.xlsx", package="SpectroPipeR"))
    quantile_color_gradient_tidy <- quantile_color_gradient %>%
      pivot_longer(cols = colnames(quantile_color_gradient)[-1],names_to = "group1",values_to = "colors")
    colnames(quantile_color_gradient_tidy)[1] <- "group2"
    quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>% select(group1,group2,colors)
    quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>% mutate(group1 = str_replace_all(group1,"condition1_",""))
    quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>% mutate(group2 = str_replace_all(group2,"condition2_",""))
    quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>%
      rowwise() %>%
      dplyr::mutate(iBAQ_quantile_comp = paste0(group1,"/",group2)) %>%
      ungroup()
    #generate quantile color vector
    quantile_color_gradient_vector <- quantile_color_gradient_tidy$colors
    names(quantile_color_gradient_vector) <- quantile_color_gradient_tidy$iBAQ_quantile_comp

    #iBAQ quantile ratio comparison color legend plot
    iBAQ_quantile_ratio_comp_legend<- ggplot(data = quantile_color_gradient_tidy %>%
                                               dplyr::mutate(group1 = factor(.data$group1,
                                                                             levels = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10")),
                                                             group2 = factor(.data$group2,
                                                                             levels = c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10"))),
                                             mapping=aes(x = .data$group1,
                                                         y = .data$group2,
                                                         fill = .data$iBAQ_quantile_comp))+
      geom_tile()+
      scale_fill_manual(values = quantile_color_gradient_vector)+
      guides(fill="none")+
      labs(title = "iBAQ quantile ratio comparison color legend")+
      theme_light(base_size = 12)+
      theme(axis.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"))

    #save cutoff plot
    ggsave_pdf_png(plot = iBAQ_quantile_ratio_comp_legend,
                   width = 5,
                   height = 5,
                   filename = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/iBAQ_quantile_ratio_comparison_legend")
    )

    # filter statistics for at least 2 peptides
    results_peca_n2 <- results_peca %>%
      dplyr::filter(.data$n>=2)
    results_peca_iBAQ_n2 <- results_peca_iBAQ %>%
      dplyr::filter(.data$n>=2) %>%
      dplyr::relocate(.data$iBAQ_quantile_comp,
                      .after = .data$slr)

    #joining spread/wide tables (slr, p, p.fdr) !!!! n>=2 !!!!! min. 2 peptides for a wide excel table
    excel_stats<- full_join(
      full_join(results_peca_n2 %>%
                  dplyr::select(.data$slr,
                                .data$PG.ProteinGroups,
                                .data$slr_ratio_meta) %>%
                  dplyr::mutate(slr_ratio_meta = paste(.data$slr_ratio_meta,"|signal_log2_ratio",sep="")) %>%
                  pivot_wider(names_from = .data$slr_ratio_meta,
                              values_from = .data$slr),
                results_peca_n2 %>%
                  dplyr::select(.data$p,
                                .data$PG.ProteinGroups,
                                .data$slr_ratio_meta) %>%
                  dplyr::mutate(slr_ratio_meta = paste(.data$slr_ratio_meta,"|raw_p_value",sep="")) %>%
                  pivot_wider(names_from = .data$slr_ratio_meta,
                              values_from = .data$p),
                by="PG.ProteinGroups"),
      results_peca_n2 %>%
        dplyr::select(.data$p.fdr,
                      .data$PG.ProteinGroups,
                      .data$slr_ratio_meta) %>%
        dplyr::mutate(slr_ratio_meta = paste(.data$slr_ratio_meta,"|adjusted_p_value",sep="")) %>%
        pivot_wider(names_from = .data$slr_ratio_meta,
                    values_from = .data$p.fdr),
      by="PG.ProteinGroups")

    #reordering in format suggest by Elke Hammer
    column_ordering_excel <- order(as.character(sapply(colnames(excel_stats)[-1],function(x)unlist(strsplit(x = x,split = "|",fixed=T))[1])))
    excel_stats_ordered <- excel_stats[,c(1,column_ordering_excel+1)]

    #write wide data
    write_csv(x = excel_stats_ordered,
              file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_WIDE_FORMAT_2_more_peptides_per_protein.csv"))


    #EXCEL: create statitsics table iBAQ quantiles workbook ====
    #conditional formating of adj.p-value
    pValueStyle_bad_iBAQ <- openxlsx::createStyle(fontColour = "grey40")
    pValueStyle_good_iBAQ <- openxlsx::createStyle(fontColour = "black", bgFill = "#F48FB1")
    raw_pValueStyle_bad_iBAQ <- openxlsx::createStyle(fontColour = "grey40")
    raw_pValueStyle_good_iBAQ <- openxlsx::createStyle(fontColour = "black")

    excel_output_iBAQ_stats <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb = excel_output_iBAQ_stats, sheetName = "statistics__min_2_Peptides")
    #write data into workbook
    openxlsx::writeDataTable(wb = excel_output_iBAQ_stats,
                             sheet =  1,
                             x = results_peca_iBAQ_n2,
                             startRow = 1,
                             startCol = 1,
                             tableStyle = "TableStyleMedium1")
    #conditional formating of slr
    openxlsx::conditionalFormatting(wb = excel_output_iBAQ_stats,
                                    sheet = 1,
                                    cols = which(colnames(results_peca_iBAQ_n2)=="slr"),
                                    rows = 1:nrow(results_peca_iBAQ_n2)+1,
                                    style = c("dodgerblue3", "white","orangered3"),
                                    rule = c(-2, 0, 2),
                                    type = "colourScale")
    #conditional formating of iBAQ quantile comparison
    for(i in 1:length(quantile_color_gradient_vector)){
      tmp_Style <- openxlsx::createStyle(bgFill = quantile_color_gradient_vector[i])

      openxlsx::conditionalFormatting(wb = excel_output_iBAQ_stats,
                                      sheet = 1,
                                      cols = which(colnames(results_peca_iBAQ_n2)=="iBAQ_quantile_comp"),
                                      rows = 1:nrow(results_peca_iBAQ_n2)+1,
                                      style = tmp_Style,
                                      rule = names(quantile_color_gradient_vector)[i],
                                      type = "contains")

    }
    #q-value formatting
    openxlsx::conditionalFormatting(wb = excel_output_iBAQ_stats,
                                    sheet = 1,
                                    cols = which(colnames(results_peca_iBAQ_n2)=="p.fdr"),
                                    rows = 1:nrow(results_peca_iBAQ_n2)+1,
                                    style = pValueStyle_good_iBAQ,
                                    rule = "<=0.05")
    openxlsx::conditionalFormatting(wb = excel_output_iBAQ_stats,
                                    sheet = 1,
                                    cols = which(colnames(results_peca_iBAQ_n2)=="p.fdr"),
                                    rows = 1:nrow(results_peca_iBAQ_n2)+1,
                                    style = pValueStyle_bad_iBAQ,
                                    rule = ">0.05")
    #p-value formatting
    openxlsx::conditionalFormatting(wb = excel_output_iBAQ_stats,
                                    sheet = 1,
                                    cols = which(colnames(results_peca_iBAQ_n2)=="p"),
                                    rows = 1:nrow(results_peca_iBAQ_n2)+1,
                                    style = raw_pValueStyle_good_iBAQ,
                                    rule = "<=0.05")
    openxlsx::conditionalFormatting(wb = excel_output_iBAQ_stats,
                                    sheet = 1,
                                    cols = which(colnames(results_peca_iBAQ_n2)=="p"),
                                    rows = 1:nrow(results_peca_iBAQ_n2)+1,
                                    style = raw_pValueStyle_bad_iBAQ,
                                    rule = ">0.05")

    # add legend
    openxlsx::addWorksheet(wb = excel_output_iBAQ_stats, sheetName = "legend",)

    excel_output_iBAQ_stats_legend <- tibble::tribble(
      ~"column", ~"description",
      "slr","signal log2-ratios on peptide basis",
      "iBAQ_quantile_comp", "iBAQ quantiles of comparison",
      "t","t of t-statistics on peptide basis",
      "score","score of t-statistics on peptide basis",
      "n","number of peptides",
      "p","raw p-value of statistics on peptide basis",
      "p.fdr","adjusted p-value (q-value) of statistics on peptide basis",
      "PG.ProteinGroups","Protein groups",
      "group1","group1 of condition comparison",
      "group2","group2 of condition comparison",
      "slr_ratio_meta","condition comparison; how the ratio is formed",
      "test","which test was used for statistics on peptide level",
      "type","which type of ratio aggregation to ProteinGroup level was used for signal log2-ratios on peptide basis",
      "significant_changed","if there is a significant change FC & q-value (cutoffs e.g.: FC = 1.5 & adjusted-p-value = 0.05)",
      "significant_changed_raw_p","if there is a significant change FC & p-value (cutoffs e.g.: FC = 1.5 & p-value = 0.05)",
      "significant_changed_fc","fold-change cutoff used for analysis",
      "significant_changed_p_value","p-value/q-value cutoff used for analysis",
      "fold_change_absolute","ablsolute fold-change",
      "fold_change_direction","fold-change direction",
      "fold_change","fold-change",
      "effect_size_method","effect size estimation method used",
      "d","effect size estimate",
      "d_pooled_SD","effect size estimate; pooled SD",
      "d_95CI_lower","effect size estimate:the lower 95% confidence interval",
      "d_95CI_upper","effect size estimate:the upper 95% confidence interval",
      "d_magnitute","a qualitative assessment of the magnitude of effect size (|d|<0.2 negligible, |d|<0.5 small, |d|<0.8 medium, otherwise large); Cohen 1992",
      "group1__mean_iBAQ","mean iBAQ intensity of group 1",
      "group_1__iBAQ_quantiles","iBAQ intensity quantile of group 1",
      "group2__mean_iBAQ","mean iBAQ intensity of group 2",
      "group_2__iBAQ_quantiles","iBAQ intensity quantile of group 2"
    )
    #write legend data into workbook
    openxlsx::writeDataTable(wb = excel_output_iBAQ_stats,
                             sheet =  2,
                             x = excel_output_iBAQ_stats_legend,
                             startRow = 1,
                             startCol = 1,
                             tableStyle = "TableStyleMedium1")

    openxlsx::saveWorkbook(wb = excel_output_iBAQ_stats,
                           overwrite = T,
                           file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_2_more_peptides_per_protein_iBAQ_quantiles.xlsx"))


    #EXCEL: create statistics table workbook ====

    #conditional formating of adj.p-value
    pValueStyle_bad <- openxlsx::createStyle(fontColour = "lightgrey", bgFill = "white")
    pValueStyle_good <- openxlsx::createStyle(fontColour = "black", bgFill = "grey")

    excel_output <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb = excel_output, sheetName = "stat_data__min_2_Peptides")
    #write data into workbook
    openxlsx::writeDataTable(wb = excel_output,sheet =  1, x = excel_stats_ordered, startRow = 1, startCol = 1, tableStyle = "TableStyleMedium1")
    #writeDataTable(wb = excel_output,sheet =  1, x = excel_stats_ordered, startRow = 1, startCol = 1, tableStyle = "none")
    #conditional formating of slr
    slr_excel_cond_format<- grep(pattern = "signal_log2_ratio",x = colnames(excel_stats_ordered))
    for(i in 1:length(slr_excel_cond_format)){
      openxlsx::conditionalFormatting(wb = excel_output,
                                      sheet = 1,
                                      cols = slr_excel_cond_format[i],
                                      rows = 1:nrow(excel_stats_ordered)+1,
                                      style = c("dodgerblue3", "white","orangered3"),
                                      rule = c(-2, 0, 2),
                                      type = "colourScale")

    }


    p_fdr_excel_cond_format<- grep(pattern = "adjusted_p_value",x = colnames(excel_stats_ordered))
    for(i in 1:length(p_fdr_excel_cond_format)){
      openxlsx::conditionalFormatting(wb = excel_output,
                                      sheet = 1,
                                      cols = p_fdr_excel_cond_format[i],
                                      rows = 1:nrow(excel_stats_ordered)+1,
                                      style = pValueStyle_good,
                                      rule = "<=0.05")
      openxlsx::conditionalFormatting(wb = excel_output,
                                      sheet = 1,
                                      cols = p_fdr_excel_cond_format[i],
                                      rows = 1:nrow(excel_stats_ordered)+1,
                                      style = pValueStyle_bad,
                                      rule = ">0.05")
    }


    #conditional formating of raw.p-value
    p_excel_cond_format<- grep(pattern = "raw_p_value",x = colnames(excel_stats_ordered))
    for(i in 1:length(p_excel_cond_format)){
      openxlsx::conditionalFormatting(wb = excel_output,
                                      sheet = 1,
                                      cols = p_excel_cond_format[i],
                                      rows = 1:nrow(excel_stats_ordered)+1,
                                      style = pValueStyle_good,
                                      rule = "<=0.05")
      openxlsx::conditionalFormatting(wb = excel_output,
                                      sheet = 1,
                                      cols = p_excel_cond_format[i],
                                      rows = 1:nrow(excel_stats_ordered)+1,
                                      style = pValueStyle_bad,
                                      rule = ">0.05")
    }

    openxlsx::saveWorkbook(wb = excel_output,
                           overwrite = T,
                           file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_WIDE_FORMAT_2_more_peptides_per_protein.xlsx"))






    # fold_change cutoff test ====
    #___________________________________
    message_function(text = "performing fold-change cutoff sensitivity analysis ...",
                     color = "blue",
                     log_file_name = log_file_name)
    #fold-change sequence
    FC_seq <- seq(1.1,5,by=0.1)

    results_peca_FC_seq <- c()
    for(i in FC_seq){
      tmp <- c()
      tmp<- results_peca %>%
        dplyr::select(.data$slr,
                      .data$PG.ProteinGroups,
                      .data$slr_ratio_meta,
                      .data$p.fdr,
                      .data$n) %>%
        dplyr::mutate(significant_changed = ifelse(test = ((abs(.data$slr) >= log2(i)) &
                                                             .data$p.fdr <= parameter$p_value_cutoff &
                                                             .data$n >= 2),
                                                   yes = "regulated (n>=2)",
                                                   no = ifelse(test = ((abs(.data$slr) >= log2(i)) &
                                                                         .data$p.fdr <= parameter$p_value_cutoff &
                                                                         .data$n == 1),
                                                               yes = "regulated (n=1)",
                                                               no = ifelse(test = ((abs(.data$slr) < log2(i)) |
                                                                                     .data$p.fdr > parameter$p_value_cutoff &
                                                                                     .data$n >= 2),
                                                                           yes = "none (n>=2)",
                                                                           no = ifelse(test = ((abs(.data$slr) < log2(i)) |
                                                                                                 .data$p.fdr > parameter$p_value_cutoff &
                                                                                                 .data$n == 1),
                                                                                       yes = "none (n=1)",
                                                                                       no = ""))))) %>%
        ungroup() %>%
        group_by(.data$significant_changed,
                 .data$slr_ratio_meta) %>%
        dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)) %>%
        ungroup() %>%
        dplyr::mutate(FC_cut = i)
      results_peca_FC_seq <- bind_rows(results_peca_FC_seq,
                                       tmp)
    }

    number_of_proteins <- length(unique(results_peca$PG.ProteinGroups))


    message_function(text = "plotting fold-change cutoff sensitivity analysis ...",
                     color = "blue",
                     log_file_name = log_file_name)

    #unique vector for comparisons ====
    comparis <- unique(results_peca$slr_ratio_meta)

    for(i in 1:length(comparis)){
      statusprogressbar(run = i,max.run = length(comparis),width = 80L)

      fc_cutoff_plot <- results_peca_FC_seq %>%
        dplyr::filter(.data$slr_ratio_meta == comparis[i]) %>%
        ggplot(mapping = aes(x = .data$FC_cut,
                             y = .data$count,
                             fill = .data$significant_changed))+
        geom_area(alpha=0.7 , size=1, colour=NA)+
        labs(title = paste(comparis[i], "\n(in total ",number_of_proteins," proteins)",sep=""),x="fold-change cutoff",
             subtitle = paste("number of regulated protein groups (adj. p-value <",parameter$p_value_cutoff,")",sep=""),
             caption=paste("pre-definded fold-change cutoff (black solid line) =",parameter$fold_change),
             fill="regulated?\n n...number of peptides")+
        theme_minimal(base_size = 13)+
        geom_vline(xintercept = parameter$fold_change)+
        scale_fill_manual(values = c("none (n>=2)"="darkgrey","none (n=1)"="lightgrey","regulated (n>=2)"="orangered4",
                                     "regulated (n=1)"="orangered1"))+
        scale_y_continuous(sec.axis = sec_axis(~(./number_of_proteins)*100, name = "% of total"))+
        theme(axis.ticks.y = element_line(color="black"),panel.grid = element_blank())

      # work around if 0 cahnged Proteins were found in cutoff plot
      if(dim(results_peca_FC_seq %>%
             group_by(.data$slr_ratio_meta,
                      .data$FC_cut) %>%
             dplyr::filter(.data$significant_changed != "none (n=1)" &
                           .data$significant_changed != "none (n>=2)") %>%
             dplyr::summarise(count = sum(.data$count)))[1] == 0){
        # blank zoom plot
        fc_cutoff_plot_zoom <- fc_cutoff_plot +
          coord_cartesian(xlim = c(1.1,parameter$fold_change+0.5))+
          labs(title = paste(comparis[i], "\n(in total ",number_of_proteins," proteins)  [Zoomed view !]",
                             sep=""))+
          theme(legend.position="none",
                axis.ticks.y = element_line(color="black"),
                panel.grid = element_blank())

      }else{
        fc_cutoff_plot_zoom <- fc_cutoff_plot + coord_cartesian(ylim = c(0,max(results_peca_FC_seq %>%
                                                                                 group_by(.data$slr_ratio_meta,
                                                                                          .data$FC_cut) %>%
                                                                                 dplyr::filter(.data$significant_changed!="none (n=1)" &
                                                                                                 .data$significant_changed!="none (n>=2)") %>%
                                                                                 dplyr::summarise(count = sum(.data$count)) %>%
                                                                                 ungroup() %>%
                                                                                 dplyr::select(.data$count))))+
          labs(title = paste(comparis[i], "\n(in total ",number_of_proteins," proteins)  [Zoomed view !]",sep=""))+
          theme(legend.position="none",axis.ticks.y = element_line(color="black"),panel.grid = element_blank())
      }



      fc_final_plot <- fc_cutoff_plot_zoom+fc_cutoff_plot
      #save cutoff plot
      ggsave_pdf_png(plot = fc_final_plot,
                     width = 15,
                     height = 5,
                     filename = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/cutoff_test.pdf",str_replace_all(comparis[i],"/","_vs_"))
      )
    }


    message("") #status bar breakline mimic


    # simple cutoff test 2 peptides -------------------------------------------

    message_function(text = "plotting fold-change simple cutoff sensitivity analysis (peptide n > 1)...",
                     color = "blue",
                     log_file_name = log_file_name)

    fc_cuts<- seq(1.1,10,by = 0.1)
    cutoff_test <- c()
    for(i in 1:length(fc_cuts)){

      cutoff_test<- rbind(cutoff_test,bind_rows(results_peca_n2 %>%
                                                  dplyr::filter(abs(.data$slr)>log2(fc_cuts[i])) %>%
                                                  group_by(.data$slr_ratio_meta) %>%
                                                  dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)) %>%
                                                  dplyr::mutate(cuts = fc_cuts[i],
                                                                q_value = "none"),
                                                results_peca_n2 %>%
                                                  dplyr::filter(abs(.data$slr)>log2(fc_cuts[i]) &
                                                                  .data$p.fdr<0.05) %>%
                                                  group_by(.data$slr_ratio_meta) %>%
                                                  dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)) %>%
                                                  dplyr::mutate(cuts = fc_cuts[i],
                                                                q_value = "0.05")))

    }

    #pdf()
    for(i in 1:length(comparis)){
      statusprogressbar(run = i,max.run = length(comparis),width = 80L)

      tmp_cutoff_plot <- cutoff_test %>%
        dplyr::filter(.data$slr_ratio_meta==comparis[i])

      fc_labels <- tibble::tibble(cuts = c(1.5,2,4),
                                  count = max(tmp_cutoff_plot$count,na.rm = T),
                                  labels= c("absFC = 1.5","absFC = 2.0","absFC = 4.0"))

      cutoff_plot <- ggplot(data = tmp_cutoff_plot,
                            mapping = aes(x = .data$cuts,
                                          y = .data$count,
                                          color = .data$q_value))+
        geom_point()+
        scale_color_brewer(palette = "Set1")+
        theme_light(base_size = 18)+
        geom_vline(xintercept = 1.5,linetype = "solid")+
        geom_vline(xintercept = 2,linetype = "dashed")+
        geom_vline(xintercept = 4,linetype = "dotted")+
        geom_label_repel(data = fc_labels,
                         mapping = aes(label = .data$labels),
                         color = "black",
                         min.segment.length = 0.01)+
        labs(title = str_replace_all(string = comparis[i],pattern = "/",replacement = "/\n"),
             subtitle = "simple cut-off plot; peptide count > 1",
             x = "absFC cuts",
             y ="protein_groups count")

      #save cutoff plot
      ggsave_pdf_png(plot = cutoff_plot,
                     width = 8,
                     height = 5.4,
                     filename = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/cutoff_simple_test__",str_replace_all(comparis[i],"/","_vs_"))
      )


    }


    message("") #status bar breakline mimic


    # Volcanoplots ------------------------------------------------------------
    message_function(text = "generating volcano plots ...",
                     color = "blue",
                     log_file_name = log_file_name)
    #filter for peptides
    results_peca_volcano <- results_peca %>%
      dplyr::filter(.data$n>=2)  #number of peptides >2

    #workaround for 0 values in p-value
    #replace with min. p-value without 0
    results_peca_volcano$p.fdr[results_peca_volcano$p.fdr==0] <- min(results_peca_volcano$p.fdr[results_peca_volcano$p.fdr!=0])
    results_peca_volcano$p[results_peca_volcano$p==0] <- min(results_peca_volcano$p[results_peca_volcano$p!=0])

    #ADJUSTED P-VALUE Volcano plots ====

    #adjusted p-value Volcano plots
    max.slr <- max(abs(results_peca_volcano$slr))
    max.log10.pval <- max(-log10(results_peca_volcano$p.fdr))

    #do volcano plots adjusted p-value
    for(i in 1:length(comparis)){
      statusprogressbar(run = i,max.run = length(comparis),width = 60L,info = "Volcano plots with adj. p-value")
      #filter data single comp.
      tmp <- results_peca_volcano %>%
        dplyr::filter(.data$slr_ratio_meta == comparis[i])
      tmp <- tmp %>%
        dplyr::mutate(euclidean_distance = sqrt((((-log10(.data$p.fdr))^2)*(.data$slr)^2)))
      tmp_top10_each <- tmp %>%
        dplyr::filter(.data$significant_changed != "none") %>%
        group_by(.data$significant_changed) %>%
        top_n(.data$euclidean_distance,
              n = 10) %>%
        ungroup()

      #setup count data matrix
      tmp.count <- bind_rows(tibble::tibble(x_pos=if_else(max.slr<3,max.slr-0.5,max.slr-2.3),
                                            y_pos=max.log10.pval+1,
                                            significant_changed="up",
                                            count = as.numeric(tmp %>%
                                                                 dplyr::filter(.data$significant_changed == "up") %>%
                                                                 dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      ),
      tibble::tibble(x_pos=0,
                     y_pos=max.log10.pval+1,
                     significant_changed="none",
                     count = as.numeric(tmp %>%
                                          dplyr::filter(.data$significant_changed == "none") %>%
                                          dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      ),
      tibble::tibble(x_pos=if_else(max.slr<3,-max.slr+0.5,-max.slr+2.3),
                     y_pos=max.log10.pval+1,
                     significant_changed="down",
                     count = as.numeric(tmp %>%
                                          dplyr::filter(.data$significant_changed == "down") %>%
                                          dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      )
      )
      #do volcano plot
      tmp_plot<- ggplot(data = tmp,
                        mapping = aes(x = .data$slr,
                                      y = -log10(.data$p.fdr),
                                      color = .data$significant_changed))+
        geom_point(alpha=0.4)+
        theme_minimal(base_size = 18)+
        labs(title = str_replace_all(string = comparis[i],pattern = "/",replacement = "/\n"),
             subtitle = paste("cutoffs: fold-change=",parameter$fold_change," / adj. p-value=",parameter$p_value_cutoff,sep=""),
             x=expression(log[2]~ratios),
             y=expression(-"log"[10]~"adj. p-value"),
             caption = "protein data filtered for equal or more 2 peptides")+
        geom_hline(yintercept = -log10(parameter$p_value_cutoff),linetype="dotted",color="black")+
        geom_vline(xintercept = log2(parameter$fold_change),linetype="dotted",color="black")+
        geom_vline(xintercept = log2(1/parameter$fold_change),linetype="dotted",color="black")+
        scale_color_manual(values = c(down = "dodgerblue3",none = "grey",up = "orangered3"))+
        coord_cartesian(xlim = c(-max.slr,max.slr),ylim = c(0,max.log10.pval+1))+
        geom_label(data=tmp.count,
                   mapping = aes(x = .data$x_pos,
                                 y = .data$y_pos,
                                 label = .data$count,
                                 fill = .data$significant_changed),
                   inherit.aes = F,size = 9,color = "white")+
        scale_fill_manual(values = c(down = "dodgerblue3",none = "grey",up = "orangered3"))+
        guides(color=F, fill=F, shape = guide_legend(override.aes = list(size=5)))+
        theme(axis.title = element_text(size=18))

      #adding labels to Top10 proteins
      tmp_plot_label_top10_each <- tmp_plot+
        geom_text_repel(data = tmp_top10_each,
                        mapping = aes(label = .data$PG.ProteinGroups),
                        min.segment.length = 0.1)

      #save plot
      ggsave_pdf_png(plot = tmp_plot+tmp_plot_label_top10_each,
                     width = 16,
                     height = 8.4,
                     filename = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/volcano_plots_adjusted_p_value__",str_replace_all(comparis[i],"/","_vs_")))
    }

    #ADJUSTED P-VALUE Volcano plots  with effect size====

    #adjusted p-value Volcano plots
    max.slr <- max(abs(results_peca_volcano$slr))
    max.log10.pval <- max(-log10(results_peca_volcano$p.fdr))

    #do volcano plots adjusted p-value
    for(i in 1:length(comparis)){
      statusprogressbar(run = i,max.run = length(comparis),width = 60L,info = "Volcano plots with adj. p-value")
      #filter data single comp.
      tmp <- results_peca_volcano %>%
        dplyr::filter(.data$slr_ratio_meta == comparis[i])
      tmp <- tmp %>%
        dplyr::mutate(euclidean_distance = sqrt((((-log10(.data$p.fdr))^2)*(.data$slr)^2)))
      tmp_top10_each <- tmp %>%
        dplyr::filter(.data$significant_changed!="none") %>%
        group_by(.data$significant_changed) %>%
        top_n(.data$euclidean_distance,
              n = 10) %>%
        ungroup()

      #setup count data matrix
      tmp.count <- bind_rows(tibble::tibble(x_pos = if_else(max.slr<3,max.slr-0.5,max.slr-2.3),
                                            y_pos = max.log10.pval+1,
                                            significant_changed = "up",
                                            count = as.numeric(tmp %>%
                                                                 dplyr::filter(.data$significant_changed == "up") %>%
                                                                 dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      ),
      tibble::tibble(x_pos=0,
                     y_pos=max.log10.pval+1,
                     significant_changed="none",
                     count = as.numeric(tmp %>%
                                          dplyr::filter(.data$significant_changed == "none") %>%
                                          dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      ),
      tibble::tibble(x_pos=if_else(max.slr<3,-max.slr+0.5,-max.slr+2.3),
                     y_pos=max.log10.pval+1,
                     significant_changed="down",
                     count = as.numeric(tmp %>%
                                          dplyr::filter(.data$significant_changed == "down") %>%
                                          dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      )
      )
      #do volcano plot
      tmp_plot<- ggplot(data = tmp,
                        mapping = aes(x = .data$slr,
                                      y = -log10(.data$p.fdr),
                                      color = .data$significant_changed,
                                      shape = .data$d_magnitute))+
        geom_point(alpha=0.4)+
        theme_minimal(base_size = 18)+
        labs(title = str_replace_all(string = comparis[i],pattern = "/",replacement = "/\n"),
             subtitle = paste("cutoffs: fold-change=",parameter$fold_change," / adj. p-value=",parameter$p_value_cutoff,sep=""),
             x=expression(log[2]~ratios),
             y=expression(-"log"[10]~"adj. p-value"),
             shape = "Cohen's eff. size\nmagnitude:",
             caption = "protein data filtered for equal or more 2 peptides")+
        geom_hline(yintercept = -log10(parameter$p_value_cutoff),linetype="dotted",color="black")+
        geom_vline(xintercept = log2(parameter$fold_change),linetype="dotted",color="black")+
        geom_vline(xintercept = log2(1/parameter$fold_change),linetype="dotted",color="black")+
        scale_color_manual(values = c(down = "dodgerblue3",none = "grey",up = "orangered3"))+
        scale_shape_manual(values = c(negligible = 3,small = 1,medium = 16,large = 15))+
        coord_cartesian(xlim = c(-max.slr,max.slr),ylim = c(0,max.log10.pval+1))+
        geom_label(data = tmp.count,
                   mapping = aes(x = .data$x_pos,
                                 y = .data$y_pos,
                                 label = .data$count,
                                 fill = .data$significant_changed),
                   inherit.aes = F,size = 9,color = "white")+
        scale_fill_manual(values = c(down = "dodgerblue3",none = "grey",up = "orangered3"))+
        guides(color=F,
               fill=F,
               shape = guide_legend(override.aes = list(size=5)))+
        theme(axis.title = element_text(size=18),
              legend.position = "bottom",
              legend.title = element_text(size = 15, vjust = 1,hjust = 1))

      #adding labels to Top10 proteins
      tmp_plot_label_top10_each <- tmp_plot+
        geom_text_repel(data = tmp_top10_each,
                        mapping = aes(label = .data$PG.ProteinGroups),
                        min.segment.length = 0.1)

      #save plot
      ggsave_pdf_png(plot = tmp_plot +
                       tmp_plot_label_top10_each,
                     width = 16,
                     height = 9.2,
                     filename = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/volcano_plots_effect_size_shape_adjusted_p_value__",str_replace_all(comparis[i],"/","_vs_")))
    }


    #RAW P-VALUE Volcano plots ====

    #p-value Volcano plots
    max.slr <- max(abs(results_peca_volcano$slr))
    max.log10.pval <- max(-log10(results_peca_volcano$p))

    #do volcano plots
    for(i in 1:length(comparis)){
      statusprogressbar(run = i,max.run = length(comparis),width = 60L,info = "Volcano plots with raw p-value")
      #filter data single comp.
      tmp <- results_peca_volcano %>%
        dplyr::filter(.data$slr_ratio_meta==comparis[i])
      tmp <- tmp %>%
        dplyr::mutate(euclidean_distance = sqrt((((-log10(.data$p))^2)*(.data$slr)^2)))
      tmp_top10_each <- tmp %>%
        dplyr::filter(.data$p<=parameter$p_value_cutoff & abs(.data$slr)>=log2(parameter$fold_change)) %>%
        group_by(.data$fold_change_direction) %>%
        dplyr::top_n(.data$euclidean_distance,
                     n = 10) %>%
        ungroup()

      #setup count data matrix
      tmp.count <- bind_rows(tibble::tibble(x_pos = if_else(max.slr<3,max.slr-0.5,max.slr-2.3),
                                            y_pos = max.log10.pval+1,
                                            significant_changed="up",
                                            count = as.numeric(tmp %>%
                                                                 dplyr::filter(.data$significant_changed_raw_p=="up") %>%
                                                                 dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      ),
      tibble::tibble(x_pos = 0,
                     y_pos = max.log10.pval+1,
                     significant_changed = "none",
                     count = as.numeric(tmp %>%
                                          dplyr::filter(.data$significant_changed_raw_p == "none") %>%
                                          dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      ),
      tibble::tibble(x_pos = if_else(max.slr<3,-max.slr+0.5,-max.slr+2.3),
                     y_pos = max.log10.pval+1,
                     significant_changed = "down",
                     count = as.numeric(tmp %>%
                                          dplyr::filter(.data$significant_changed_raw_p == "down") %>%
                                          dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)))
      )
      )
      #do volcano plot
      tmp_plot<- ggplot(data = tmp,aes(x = .data$slr,
                                       y = -log10(.data$p),
                                       color = .data$significant_changed_raw_p))+
        geom_point(alpha=0.4)+
        theme_minimal(base_size = 18)+
        labs(title = str_replace_all(string = comparis[i],pattern = "/",replacement = "/\n"),
             subtitle = paste("cutoffs: fold-change=",parameter$fold_change," / p-value=",parameter$p_value_cutoff,sep=""),
             x=expression(log[2]~ratios),
             y=expression(-"log"[10]~"p-value"),
             caption = "protein data filtered for equal or more 2 peptides")+
        geom_hline(yintercept = -log10(parameter$p_value_cutoff),linetype="dotted",color="black")+
        geom_vline(xintercept = log2(parameter$fold_change),linetype="dotted",color="black")+
        geom_vline(xintercept = log2(1/parameter$fold_change),linetype="dotted",color="black")+
        scale_color_manual(values = c(down = "dodgerblue3",none = "grey",up = "orangered3"))+
        coord_cartesian(xlim = c(-max.slr,max.slr),ylim = c(0,max.log10.pval+1))+
        geom_label(data = tmp.count,
                   mapping = aes(x = .data$x_pos,
                                 y = .data$y_pos,
                                 label = .data$count,
                                 fill = .data$significant_changed),
                   inherit.aes = F,size = 9,color = "white")+
        scale_fill_manual(values = c(down = "dodgerblue3",none = "grey",up = "orangered3"))+
        guides(color=F, fill=F)+
        theme(axis.title = element_text(size=18))

      # add new plots for cohens d

      #adding labels to Top10 proteins
      tmp_plot_label_top10_each <- tmp_plot+
        geom_text_repel(data = tmp_top10_each,
                        mapping = aes(label=.data$PG.ProteinGroups),
                        min.segment.length = 0.1)

      #save plot
      ggsave_pdf_png(plot = tmp_plot+tmp_plot_label_top10_each,
                     width = 16,
                     height = 8.4,
                     filename = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/volcano_plots_raw_p_value__",str_replace_all(comparis[i],"/","_vs_")))

    }

    #perform ratio comparison benchmark ====
    ## use Peptide-Ratios as ground truth and compare them with Protein ratios
    ## use median protein intensity over replicates to generate protein wise ratios
    message("")

    #adaption of theme size of ratio comparison benchmark plots
    if(ncol(condition_comparisons)==1){
      theme_size = 10
    }else{
      theme_size = 18
    }

    #add detection with at least 2 peptides ====
    detection_percentage <- SpectroPipeR_data_quant$PG_2_peptides_ID_raw %>%
      dplyr::select(.data$R.Condition,
                    .data$PG.ProteinGroups,
                    .data$present_replicate_percentage)

    detection_percentage_final <- c()
    for(i in 1:ncol(condition_comparisons)){
      detection_percentage_final <- bind_rows(detection_percentage_final,
                                              full_join(detection_percentage %>%
                                                          dplyr::filter(.data$R.Condition == condition_comparisons[1,i]) %>%
                                                          dplyr::rename(group1 = .data$R.Condition,
                                                                        group1_present_replicate_percentage = .data$present_replicate_percentage) %>%
                                                          dplyr::select(.data$PG.ProteinGroups,
                                                                        everything()),
                                                        detection_percentage %>%
                                                          dplyr::filter(.data$R.Condition == condition_comparisons[2,i]) %>%
                                                          dplyr::rename(group2 = .data$R.Condition,
                                                                        group2_present_replicate_percentage = .data$present_replicate_percentage) %>%
                                                          dplyr::select(.data$PG.ProteinGroups,everything()),
                                                        by = "PG.ProteinGroups") %>%
                                                dplyr::mutate(slr_ratio_meta = paste(condition_comparisons[1,i],"/",condition_comparisons[2,i],sep="")) %>%
                                                dplyr::mutate(group1 = condition_comparisons[1,i],
                                                              group2= condition_comparisons[2,i])
      )

    }

    # polish detection_percentage_final NA --> 0

    detection_percentage_final <- detection_percentage_final %>%
      dplyr::mutate(group1_present_replicate_percentage = ifelse(test = is.na(.data$group1_present_replicate_percentage),yes = 0,no = .data$group1_present_replicate_percentage),
                    group2_present_replicate_percentage = ifelse(test = is.na(.data$group2_present_replicate_percentage),yes = 0,no = .data$group2_present_replicate_percentage))


    #perform signal to noise comparison ====
    message_function(text = "condition-comparison-wise signal to noise comparison ...",
                     color = "blue",
                     log_file_name = log_file_name)

    S2N_data <- SpectroPipeR_data_quant$data_input_normalized %>%
      group_by(.data$R.Condition,.data$PG.ProteinGroups) %>%
      dplyr::summarise(median_S2N = median(.data$EG.SignalToNoise, na.rm=T))

    # generate groups
    S2N_data_final <- c()
    for(i in 1:ncol(condition_comparisons)){

      S2N_data_final <- bind_rows(S2N_data_final,
                                  full_join(S2N_data %>%
                                              dplyr::filter(.data$R.Condition == condition_comparisons[1,i]) %>%
                                              dplyr::rename(group1 = .data$R.Condition, group1_median_S2N = .data$median_S2N) %>%
                                              dplyr::select(.data$PG.ProteinGroups,everything()),
                                            S2N_data %>%
                                              dplyr::filter(.data$R.Condition==condition_comparisons[2,i]) %>%
                                              dplyr::rename(group2 = .data$R.Condition, group2_median_S2N = .data$median_S2N) %>%
                                              dplyr::select(.data$PG.ProteinGroups,everything()), by = "PG.ProteinGroups") %>%
                                    dplyr::mutate(slr_ratio_meta = paste(.data$group1,"/",.data$group2,sep=""))
      )

    }

    # add stats
    S2N_data_final <- full_join(S2N_data_final,
                                results_peca,
                                by = c("group1","group2","slr_ratio_meta","PG.ProteinGroups"))
    # sorting
    S2N_data_final <- S2N_data_final %>%
      dplyr::select(.data$PG.ProteinGroups,
                    .data$group1,
                    .data$group2,
                    .data$slr_ratio_meta,
                    everything())

    # filter for 2 peptides
    S2N_data_final_n2 <- S2N_data_final %>%
      dplyr::filter(.data$n>1)

    # choose label
    S2N_data_final_n2_labeled <- S2N_data_final_n2 %>%
      dplyr::filter(.data$group2_median_S2N<1 &
                      .data$group1_median_S2N<1 &
                      .data$significant_changed!="none")


    # signal to noise plot ====
    signal_to_noise_plot <- ggplot(data = S2N_data_final_n2,
                                   mapping = aes(x = log2(.data$group2_median_S2N),
                                                 y = log2(.data$group1_median_S2N),
                                                 color = .data$significant_changed))+
      geom_point(alpha = 0.6,mapping = aes(size = n))+
      facet_wrap(~.data$slr_ratio_meta,ncol=5)+
      scale_size_continuous(range = c(0.2,8),breaks = c(2,4,6,8,10), trans = "identity")+
      scale_color_manual(values = c(down = "dodgerblue3",
                                    none = "grey",
                                    up = "orangered3"))+
      theme_light(base_size = theme_size)+
      geom_abline(intercept = 0,slope = 1)+
      geom_vline(xintercept = 0, linetype = "dashed")+
      geom_hline(yintercept = 0, linetype = "dashed")+
      labs(title = "Signal to noise comparison",
           subtitle = paste("significant changed: FC >", parameter$fold_change," & q-value < ",parameter$p_value_cutoff ,sep=""),
           color = "sign. changed (q-value)",
           x = expression(log[2]~"median group2 S/N"),
           y = expression(log[2]~"median group1 S/N"),
           caption = "labeled ProteinGroups: significant changed with a median S/N in both groups below 1")+
      geom_label_repel(data = S2N_data_final_n2_labeled,
                       mapping = aes(label = .data$PG.ProteinGroups),
                       show.legend = F,
                       min.segment.length = 0.001)+
      guides(fill = F)+
      theme(legend.position = "bottom")



    #save scatter plot
    message_function(text = "...signal to noise: save scatter plot...",
                     color = "blue",
                     log_file_name = log_file_name)

    # setup scatter plot dimensions for saving plot
    max_dimension_scatter <- ncol(condition_comparisons)/5


    if(ncol(condition_comparisons)<5){
      scatter_height <- 7+2
      scatter_width <- (max_dimension_scatter*5)*7
    }else{
      if(round(max_dimension_scatter)-max_dimension_scatter<0){
        scatter_height <- 7*(round(max_dimension_scatter)+1)+2
        scatter_width <- 35
      }else{
        scatter_height <- 7*(round(max_dimension_scatter))+2
        scatter_width <- 35
      }
    }

    # save scatter plot
    ggsave_pdf_png(filename = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Signal_to_noise__scatter_plot"),
                   plot = signal_to_noise_plot,
                   limitsize = F,
                   height = scatter_height,
                   width = scatter_width)




    message_function(text = "condition-comparison-wise comparison of peptide-int.-ratios vs. protein-int.-ratios ...",
                     color = "blue",
                     log_file_name = log_file_name)

    # get protein data
    protein_data<- SpectroPipeR_data_quant$protein_data

    # generate median protein intensity over replicates

    median_protein_data<- protein_data %>%
      group_by(.data$PG.ProteinGroups,.data$R.Condition) %>%
      dplyr::summarise(median_protein_intensity = median(.data$protein_intensity, na.rm=T))

    message_function(text = "...calculating protein ratios...",
                     color = "blue",
                     log_file_name = log_file_name)

    median_protein_data_ratios <- c()
    for(i in 1:ncol(condition_comparisons)){
      median_protein_data_ratios <- bind_rows(median_protein_data_ratios,
                                              median_protein_data %>%
                                                dplyr::filter(.data$R.Condition%in%c(condition_comparisons[1,i],condition_comparisons[2,i])) %>%
                                                group_by(.data$PG.ProteinGroups) %>%
                                                dplyr::summarise(median_protein_data_ratio = .data$median_protein_intensity[.data$R.Condition==condition_comparisons[1,i]]/
                                                                   .data$median_protein_intensity[.data$R.Condition==condition_comparisons[2,i]]) %>%
                                                dplyr::mutate(log2_median_protein_data_ratio = log2(.data$median_protein_data_ratio),
                                                              group1 = condition_comparisons[1,i],
                                                              group2 = condition_comparisons[2,i],
                                                              slr_ratio_meta = paste(condition_comparisons[1,i],"/",condition_comparisons[2,i],sep=""))
      )
    }

    message_function(text = "...combining ratio tables...",
                     color = "blue",
                     log_file_name = log_file_name)

    combined_pep_prot_ratios <- full_join(results_peca,
                                          median_protein_data_ratios,
                                          by = c("PG.ProteinGroups", "group1", "group2", "slr_ratio_meta"))

    message_function(text = "...generating ratio-ratios: protein_ratios/peptide_ratios...",
                     color = "blue",
                     log_file_name = log_file_name)

    # add log2 transformation ratios
    combined_pep_prot_ratios <- combined_pep_prot_ratios %>%
      dplyr::mutate(log2_protein_ratios__peptide_ratios = log2((2^.data$log2_median_protein_data_ratio)/(2^.data$slr)))


    # add over and under-estimation tag

    combined_pep_prot_ratios <- combined_pep_prot_ratios %>%
      group_by(.data$PG.ProteinGroups,
               .data$slr_ratio_meta) %>%
      dplyr::mutate(protein_estimation_comment = ratio_judgement_function(log2_peptide_ratio = .data$slr,
                                                                          log2_median_protein_data_ratio = .data$log2_median_protein_data_ratio,
                                                                          log2_protein_ratios__peptide_ratios = .data$log2_protein_ratios__peptide_ratios,
                                                                          log2_cut = 1)) %>%
      ungroup()

    #View(combined_pep_prot_ratios %>% select(PG.ProteinGroups,n,slr,log2_median_protein_data_ratio,log2_protein_ratios__peptide_ratios,protein_estimation_comment))

    # filter for at least 2 peptides
    message_function(text = "...filter for at least 2 peptides...",
                     color = "blue",
                     log_file_name = log_file_name)

    combined_pep_prot_ratios_n2 <- combined_pep_prot_ratios %>%
      dplyr::filter(.data$n>1)


    # add protein intensity for group1 and group2
    message_function(text = "...adding of protein intensity to table...",
                     color = "blue",
                     log_file_name = log_file_name)

    combined_pep_prot_ratios_n2 <- left_join(combined_pep_prot_ratios_n2,
                                             median_protein_data %>%
                                               dplyr::rename(group1_median_protein_intensity = .data$median_protein_intensity),
                                             by = c("group1"="R.Condition","PG.ProteinGroups"))

    combined_pep_prot_ratios_n2 <- left_join(combined_pep_prot_ratios_n2,
                                             median_protein_data %>%
                                               dplyr::rename(group2_median_protein_intensity = .data$median_protein_intensity),
                                             by = c("group2"="R.Condition","PG.ProteinGroups"))

    # add signal to noise
    message_function(text = "...add signal to noise per group...",
                     color = "blue",
                     log_file_name = log_file_name)

    combined_pep_prot_ratios_n2 <- left_join(combined_pep_prot_ratios_n2,
                                             S2N_data_final_n2 %>%
                                               dplyr::select(.data$PG.ProteinGroups,
                                                             .data$group1,
                                                             .data$group2,
                                                             .data$slr_ratio_meta,
                                                             .data$group1_median_S2N,
                                                             .data$group2_median_S2N),
                                             by = c("PG.ProteinGroups","group1","group2","slr_ratio_meta"))
    # add detection
    message_function(text = "...add detection with selected q-value cutoff with at least 2 peptides per replicate...",
                     color = "blue",
                     log_file_name = log_file_name)

    combined_pep_prot_ratios_n2 <- left_join(combined_pep_prot_ratios_n2,
                                             detection_percentage_final,
                                             by = c("PG.ProteinGroups","group1","group2","slr_ratio_meta"))




    # add direction indicator (peptide ratio and protein ratio having same direction) YES/NO
    message_function(text = "...add direction comparison for protein or peptide condition comp. ratio...",
                     color = "blue",
                     log_file_name = log_file_name)

    combined_pep_prot_ratios_n2 <- combined_pep_prot_ratios_n2 %>%
      dplyr::mutate(direction_of_protein_ratios_and_peptide_ratios = ifelse(test = .data$log2_median_protein_data_ratio>0 & .data$slr >0,
                                                                            yes = "same direction - up",
                                                                            no = ifelse(test = .data$log2_median_protein_data_ratio<0 & .data$slr <0,
                                                                                        yes = "same direction - down",
                                                                                        no = "contrary directions")))

    message_function(text = "...save table of stat. significant with poor signal to noise...",
                     color = "blue",
                     log_file_name = log_file_name)

    # filter table of stat. significant with poor signal to noise
    combined_pep_prot_ratios_n2_poor_S2N_significant <- combined_pep_prot_ratios_n2 %>%
      dplyr::filter(.data$group2_median_S2N<1 & .data$group1_median_S2N<1 & .data$significant_changed!="none")

    # save table of stat. significant with poor signal to noise
    write_csv(x = combined_pep_prot_ratios_n2_poor_S2N_significant,
              file = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Signal_to_noise__significant_proteins_poor_with_signal_to_noise.csv"),
              col_names = T)


    # generate counts > 2fold
    message_function(text = "...counting protein which having a 2fold difference...",
                     color = "blue",
                     log_file_name = log_file_name)

    #workaround for large log2_protein_ratios__peptide_ratios
    if(max(abs(combined_pep_prot_ratios_n2$slr),na.rm = T)>10){
      combined_pep_prot_ratios_counts <- combined_pep_prot_ratios_n2 %>%
        dplyr::filter(!is.na(.data$protein_estimation_comment)) %>%
        group_by(.data$slr_ratio_meta,
                 .data$protein_estimation_comment) %>%
        dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)) %>%
        dplyr::mutate(y = ifelse(test = .data$protein_estimation_comment=="protein ratio to high",
                                 yes = 0,
                                 no = ifelse(.data$protein_estimation_comment=="OK",-2,-4)),
                      x = max(combined_pep_prot_ratios_n2$slr,na.rm=T)-1) %>%
        ungroup()
    }else{
      combined_pep_prot_ratios_counts <- combined_pep_prot_ratios_n2 %>%
        dplyr::filter(!is.na(.data$protein_estimation_comment)) %>%
        group_by(.data$slr_ratio_meta,
                 .data$protein_estimation_comment) %>%
        dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)) %>%
        mutate(y = ifelse(test = .data$protein_estimation_comment=="protein ratio to high",
                          yes = -1,
                          no = ifelse(.data$protein_estimation_comment=="OK",-2,-3)),
               x = max(combined_pep_prot_ratios_n2$slr,na.rm=T)-1) %>%
        ungroup()
    }



    # get Top15 changed in each direction
    message_function(text = "...select Top15 over- or under-estimated proteins...",
                     color = "blue",
                     log_file_name = log_file_name)

    combined_pep_prot_ratios_n2_top15 <- combined_pep_prot_ratios_n2 %>%
      dplyr::filter(.data$protein_estimation_comment!="OK") %>%
      group_by(.data$slr_ratio_meta,
               .data$protein_estimation_comment) %>%
      dplyr::top_n(n = 15,wt = abs(.data$log2_protein_ratios__peptide_ratios)) %>%
      ungroup()


    # quantity_bench_scatter plot ====
    quantity_bench_scatter<- ggplot(data = combined_pep_prot_ratios_n2,
                                    mapping = aes(x = .data$slr,
                                                  y = .data$log2_median_protein_data_ratio,
                                                  color = .data$log2_protein_ratios__peptide_ratios))+
      geom_point(alpha = 0.6,
                 mapping = aes(size = .data$n))+
      scale_size_continuous(range = c(0.2,8),breaks = c(2,4,6,8,10), trans = "identity")+
      facet_wrap(~.data$slr_ratio_meta,ncol=5)+
      scale_color_gradient2(low = "dodgerblue3",
                            high = "orangered3",
                            mid = "grey",
                            midpoint = 0,
                            oob = squish,
                            limits = c(-1,1))+
      theme_light(base_size = theme_size)+
      geom_abline(intercept = 0,slope = 1)+
      geom_abline(intercept = -1,slope = 1, linetype = "dashed")+
      geom_abline(intercept = 1,slope = 1, linetype = "dashed")+
      geom_text_repel(data = combined_pep_prot_ratios_n2_top15,
                      mapping = aes(label = .data$PG.ProteinGroups))+
      labs(title = "Comparison of peptide int. ratios vs. protein int. ratios",
           subtitle = paste("pep. ratio aggr.: ", parameter$type_slr,"; protein int. alg.: ",parameter$protein_intensity_estimation,sep=""),
           color = expression(log[2]~"(protein-int. ratios/\npeptide-int. ratios)"),
           x = expression(log[2]~"peptide int. ratios"),
           y = expression(log[2]~"protein int. ratios"),
           caption = "proteins with 2 or more peptides")+
      geom_label(data = combined_pep_prot_ratios_counts,
                 mapping = aes(.data$x,
                               .data$y,
                               label = paste("count:",.data$count),
                               fill = .data$protein_estimation_comment),
                 inherit.aes = F, color = "white")+
      scale_fill_manual(values = c("protein ratio to low" = "dodgerblue3",
                                   "protein ratio to high" = "orangered3",
                                   "OK" = "grey"))+
      guides(fill = F, color = guide_colourbar(barwidth = 10, barheight = 0.5))+
      theme(legend.position = "bottom")

    #save scatter plot
    message_function(text = "...protein int. benchmark: save scatter plot...",
                     color = "blue",
                     log_file_name = log_file_name)

    # setup scatter plot dimensions for saving plot
    max_dimension_scatter <- ncol(condition_comparisons)/5



    if(ncol(condition_comparisons)<5){
      scatter_height <- 7+2
      scatter_width <- (max_dimension_scatter*5)*7
    }else{
      if(round(max_dimension_scatter)-max_dimension_scatter<0){
        scatter_height <- 7*(round(max_dimension_scatter)+1)+2
        scatter_width <- 35
      }else{
        scatter_height <- 7*(round(max_dimension_scatter))+2
        scatter_width <- 35
      }
    }

    # save scatter plot
    ggsave_pdf_png(filename = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__scatter_plot"),
                   plot = quantity_bench_scatter,
                   limitsize = F,
                   height = scatter_height,
                   width = scatter_width)




    # quantity_bench_histogram plot ====
    #set bin_width
    bin_width <- 0.1

    # estimate max counts
    max_counts_in_hist<- max(combined_pep_prot_ratios_n2 %>%
                               group_by(.data$slr_ratio_meta) %>%
                               dplyr::summarise(max_counts = max(hist(.data$log2_protein_ratios__peptide_ratios,
                                                                      breaks=seq(-10,10,by=bin_width),
                                                                      plot=FALSE)$counts)) %>%
                               dplyr::select(.data$max_counts))


    # add data for hist label
    max_ratio_ratio <- max(abs(combined_pep_prot_ratios_n2$log2_protein_ratios__peptide_ratios),na.rm=T)
    hist_x_high <- if_else(max_ratio_ratio<3,max_ratio_ratio,max_ratio_ratio-1)

    combined_pep_prot_ratios_counts <- combined_pep_prot_ratios_counts %>%
      dplyr::mutate(hist_y = ifelse(test = .data$protein_estimation_comment=="protein ratio to high",
                                    yes = max_counts_in_hist,
                                    no = ifelse(.data$protein_estimation_comment=="OK",
                                                max_counts_in_hist*0.8,
                                                max_counts_in_hist*0.6)),
                    hist_x = hist_x_high)

    # histogramm plot
    quantity_bench_hist <- ggplot(data = combined_pep_prot_ratios_n2,
                                  mapping = aes(.data$log2_protein_ratios__peptide_ratios,
                                                fill = .data$protein_estimation_comment))+
      geom_histogram(binwidth = bin_width)+
      facet_wrap(~.data$slr_ratio_meta,ncol=5)+
      geom_vline(xintercept = 0)+
      geom_vline(xintercept = 1,linetype = "dashed")+
      geom_vline(xintercept = -1,linetype = "dashed")+
      theme_light(base_size = theme_size)+
      coord_cartesian(xlim = c(-max_ratio_ratio,max_ratio_ratio))+
      labs(title = "Comparison of peptide int. ratios vs. protein int. ratios",
           fill = "protein estimation\njudgement",
           subtitle = paste("pep. ratio aggr.: ",
                            parameter$type_slr,
                            "; protein int. alg.: ",
                            parameter$protein_intensity_estimation,
                            sep=""),
           y = "count",
           x = expression(log[2]~"protein int. ratios/peptide int. ratios"),
           caption = "proteins with 2 or more peptides")+
      scale_fill_manual(values = c("protein ratio to low" = "dodgerblue3",
                                   "protein ratio to high" = "orangered3",
                                   "OK" = "grey"))+
      geom_label(data = combined_pep_prot_ratios_counts,
                 mapping = aes(x = .data$hist_x,
                               y = .data$hist_y,
                               label = paste("N =",.data$count)),
                 color = "white", size = 6,show.legend = F)

    message_function(text = "...protein int. benchmark: save histogram plot...",
                     color = "blue",
                     log_file_name = log_file_name)

    # setup histogram plot dimensions for saving plot
    max_dimension_hist <- ncol(condition_comparisons)/5

    if(ncol(condition_comparisons)<5){
      hist_height <- 5
      hist_width <- (ncol(condition_comparisons)*5)+1
    }else{
      if(round(max_dimension_scatter)-max_dimension_scatter<0){
        hist_height <- 5*(round(max_dimension_scatter)+1)
        hist_width <- 25
      }else{
        hist_height <- 5*(round(max_dimension_scatter))
        hist_width <- 25
      }
    }
    # save hist plot
    ggsave_pdf_png(filename = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__histogram_plot"),
                   plot = quantity_bench_hist,
                   limitsize = F,
                   height = hist_height,
                   width = hist_width)

    # quantity_bench_MA_like plot ====
    MA_like_plot_data <- combined_pep_prot_ratios_n2 %>%
      dplyr::select(.data$PG.ProteinGroups,
                    .data$group1,
                    .data$group2,
                    .data$n,
                    .data$slr_ratio_meta,
                    .data$log2_protein_ratios__peptide_ratios,
                    .data$group2_median_protein_intensity,
                    .data$group1_median_protein_intensity) %>%
      pivot_longer(c("group1_median_protein_intensity",
                     "group2_median_protein_intensity"),
                   names_to = "group_protein_int",
                   values_to = "protein_intensity") %>%
      mutate(group = ifelse(grepl(pattern =  "group1",x = .data$group_protein_int),
                            .data$group1,
                            .data$group2))

    # add topN highlights
    MA_like_plot_data <- left_join(MA_like_plot_data,combined_pep_prot_ratios_n2_top15 %>%
                                     dplyr::select(.data$slr_ratio_meta,
                                                   .data$PG.ProteinGroups) %>%
                                     dplyr::mutate(highlighted = TRUE),
                                   by = c("PG.ProteinGroups", "slr_ratio_meta"))



    quantity_bench_hist_MA_like_plot <- ggplot(data = MA_like_plot_data,
                                               mapping = aes(x = .data$protein_intensity,
                                                             y = .data$log2_protein_ratios__peptide_ratios,
                                                             color = .data$log2_protein_ratios__peptide_ratios))+
      geom_point(mapping = aes(size = .data$n))+
      scale_size_continuous(range = c(0.2,8),breaks = c(2,4,6,8,10), trans = "identity")+
      scale_x_log10()+
      facet_wrap(.data$slr_ratio_meta~.data$group,ncol=2)+
      scale_color_gradient2(low = "dodgerblue3",
                            high = "orangered3",
                            mid = "grey",
                            midpoint = 0,oob = squish,limits = c(-1,1))+
      theme_light(base_size = 18)+
      geom_hline(yintercept = 0)+
      geom_hline(yintercept = -1, linetype = "dashed")+
      geom_hline(yintercept = 1, linetype = "dashed")+
      geom_text_repel(data = MA_like_plot_data %>%
                        dplyr::filter(.data$highlighted==TRUE),
                      mapping = aes(label = .data$PG.ProteinGroups))+
      labs(title = "Comparison of protein int. and peptide int./protein int. ratios",
           subtitle = paste("pep. ratio aggr.: ", parameter$type_slr,"; protein int. alg.: ",parameter$protein_intensity_estimation,sep=""),
           color = expression(log[2]~"(protein-int. ratios/\npeptide-int. ratios)"),
           y = expression(log[2]~"protein-int. ratios/peptide-int. ratios"),
           x = expression("protein int."),
           caption = "proteins with 2 or more peptides")+
      guides(fill = F, color = guide_colourbar(barwidth = 10, barheight = 0.5))+
      theme(legend.position = "bottom")
    #save MA like plot
    ggsave_pdf_png(filename = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__MA_like_plot"),
                   plot = quantity_bench_hist_MA_like_plot,
                   limitsize = F,
                   height = 7*ncol(condition_comparisons),
                   width = 15)


    message_function(text = "...protein int. benchmark: save table...",
                     color = "blue",
                     log_file_name = log_file_name)

    write_csv(x = combined_pep_prot_ratios_n2,
              file = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__table.csv"),
              col_names = T)
    write_csv(x = combined_pep_prot_ratios_n2_top15,
              file = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__table_top15.csv"),
              col_names = T)

    # write_csv(x = combined_pep_prot_ratios_counts,
    #           file = paste("05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__table_counts.csv",sep=""),
    #           col_names = T)

    # barplot
    quantity_bench_count_plot <- ggplot(data = combined_pep_prot_ratios_counts,
                                        mapping = aes(.data$slr_ratio_meta,
                                                      .data$count,
                                                      fill = .data$protein_estimation_comment))+
      geom_bar(stat = "identity")+
      coord_flip()+
      theme_light(base_size = 12)+
      scale_fill_manual(values = c("protein ratio to low" = "dodgerblue3",
                                   "protein ratio to high" = "orangered3",
                                   "OK" = "grey"))+
      labs(x = "comparisons", fill = "protein estimation judgement",
           title = "Protein-intensity estimation: Counts plot", caption = "2fold cutoff")+
      theme(legend.position = "bottom")

    # save bar plot
    ggsave_pdf_png(filename = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__barplot"),
                   plot = quantity_bench_count_plot,
                   limitsize = F,
                   height = (0.8*ncol(condition_comparisons))+2,
                   width = 8)



    # protein intensity estimation bench: gradient of difference ====
    message_function(text = "...counting protein: gradient of difference...",
                     color = "blue",
                     log_file_name = log_file_name)

    protein_off_seq <- seq(1.1,4,by = 0.1)

    #setup matrix with 2 peptides only
    combined_pep_prot_ratios_2n <- combined_pep_prot_ratios %>%
      dplyr::filter(.data$n>1)


    combined_pep_prot_ratios_counts_sequence <- c()
    for(i in 1:length(protein_off_seq)){
      combined_pep_prot_ratios_counts_sequence <- bind_rows(combined_pep_prot_ratios_counts_sequence,
                                                            combined_pep_prot_ratios_n2 %>%
                                                              dplyr::filter(!is.na(.data$protein_estimation_comment)) %>%
                                                              group_by(.data$PG.ProteinGroups,
                                                                       .data$slr_ratio_meta) %>%
                                                              dplyr::mutate(protein_estimation_comment = ratio_judgement_function(log2_peptide_ratio = .data$slr,
                                                                                                                                  log2_median_protein_data_ratio = .data$log2_median_protein_data_ratio,
                                                                                                                                  log2_protein_ratios__peptide_ratios = .data$log2_protein_ratios__peptide_ratios,
                                                                                                                                  log2_cut = log2(protein_off_seq[i]))) %>%
                                                              dplyr::filter(!is.na(.data$protein_estimation_comment)) %>%
                                                              group_by(.data$slr_ratio_meta,
                                                                       .data$protein_estimation_comment) %>%
                                                              dplyr::summarise(count = n_distinct(.data$PG.ProteinGroups)) %>%
                                                              dplyr::mutate(FC_cutoff = protein_off_seq[i]) %>%
                                                              ungroup()
      )



    }

    # total protein count
    total_protein_count<- length(unique(combined_pep_prot_ratios_n2$PG.ProteinGroups))
    # area plot
    quantity_bench_area_plot <- ggplot(data = combined_pep_prot_ratios_counts_sequence,
                                       mapping = aes(x = .data$FC_cutoff,
                                                     y = .data$count,
                                                     fill=.data$protein_estimation_comment))+
      geom_area()+
      facet_wrap(~.data$slr_ratio_meta)+
      scale_y_continuous(sec.axis = sec_axis(~(./total_protein_count)*100, name = "%"))+
      theme_light(base_size = theme_size)+
      scale_fill_manual(values = c("protein ratio to low" = "dodgerblue3",
                                   "protein ratio to high" = "orangered3",
                                   "OK" = "grey"))+
      geom_vline(xintercept = c(1.1,1.5,2,4),linetype = "dashed")+
      labs(x = "FC cut", fill = "protein estimation\njudgement",
           title = "Protein-intensity estimation: area plot",
           subtitle = "sequential cutoff (Proteins with at least 2 peptides)",
           caption = "dashed lines at 1.1, 1.5, 2.0, 4.0 FC")+
      theme(legend.position = "bottom")+
      geom_label_repel(data = combined_pep_prot_ratios_counts_sequence %>%
                         dplyr::filter(.data$FC_cutoff%in%c(1.1,1.5,2,4)) %>%
                         group_by(.data$slr_ratio_meta,
                                  .data$FC_cutoff) %>%
                         mutate(y_lab = ifelse(test = .data$protein_estimation_comment == "protein ratio to low",
                                               yes = .data$count[.data$protein_estimation_comment == "protein ratio to low"],
                                               no = ifelse(test = .data$protein_estimation_comment == "protein ratio to high",
                                                           yes = sum(.data$count[.data$protein_estimation_comment %in% c("protein ratio to low","protein ratio to high")]),
                                                           no = sum(.data$count)))),
                       mapping = aes(y = .data$y_lab,
                                     label = .data$count),
                       segment.color = "black",
                       min.segment.length = 0.01,
                       color = "white",show.legend = FALSE)

    message_function(text = "...protein int. benchmark: gradient of difference area plots...",
                     color = "blue",
                     log_file_name = log_file_name)
    # save area plot
    ggsave_pdf_png(filename = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__FC_gradient_area_plot"),
                   plot = quantity_bench_area_plot,
                   limitsize = F,
                   height = scatter_height*0.7,
                   width = scatter_width*0.7)


    # add meta data table for column description ====

    meta_data_column_explained_tibble<- tibble::tibble(parameter = c("slr",
                                                                     "t",
                                                                     "score",
                                                                     "n",
                                                                     "p",
                                                                     "p.fdr",
                                                                     "PG.ProteinGroups",
                                                                     "group1",
                                                                     "group2",
                                                                     "slr_ratio_meta",
                                                                     "test",
                                                                     "type",
                                                                     "significant_changed",
                                                                     "significant_changed_raw_p",
                                                                     "significant_changed_fc",
                                                                     "significant_changed_p_value",
                                                                     "fold_change_absolute",
                                                                     "fold_change_direction",
                                                                     "fold_change",
                                                                     "median_protein_data_ratio",
                                                                     "log2_median_protein_data_ratio",
                                                                     "log2_protein_ratios__peptide_ratios",
                                                                     "protein_estimation_comment",
                                                                     "group1_median_protein_intensity",
                                                                     "group2_median_protein_intensity",
                                                                     "group1_median_S2N",
                                                                     "group2_median_S2N",
                                                                     "group1_present_replicate_percentage",
                                                                     "group2_present_replicate_percentage",
                                                                     "direction_of_protein_ratios_and_peptide_ratios",
                                                                     "effect_size_method",
                                                                     "d",
                                                                     "d_pooled_SD",
                                                                     "d_95CI_lower",
                                                                     "d_95CI_upper",
                                                                     "d_magnitute"),
                                                       description= c("signal log2-ratios on peptide basis",#slr
                                                                      "t of t-statistics on peptide basis",#t
                                                                      "score of t-statistics on peptide basis",#score
                                                                      "number of peptides",#n
                                                                      "raw p-value of statistics on peptide basis",#p
                                                                      "adjusted p-value (q-value) of statistics on peptide basis",#p.fdr
                                                                      "Protein groups",#PG.ProteinGroups
                                                                      "group1 of condition comparison",#group1
                                                                      "group2 of condition comparison",#group2
                                                                      "how the ratio is formed", # slr_ratio_meta
                                                                      "which test was used for statistics on peptide level",#test
                                                                      "which type of ratio aggregation to ProteinGroup level was used for signal log2-ratios on peptide basis",#type
                                                                      "if there is a significant change FC & q-value (cutoffs e.g.: FC = 1.5 & adjusted-p-value = 0.05)",#significant_changed
                                                                      "if there is a significant change FC & p-value(cutoffs e.g.: FC = 1.5 & p-value = 0.05)",#significant_changed_raw_p
                                                                      "fold-change cutoff used for analysis",#significant_changed_fc
                                                                      "p-value/q-value cutoff used for analysis",#significant_changed_p_value
                                                                      "ablsolute fold-change",#fold_change_absolute
                                                                      "fold-change direction",#fold_change_direction
                                                                      "fold-change",#fold_change
                                                                      "ratio of protein data (median protein intensity per condtion; ratio between 2 condtion medians)",#median_protein_data_ratio
                                                                      "log2 ratio of protein data (median protein intensity per condtion; ratio between 2 condtion medians)",#log2_median_protein_data_ratio
                                                                      "log2 RATIOprotein/RATIOpeptide per comparison and ProteinGroup",#log2_protein_ratios__peptide_ratios
                                                                      "is the protein intensity estimate ratio to high or to low in comparsion to the peptide ratio",#protein_estimation_comment
                                                                      "median protein intensity estimate of group 1",#group1_median_protein_intensity
                                                                      "median protein intensity estimate of group 2",#group2_median_protein_intensity
                                                                      "group1 median over signal-to-noise of ions",#group1_median_S2N
                                                                      "group2 median over signal-to-noise of ions",#group2_median_S2N
                                                                      "group1 present with selected q-value cutoff with at least 2 peptides per replicate: percentage of replicates per condition; NA = only 1 peptide per replicate >> may have globally 2 or more peptide",#group1_present_replicate_percentage
                                                                      "group2 present with selected q-value cutoff with at least 2 peptides per replicate: percentage of replicates per condition; NA = only 1 peptide per replicate >> may have globally 2 or more peptide",#group2_present_replicate_percentage
                                                                      "does the protein ratio and peptide ratio point into the same direction?",#direction_of_protein_ratios_and_peptide_ratios
                                                                      "effect size estimation method used (mean scaled peptides intensities are used as input)",
                                                                      "effect size estimate",
                                                                      "within-groups standard diviation",
                                                                      "the lower 95% confidence interval",
                                                                      "the upper 95% confidence interval",
                                                                      "a qualitative assessment of the magnitude of effect size (|d|<0.2 negligible, |d|<0.5 small, |d|<0.8 medium, otherwise large); Cohen 1992"
                                                       ))
    # write column annotation table for Protein benchmark table
    readr::write_csv(x = meta_data_column_explained_tibble,
                     file = paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__table_annotation.csv"),
                     col_names = T)


    #return results

    output_list <-  list(stat_results = results_peca,
                         stat_column_description = meta_data_column_explained_tibble,
                         stats_results_iBAQ_quantiles = results_peca_iBAQ,
                         stat_results_filtered = results_peca_filtered)

    message_function(text = paste0("statistics module done --> please check outputs in folder: ", out_folder,"/","06_statistics/"),
                     color = "blue",
                     log_file_name = log_file_name)

    class(output_list) <- "SpectroPipeR_statistics"
    return(output_list)
    }#end condition == NULL

  }#end condition == NULL





}


