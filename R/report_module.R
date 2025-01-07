#' SpectroPipeR: report module
#' @description
#' Function for generating an interactive analysis report (HTML file). This is the final part of the SpectroPipeR pipeline where the analysis modules output serves as input.
#'
#' @param SpectroPipeR_data  (mandatory) is the SpectroPipeR_data list object from read_spectronaut_module() object e.g. `SpectroPipeR_data` please see example below
#' @param SpectroPipeR_data_quant (mandatory) is the SpectroPipeR_data_quant list object from norm_quant_module() object e.g. `SpectroPipeR_data_quant` please see example below
#' @param SpectroPipeR_data_stats (optional) is the SpectroPipeR_data_quant list object from statistics_module() object e.g. `SpectroPipeR_data_stats` please see example below; if you have only 1 replicate data set a statistical analysis is not possible so leave this parameter to NULL a report without statitsics will be generated
#' @param open_rendered_report logical - if rendered report should be opened (TRUE or FALSE); default = FALSE
#'
#' @returns generates and exports an interactive standalone html report to the output folder
#' @export
#'
#' @import tidyverse
#' @importFrom dplyr filter
#' @importFrom stats median sd lm as.formula quantile
#' @importFrom DT formatStyle datatable
#' @importFrom kableExtra kable_styling
#' @import readxl
#' @import htmltools
#' @import htmlwidgets
#' @import quarto
#' @importFrom methods is
#' @import rmarkdown
#' @import knitr
#' @import readxl
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
#'                                  package="SpectroPipeR")
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
#'SpectroPipeR_MVA <- MVA_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
#'           HCPC_analysis = FALSE)
#'
#'# step 4: statistics module
#'SpectroPipeR_data_stats <- statistics_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
#'                                             condition_comparisons = cbind(c("HYE mix A",
#'                                                                             "HYE mix B")))
#'
#'# step 5: report module
#'SpectroPipeR_report_module(SpectroPipeR_data = SpectroPipeR_data,
#'                           SpectroPipeR_data_quant = SpectroPipeR_data_quant,
#'                           SpectroPipeR_data_stats = SpectroPipeR_data_stats)
#'}


SpectroPipeR_report_module <- function(SpectroPipeR_data = NULL,
                                 SpectroPipeR_data_quant = NULL,
                                 SpectroPipeR_data_stats = NULL,
                                 open_rendered_report = FALSE){


  #output folder
  out_folder <- SpectroPipeR_data$parameter$output_folder
  #time stamp of log file
  time_stamp_log_file <- SpectroPipeR_data$time_stamp_log_file

  #sample length
  sample_length <- SpectroPipeR_data$sample_length

  #check parameters
  #parameter_check(parameter)
  #TODO: if paramters added like batch adjusting !!! switch to norm quant parameters
  parameter <- SpectroPipeR_data_quant$parameter
  log_file_name <- SpectroPipeR_data_quant$parameter$log_file_name

  message_function(text = "",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "# REPORT MODULE",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "#*****************************************",
                   color = "white",
                   log_file_name = log_file_name)
  message_function(text = "",
                   color = "white",
                   log_file_name = log_file_name)

  message_function(text = "generating methods part ...",
                   color = "blue",
                   log_file_name = log_file_name)


# convert parameters to tibble --------------------------------------------
parameter_tibble <- as_tibble(parameter) %>%
  pivot_longer(cols = names(parameter),
               names_to = "parameter",
               values_to = "value",
               values_transform = list(value = as.character))


# project sample overview -------------------------------------------------
project_samples <- SpectroPipeR_data$spectronaut_output %>%
    dplyr::distinct(.data$R.FileName,
                    .data$R.Condition,
                    .data$R.Replicate)

project_samples_aggregated <- project_samples %>%
  dplyr::group_by(.data$R.Condition) %>%
  dplyr::summarise(sample_count = n_distinct(.data$R.FileName))

# protein/ion counts ----------------------------------------------------------

protein_counts <- readr::read_csv(paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count__strippedPEP_without_Qvalue_cut.csv"),show_col_types = FALSE)
protein_counts_condense<- protein_counts %>%
  dplyr::distinct(.data$protein_count,
                  .data$number_of_peptides)

protein_count_plot_all <- paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count__strippedPEP_without_Qvalue_cut.png")
protein_count_plot_filtered <- paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/protein_count__strippedPEP.png")
ion_id_rate <- paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/ID_counts_plot_ion_filter.png")

# ON/OFF analysis ---------------------------------------------------------
PG_2_peptides_ID_wide_binary <- readr::read_csv(file = paste0(out_folder,"/","02_ID_rate/",
                                                 sample_length,
                                                 "_sample_analysis/Detected_ProteinGroups__UpSetR__plot__binary_coded.csv"),
                                                show_col_types = FALSE
                                                )
UpSetR_PG_2_peptides_ID_wide_binary <- paste0(out_folder,"/","02_ID_rate/",sample_length,"_sample_analysis/Detected_ProteinGroups__UpSetR__plot.png")


# normalization of data ---------------------------------------------------
normalization_boxplot <- list.files(paste0(out_folder,"/","03_normalization/",sample_length,"_sample_analysis"),
                              pattern = "normalization_boxplot.png",full.names = T)

normalization_factors <- read_csv(paste0(out_folder,"/","03_normalization/",sample_length,"_sample_analysis/Median_normalization_factors.csv"),show_col_types = FALSE)

protein_normalization_plot <- list.files(paste0(out_folder,"/","03_normalization/",sample_length,"_sample_analysis"),
                                         pattern = "normalization_factor_plot.png",full.names = T)

protein_normalization_boxplot <- list.files(paste0(out_folder,"/","03_normalization/",sample_length,"_sample_analysis"),
                                         pattern = "normalization_factor_BOXplot.png",full.names = T)


# cumulative CV -----------------------------------------------------------

cumsum_CV_20percent <- read_csv(paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/CV_20percent_cumulative_frequency.csv"),show_col_types = FALSE)
cumsum_CV_plot <-  list.files(paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis"),
                              pattern = "CV_cumulative_frequency_plot.png",full.names = T)



# was batch adjusting performed -------------------------------------------

batch_adjusting <- SpectroPipeR_data_quant$parameter$batch_adjusting


# get MaxLFQ distribution plot --------------------------------------------

MaxLFQ_intensity_distribution_plot <- paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/MaxLFQ_protein_intensity_boxplot.png")

# PCA analysis & correlation ------------------------------------------------

pca_1_2_plot_cond <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_protein_level_conditions_marked.png")
pca_1_2_plot_rep <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_protein_level_replicates_marked.png")
pca_1_2_plot_measurement_order <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_protein_level_measurement_order.png")
PCA_3d_plot_protein <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/3D_PCA_plot__protein_level.html")
pca_overview_plot <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plots.png")
pca_1_5_plot_pep <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__peptide_level.png")
pca_1_5_plot_prot <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__protein_level.png")
pca_1_5_plot_pep_measurement_order <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__peptide_level_measurement_order.png")
pca_1_5_plot_prot_measurement_order <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_plot_1st_to_5th_dimension__protein_level_measurement_order.png")
pca_contrib <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/PCA_individuals_contribution_plot.png")
correlation_plot <- paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis/correlation_plots.png")


# UMAP lookup -------------------------------------------------------------

umap_plot <- list.files(paste0(out_folder,"/","04_multivariate_analysis/",sample_length,"_sample_analysis"),
                        pattern = "UMAP_plot_protein_level_conditions_marked.png",full.names = T)

# CV ----------------------------------------------------------------------

cv_plot <- paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/CV_vs_intensity_plot.png")


# statistics --------------------------------------------------------------
#if statistics is provided
if(!is.null(SpectroPipeR_data_stats)){
  #Protein_intensity_benchmark__barplot
  Protein_intensity_benchmark__barplot <- paste0(out_folder,"/","05_processed_data/",sample_length,"_sample_analysis/Protein_intensity_benchmark__barplot.png")

  #load stats
  statistics <- readr::read_csv(file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis.csv"),show_col_types = FALSE)
  statistics_2_pep <- statistics %>%
    dplyr::filter(.data$n>=2) %>%
    dplyr::select(.data$PG.ProteinGroups,
                  .data$slr_ratio_meta,
                  .data$slr,
                  .data$p.fdr,
                  .data$p,
                  .data$n,
                  .data$fold_change)
  statistics_comparison <- statistics %>%
    dplyr::distinct(.data$group1,
                    .data$group2,
                    .data$slr_ratio_meta)

  #load stats iBAQ
  statistics_iBAQ <- readr::read_csv(file = paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis/statistical_analysis_iBAQ_added.csv"),show_col_types = FALSE)
  statistics_iBAQ_2_pep <- statistics_iBAQ %>%
    dplyr::filter(.data$n>=2) %>%
    dplyr::select(.data$PG.ProteinGroups,
                  .data$slr_ratio_meta,
                  .data$slr,
                  .data$p.fdr,
                  .data$p,
                  .data$n,
                  .data$fold_change,
                  .data$iBAQ_quantile_comp)
  #load quantile color gradient
  quantile_color_gradient <- readxl::read_excel(path = system.file("extdata", "iBAQ_dynamic_range_2D_Color_gradient.xlsx", package="SpectroPipeR"))
  quantile_color_gradient_tidy <- quantile_color_gradient %>%
    pivot_longer(cols = colnames(quantile_color_gradient)[-1],
                 names_to = "group1",
                 values_to = "colors")
  colnames(quantile_color_gradient_tidy)[1] <- "group2"
  quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>%
    dplyr::select(.data$group1,
                  .data$group2,
                  .data$colors)
  quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>%
    dplyr::mutate(group1 = str_replace_all(.data$group1,"condition1_",""))
  quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>%
    dplyr::mutate(group2 = str_replace_all(.data$group2,"condition2_",""))
  quantile_color_gradient_tidy <- quantile_color_gradient_tidy %>%
    dplyr::rowwise() %>%
    dplyr::mutate(iBAQ_quantile_comp = paste0(.data$group1,"/",.data$group2))

  #generate quantile color vector
  quantile_color_gradient_vector <- quantile_color_gradient_tidy$colors
  names(quantile_color_gradient_vector) <- quantile_color_gradient_tidy$iBAQ_quantile_comp

  quantile_color_gradient_plot <- list.files(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"),
                                             pattern = "iBAQ_quantile_ratio_comparison_legend.png",full.names = T)
  # volcano plots -----------------------------------------------------------

  volcano_plots_pFDR <- list.files(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"),
                                   pattern = "volcano_plots_adjusted_p_value",full.names = T)
  volcano_plots_pFDR <- volcano_plots_pFDR[grep(pattern = "png",x = volcano_plots_pFDR)]

  volcano_plots_effSize_pFDR <- list.files(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"),
                                           pattern = "volcano_plots_effect_size_shape_adjusted_p_value",full.names = T)
  volcano_plots_effSize_pFDR <- volcano_plots_effSize_pFDR[grep(pattern = "png",x = volcano_plots_effSize_pFDR)]

  volcano_plots_pRAW <- list.files(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"),
                                   pattern = "volcano_plots_raw_p_value",full.names = T)
  volcano_plots_pRAW <- volcano_plots_pRAW[grep(pattern = "png",x = volcano_plots_pRAW)]


  # cutoff test -------------------------------------------------------------

  cutoff_simple_test_plot <- list.files(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis"),
                                        pattern = "cutoff_simple_test",full.names = T)
  cutoff_simple_test_plot <- cutoff_simple_test_plot[grep(pattern = "png",x = cutoff_simple_test_plot)]

  cutoff_test_plot <- list.files(paste0(out_folder,"/","06_statistics/",sample_length,"_sample_analysis",sep = ""),
                                 pattern = "cutoff_test",full.names = T)
  cutoff_test_plot <- cutoff_test_plot[grep(pattern = "png",x = cutoff_test_plot)]

}else{
  statistics <-  NA
  statistics_2_pep <-  NA
  statistics_iBAQ_2_pep <-  NA
  Protein_intensity_benchmark__barplot <-  NA
  quantile_color_gradient_plot <-  NA
  quantile_color_gradient_tidy <-  NA
  statistics_comparison <-  NA
  volcano_plots_pFDR <-  NA
  volcano_plots_effSize_pFDR <-  NA
  volcano_plots_pRAW <-  NA
  cutoff_simple_test_plot <-  NA
  cutoff_test_plot <-  NA
}

# params setting for quarto rendering -------------------------------------

params_quarto <- list(
               output_folder = out_folder,
               sample_length = sample_length,
               parameter_tibble = parameter_tibble,
               batch_adjusting = batch_adjusting,
               protein_intensity_estimation = parameter$protein_intensity_estimation,
               filter_oxidized_peptides = parameter$filter_oxidized_peptides,
               stat_test = parameter$stat_test,
               project_samples_aggregated = project_samples_aggregated,
               project_samples = project_samples,
               protein_counts = protein_counts,
               protein_counts_condense = protein_counts_condense,
               protein_count_plot_all = protein_count_plot_all,
               protein_count_plot_filtered = protein_count_plot_filtered,
               MaxLFQ_intensity_distribution_plot = MaxLFQ_intensity_distribution_plot,
               cumsum_CV_20percent = cumsum_CV_20percent,
               cumsum_CV_plot = cumsum_CV_plot,
               ion_id_rate = ion_id_rate,
               PG_2_peptides_ID_wide_binary = PG_2_peptides_ID_wide_binary,
               UpSetR_PG_2_peptides_ID_wide_binary = UpSetR_PG_2_peptides_ID_wide_binary,
               normalization_factors = normalization_factors,
               normalization_boxplot = normalization_boxplot,
               protein_normalization_plot = protein_normalization_plot,
               protein_normalization_boxplot = protein_normalization_boxplot,
               pca_1_2_plot_cond = pca_1_2_plot_cond,
               pca_1_2_plot_rep = pca_1_2_plot_rep,
               pca_1_2_plot_measurement_order = pca_1_2_plot_measurement_order,
               PCA_3d_plot_protein = PCA_3d_plot_protein,
               pca_overview_plot = pca_overview_plot,
               pca_1_5_plot_pep = pca_1_5_plot_pep,
               pca_1_5_plot_prot = pca_1_5_plot_prot,
               pca_1_5_plot_pep_measurement_order = pca_1_5_plot_pep_measurement_order,
               pca_1_5_plot_prot_measurement_order = pca_1_5_plot_prot_measurement_order,
               pca_contrib = pca_contrib,
               umap_plot = umap_plot,
               correlation_plot = correlation_plot,
               cv_plot = cv_plot,
               statistics = statistics,
               statistics_2_pep = statistics_2_pep,
               statistics_iBAQ_2_pep = statistics_iBAQ_2_pep,
               Protein_intensity_benchmark__barplot = Protein_intensity_benchmark__barplot,
               quantile_color_gradient_plot = quantile_color_gradient_plot,
               quantile_color_gradient_tidy = quantile_color_gradient_tidy,
               statistics_comparison = statistics_comparison,
               volcano_plots_pFDR = volcano_plots_pFDR,
               volcano_plots_effSize_pFDR = volcano_plots_effSize_pFDR,
               volcano_plots_pRAW = volcano_plots_pRAW,
               cutoff_simple_test_plot = cutoff_simple_test_plot,
               cutoff_test_plot = cutoff_test_plot
               )


# render quarto html report -----------------------------------------------
message_function(text = "render HTML report ... this might take a while",
                 color = "blue",
                 log_file_name = log_file_name)

#get current working directory
wd_current <- getwd()

#copy quarto dependent files
#copy quarto document with and without stats
if(!is.null(SpectroPipeR_data_stats)){
  file.copy(from = system.file("extdata", "DIA_MS_analysis_report_Master.qmd", package="SpectroPipeR"),
            to = paste0(out_folder,"/","DIA_MS_analysis_report_Master.qmd"),
            overwrite = T)
}else{
  file.copy(from = system.file("extdata", "DIA_MS_analysis_report_Master_wo_stats.qmd", package="SpectroPipeR"),
            to = paste0(out_folder,"/","DIA_MS_analysis_report_Master.qmd"),
            overwrite = T)
}
file.copy(from = system.file("extdata", "SpectroPipeR_hexbin_logo.png", package="SpectroPipeR"),
          to = paste0(out_folder,"/","SpectroPipeR_hexbin_logo.png"),
          overwrite = T)
file.copy(from = system.file("extdata", "styles.css", package="SpectroPipeR"),
          to = paste0(out_folder,"/","styles.css"),
          overwrite = T)

#set working directory to render report
setwd(out_folder)

quarto::quarto_render(input = "DIA_MS_analysis_report_Master.qmd",
                      #execute_dir = out_folder,
                      execute_params = params_quarto,
                      output_file = paste0(time_stamp_log_file,
                                           "_SpectroPipeR_report.html")
)
#set wd back to initial
setwd(wd_current)

#remove temp files
file.remove(paste0(out_folder,"/","DIA_MS_analysis_report_Master.qmd"))
file.remove(paste0(out_folder,"/","SpectroPipeR_hexbin_logo.png"))
file.remove(paste0(out_folder,"/","styles.css"))


# open report
if(open_rendered_report==TRUE){
  browseURL(url = paste0(out_folder,"/",time_stamp_log_file,
                         "_SpectroPipeR_report.html"))
}

#set working directory back to used
#setwd(wd_current)

message_function(text = "render HTML report ... DONE!",
                 color = "blue",
                 log_file_name = log_file_name)
# move and delete temp file
#file.copy(from = "single_module_R_scripts/DIA_MS_analysis_report_Master.html",to = "DIA_MS_analysis_report.html",overwrite = T)
#file.remove(from = "single_module_R_scripts/DIA_MS_analysis_report_Master.html")




} # end of function

