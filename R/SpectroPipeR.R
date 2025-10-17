#' SpectroPipeR: complete Spectronaut analysis with one function
#' @description
#' Function for performing the whole SpectroPipeR analysis workflow.
#'
#' @param file location (path) of Spectronaut output report;
#'   you should use the `Spectronaut_export_scheme()` function for getting a SpectroPipeR report scheme encompassing all mandatory columns
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
#' | protein_intensity_estimation | **default = "MaxLFQ"** - _character_ |
#' |                              |  Hi3 = Hi3 protein intensity estimation |
#' |                              |  MaxLFQ = MaxLFQ protein intensity estimation |
#' |                              |  directLFQ = directLFQ protein intensity estimation |
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
#' @param report_copy if TRUE --> copy Spectronaut input report to SpectroPipeR project folder 01_input_data
#' @param ID_condition_filtering TRUE or FALSE if a condition-wise filtering should be performed
#' @param ID_condition_filtering_percent (numerical value ranging from 0 - 1, default = 0.5) define the proportion for the condition-wise ID filtering
#' @param max_chars_file_name_capping integer, (default = 25) number of max characters used for raw file name presentation; must be adjusted if function
#' @param batch_adjusting logical - if batch adjusting with ComBat (sva package) should be performed; default = FALSE
#' @param batch_adjusting_column character - column name in sample__batch_meta_data_file, which should be used for assigning the samples to batches
#' @param sample__batch_meta_data_file character - sample batch file; tab-delimited txt-file, containing "R.FileName" column e.g. sample__batch_meta_data_file = "Sample_MetaData_Batches.txt"
#'
#' _example table for batch meta data:_
#' | <u> __R.FileName__ </u> | <u> __digest_batch__ </u> |
#' |:------------------------|:------------------------|
#' |20230403_TIMSTOF_1_S1-B11_1_6690|1|
#' |20230403_TIMSTOF_2_S1-G11_1_6695|1|
#' |20230403_TIMSTOF_3_S1-E7_1_6661|1|
#' |20230403_TIMSTOF_4_S1-A9_1_6673|2|
#' |20230403_TIMSTOF_7_S1-D8_1_6668|2|
#' |20230403_TIMSTOF_9_S1-D3_1_6627|2|
#'
#' **A good starting point for the generation of the table is the '*_ConditionSetup.tsv' in your Spectronaut Pipeline Report export folder**
#'
#' @param skipping_MaxLFQ_median_norm logical - if median normalization after MaxLFQ calculation should be skipped; default = FALSE; applied only if MaxLFQ protein estimation is selected
#' @param number_of_cores_adjusting numeric - number of processor cores used for batch or covariate adjustment
#' @param covariate_adjusting_formula character - provide a formula passed to lm() for covariate adjustment e.g. "log10_peptide_intensity ~ log10(CRP)+log10(age)+as.factor(sex)"; you may also use ns() function e.g. "log10_peptide_intensity ~ ns(age, df=3)"
#' @param covariate_adjusting_meta_data_file covariate meta csv file, containing "R.FileName; age; sex;..."; you may find a start file in the 02_ID_rate folder > file_list.csv column e.g. covariate_adjusting_meta_data_file = "covariate_MetaData_file.csv"
#' _example table for covariate meta data:_
#' | <u> __R.FileName__ </u> | <u> __R.Condition__ </u>| <u> __sex__ </u>| <u> __CRP__ </u>|
#' |:------------------------|:----------------|:----------------|:----------------|
#' |20230403_TIMSTOF_1_S1-B11_1_6690|heathy|1|3.8|
#' |20230403_TIMSTOF_2_S1-G11_1_6695|heathy|2|5.1|
#' |20230403_TIMSTOF_3_S1-E7_1_6661|heathy|1|1.2|
#' |20230403_TIMSTOF_4_S1-A9_1_6673|cancer|1|50.2|
#' |20230403_TIMSTOF_7_S1-D8_1_6668|cancer|2|30.8|
#' |20230403_TIMSTOF_9_S1-D3_1_6627|cancer|2|64.1|
#'
#' @param HCPC_analysis boolean; should a HCPC be performed or not
#' @param costum_colors if you would like to use your own colors for condition coloring please provide a named color vector (e.g. c(condition1 = "black", condition2 = "grey")); names should have the same naming and length like the conditions set in Spectronaut
#' @param condition_comparisons condition comparisons for pairwise- comparison; e.g. condition_comparisons <- cbind(c("condition1","control"),c("condition3","control") )
#' @param number_of_cores_statistics number of processor cores to be used for the calculations default = 2;
#'
#' `parallel::detectCores()-2` for faster processing (will detect the number of cores in the system and use nearly all cores)
#'
#' @param build_HTML_report boolean; if a HTML report of the analysis should be generated or not
#'
#' @return SpectroPipeR list object containing tables and plots of the analysis in addition to the automatically saved tables and plots.
#'  For the description of the generated figures and tables please read the manual & vignettes
#'
#'  **The SpectroPipeR list element contains:**
#'  - SpectroPipeR_data
#'  - SpectroPipeR_data_quant
#'  - SpectroPipeR_data_MVA
#'  - SpectroPipeR_data_stats
#'
#' **SpectroPipeR_data:**
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
#' **SpectroPipeR_data_quant:**
#' | <u> __list element__ </u> | <u> __description__ </u>              |
#' |:--------------------------|:--------------------------------------|
#' | data_input_normalized     | *tibble:* Spectronaut report tibble provided for the analysis |
#' | MedianNormalizationFactor | *tibble:* ion normalization factor table |
#' | MedianNormalizationFactor_outlier | *tibble:* table containing detected norm. outliers on ion level |
#' | NormFactor_plot            | *ggplot2 plot:* norm. factor plots |
#' | iBAQ_intensities           | *tibble:* table containing the iBAQ int. |
#' | iBAQ_intensities_summary   | *tibble:* table containing the per condition summarized iBAQ int. |
#' | protein_data      | *tibble:* protein intensity table (e.g. Hi3 or MaxLFQ, directLFQ)|
#' | PG_2_peptides_ID_raw      | *tibble:* with protein groups with at least 2 peptides with peptide
#' |                           | and replicate count |
#' | protein_data_normalization_factor| *tibble:* normalization factor table for protein int. data |
#' | peptide_intensity_filtered_2pep_hi3| *tibble:* if Hi3 protein int. was selected a table containing |
#' |                            | the peptides and intensities used for Hi3 protein intensity calculation |
#' | peptide_intensity             | *tibble:* peptides intensity table based on norm. ion intensity |
#' | parameter                 | *list:* parameters provided by the user updated with the |
#' |                           | norm_quant_module() parameters|
#' | CV_cumulative_frequency   | *tibble:* cumulative frequency table on peptide |
#' |                           | and protein intensity level |
#' | sample_length             | *numberical value:* number of samples in the provided Spectronaut report |
#'
#' **SpectroPipeR_data_MVA:**
#' | <u> __list element__ </u> | <u> __description__ </u>              |
#' |:--------------------------|:--------------------------------------|
#' | PCA_peptide_intensity          | *PCA list element:* PCA list element of peptide int. |
#' | PCA_protein_intensity          | *PCA list element:* PCA list element of protein int. |
#' | UMAP_protein_intensity         | *umap element:* UMAP element of protein int. |
#' | peptide_intensity_correlation  | *matrix:* Spearman correlation scores of peptide int. |
#' | protein_intensity_correlation  | *matrix:* Spearman correlation scores of protein int. |
#'
#' **SpectroPipeR_data_stats:**
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
#' @details
#' **batch adjustment**
#'
#' Batch effects refer to systematic differences between batches (groups) of samples in high-throughput experiments.
#' These differences can arise due to various factors, such as batch variations in sample preparation, handling, processing
#' procedures and measurement orders. Batch effects can obscure the true biological signal and lead to incorrect conclusions
#' if not properly accounted for.
#' In the SpectroPipeR pipeline, the ComBat tool was employed to adjust for batch effects in the datasets where the batch
#' covariate was known. ComBat utilizes the methodology described in [Johnson et al. 2007](https://pubmed.ncbi.nlm.nih.gov/16632515/).
#' It uses an empirical Bayes (EB) framework for adjusting data for batch effects that is robust to outliers in small sample sizes
#' and performs comparable to existing methods for large samples. [Johnson et al. 2007:](https://pubmed.ncbi.nlm.nih.gov/16632515/)
#' This method incorporates systematic batch biases common across genes in making adjustments, assuming that phenomena resulting in
#' batch effects often affect many genes in similar ways (i.e. increased expression, higher variability, etc). Specifically, the
#' the L/S model parameters are estimated that represent the batch effects by pooling information across peptides in each
#' batch to shrink the batch effect parameter estimates toward the overall mean of the batch effect estimates (across genes).
#' These EB estimates are then used to adjust the data for batch effects, providing more robust adjustments for the batch effect on each peptide.
#' In SpectroPipeR a parametric ComBAT emperical Bayes adjustment is implemented by utilizing the sva-package.
#'
#' **covariate adjustment**
#'
#' If a covariate adjustment of peptide intensity data was performed using the users input formula,
#' a linear mixed model (LMM) was calculated based on that formula per peptide and the outcoming
#' residuals were added to the mean peptide intensity over the samples. This means that the adjusted
#' peptide intensities retain their intensity level (low intense peptides keep their low intensity and
#' high intense ions keep their higher intensity).
#'
#'
#' @export
#' @md
#'
#' @examples
#' \donttest{
#'# load library
#'library(SpectroPipeR)
#'
#'# use default parameters list
#'params <- list(output_folder = "../SpectroPipeR_test_folder")
#'
#'# example input file
#'example_file_path <- system.file("extdata",
#'                                 "SN_test_HYE_mix_file.tsv",
#'                                 package="SpectroPipeR")
#'# perform the analysis
#'SpectroPipeR_analysis <- SpectroPipeR(file = example_file_path,
#'                                      parameter = params,
#'                                      condition_comparisons = cbind(c("HYE mix A","HYE mix B"))
#'                                      )
#'}

SpectroPipeR <- function( file = "",
                          parameter = list(),
                          max_chars_file_name_capping = 25,
                          ID_condition_filtering = FALSE,
                          ID_condition_filtering_percent = 0.5,
                          batch_adjusting = FALSE,
                          sample__batch_meta_data_file = NULL,
                          batch_adjusting_column = "",
                          number_of_cores_adjusting = parallel::detectCores() - 2,
                          covariate_adjusting_formula = "",
                          covariate_adjusting_meta_data_file = "",
                          skipping_MaxLFQ_median_norm = FALSE,
                          HCPC_analysis = FALSE,
                          costum_colors = NULL,
                          condition_comparisons = NULL,
                          number_of_cores_statistics = 2,
                          build_HTML_report = TRUE,
                          report_copy = FALSE
                      ){
  #switch OFF warnings, font issue; In grid.Call(C_textBounds, as.graphicsAnnot(x$label),  ... :font metrics unknown for character 0xa
  options( warn = -1 )

  # step 1: load Spectronaut data module
  SpectroPipeR_data <- read_spectronaut_module(file = file,
                                               parameter = parameter,
                                               print.plot = FALSE,
                                               report_copy = report_copy,
                                               max_chars_file_name_capping = max_chars_file_name_capping,
                                               ID_condition_filtering = ID_condition_filtering,
                                               ID_condition_filtering_percent = ID_condition_filtering_percent)

  # early error messages for condition_comparisons parameter error
  if(!is.null(condition_comparisons)){
    log_file_name <- SpectroPipeR_data$log_file_name
    # test if all conditions are in input
    if(sum(unique(SpectroPipeR_data$spectronaut_output$R.Condition)%in%
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
  }


  # step 2: normalize & quantification module
  SpectroPipeR_data_quant <- norm_quant_module(SpectroPipeR_data = SpectroPipeR_data,
                                               print.plot = FALSE,
                                               batch_adjusting = batch_adjusting,
                                               sample__batch_meta_data_file = sample__batch_meta_data_file,
                                               batch_adjusting_column = batch_adjusting_column,
                                               number_of_cores_adjusting = number_of_cores_adjusting,
                                               skipping_MaxLFQ_median_norm = skipping_MaxLFQ_median_norm,
                                               covariate_adjusting_formula = covariate_adjusting_formula,
                                               covariate_adjusting_meta_data_file = covariate_adjusting_meta_data_file,
                                               costum_colors = costum_colors)

  # step 3: MVA module
  SpectroPipeR_data_MVA <- MVA_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
             costum_colors = costum_colors,
             HCPC_analysis = FALSE)

  # step 4: statistics module
  SpectroPipeR_data_stats <- statistics_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
                                               number_of_cores = number_of_cores_statistics,
                                               condition_comparisons = condition_comparisons)

  # step 5: report module
  if(build_HTML_report == TRUE){
    SpectroPipeR_report_module(SpectroPipeR_data = SpectroPipeR_data,
                               SpectroPipeR_data_quant = SpectroPipeR_data_quant,
                               SpectroPipeR_data_stats = SpectroPipeR_data_stats)
  }

  #return SpectroPipeR analysis as list element
  output_list <- list(SpectroPipeR_data = SpectroPipeR_data,
                 SpectroPipeR_data_quant = SpectroPipeR_data_quant,
                 SpectroPipeR_data_MVA = SpectroPipeR_data_MVA,
                 SpectroPipeR_data_stats = SpectroPipeR_data_stats)
  class(output_list) <- "SpectroPipeR"
  return(output_list)

  #switch warnings ON again
  options( warn = 0 )

}
