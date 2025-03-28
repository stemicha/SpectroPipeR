% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpectroPipeR.R
\name{SpectroPipeR}
\alias{SpectroPipeR}
\title{SpectroPipeR: complete Spectronaut analysis with one function}
\usage{
SpectroPipeR(
  file = "",
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
)
}
\arguments{
\item{file}{location (path) of Spectronaut output report;
you should use the \code{Spectronaut_export_scheme()} function for getting a SpectroPipeR report scheme encompassing all mandatory columns}

\item{parameter}{\strong{mandatory parameter list element}

\emph{table of list elements:}\tabular{ll}{
   \if{html}{\out{<u>}} \strong{parameter} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{description} \if{html}{\out{</u>}} \cr
   output_folder \tab \strong{mandatory !!!} - \emph{character} - output folder path (abs.) \cr
   ion_q_value_cutoff \tab \strong{default = 0.01} - \emph{numeric} - Q-value used in Spectronaut analysis: Biognosys \cr
    \tab default is 0.01 = 1\% error rate \cr
   id_drop_cutoff \tab \strong{default = 0.3} - \emph{numeric} - value between 0-1 (1 = 100\%); xx percent lower \cr
    \tab than median of ion ID rate => outlier \cr
   normalization_method \tab \strong{default = "median"} - \emph{character} - "median" or Spectronaut - auto-detection \cr
    \tab is per default ON, meaning if normalization was performed in Spectronaut \cr
    \tab this will be detected and preferred over parameter setting here; \cr
    \tab median normalization is the fallback option \cr
   normalization_factor_cutoff_outlier \tab \strong{default = 4} - \emph{numeric} - median off from global median \cr
    \tab (4 means abs. 4fold off) \cr
   filter_oxidized_peptides \tab \strong{default = TRUE} \emph{logical} - if oxidized peptides should be removed before \cr
    \tab peptide quantification \cr
   protein_intensity_estimation \tab \strong{default = "MaxLFQ"} - \emph{character} - Hi3 = Hi3 protein intensity estimation, \cr
    \tab MaxLFQ = MaxLFQ protein intensity estimation \cr
   stat_test \tab \strong{default = "rots"} - \emph{character} - choose statistical test: "rots" = reproducibility \cr
    \tab optimized test statistics, "modt" = moderate t-test (lmfit, eBayes), \cr
    \tab "t" = t-test \cr
   type_slr \tab \strong{default = "median"} - \emph{character} - choose ratio aggregation method: \cr
    \tab "median" or "tukey" is used when calculating protein values \cr
   fold_change \tab \strong{default = 1.5} - \emph{numeric} - fold-change used as cutoff e.g. 1.5 \cr
   p_value_cutoff \tab \strong{default = 0.05} - \emph{numeric} - p-value used as cutoff e.g. 0.05 \cr
   paired \tab \strong{default = FALSE} - \emph{logical} - Should paired statistics be applied? \cr
}
\tabular{l}{
   \if{html}{\out{<u>}}example parameters list (default)\if{html}{\out{</u>}}: \cr
   params <- list(output_folder = "../Spectronaut_example", \cr
   ion_q_value_cutoff = 0.01, \cr
   id_drop_cutoff = 0.3, \cr
   normalization_method = "median", \cr
   normalization_factor_cutoff_outlier = 4, \cr
   filter_oxidized_peptides = T, \cr
   protein_intensity_estimation = "MaxLFQ", \cr
   stat_test = "rots", \cr
   type_slr = "median", \cr
   fold_change = 1.5, \cr
   p_value_cutoff = 0.05, \cr
   paired = FALSE \cr
   ) \cr
}}

\item{max_chars_file_name_capping}{integer, (default = 25) number of max characters used for raw file name presentation; must be adjusted if function}

\item{ID_condition_filtering}{TRUE or FALSE if a condition-wise filtering should be performed}

\item{ID_condition_filtering_percent}{(numerical value ranging from 0 - 1, default = 0.5) define the proportion for the condition-wise ID filtering}

\item{batch_adjusting}{logical - if batch adjusting with ComBat (sva package) should be performed; default = FALSE}

\item{sample__batch_meta_data_file}{character - sample batch file; tab-delimited txt-file, containing "R.FileName" column e.g. sample__batch_meta_data_file = "Sample_MetaData_Batches.txt"

\emph{example table for batch meta data:}\tabular{ll}{
   \if{html}{\out{<u>}} \strong{R.FileName} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{digest_batch} \if{html}{\out{</u>}} \cr
   20230403_TIMSTOF_1_S1-B11_1_6690 \tab 1 \cr
   20230403_TIMSTOF_2_S1-G11_1_6695 \tab 1 \cr
   20230403_TIMSTOF_3_S1-E7_1_6661 \tab 1 \cr
   20230403_TIMSTOF_4_S1-A9_1_6673 \tab 2 \cr
   20230403_TIMSTOF_7_S1-D8_1_6668 \tab 2 \cr
   20230403_TIMSTOF_9_S1-D3_1_6627 \tab 2 \cr
}


\strong{A good starting point for the generation of the table is the '*_ConditionSetup.tsv' in your Spectronaut Pipeline Report export folder}}

\item{batch_adjusting_column}{character - column name in sample__batch_meta_data_file, which should be used for assigning the samples to batches}

\item{number_of_cores_adjusting}{numeric - number of processor cores used for batch or covariate adjustment}

\item{covariate_adjusting_formula}{character - provide a formula passed to lm() for covariate adjustment e.g. "log10_peptide_intensity ~ log10(CRP)+log10(age)+as.factor(sex)"; you may also use ns() function e.g. "log10_peptide_intensity ~ ns(age, df=3)"}

\item{covariate_adjusting_meta_data_file}{covariate meta csv file, containing "R.FileName; age; sex;..."; you may find a start file in the 02_ID_rate folder > file_list.csv column e.g. covariate_adjusting_meta_data_file = "covariate_MetaData_file.csv"
\emph{example table for covariate meta data:}\tabular{llll}{
   \if{html}{\out{<u>}} \strong{R.FileName} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{R.Condition} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{sex} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{CRP} \if{html}{\out{</u>}} \cr
   20230403_TIMSTOF_1_S1-B11_1_6690 \tab heathy \tab 1 \tab 3.8 \cr
   20230403_TIMSTOF_2_S1-G11_1_6695 \tab heathy \tab 2 \tab 5.1 \cr
   20230403_TIMSTOF_3_S1-E7_1_6661 \tab heathy \tab 1 \tab 1.2 \cr
   20230403_TIMSTOF_4_S1-A9_1_6673 \tab cancer \tab 1 \tab 50.2 \cr
   20230403_TIMSTOF_7_S1-D8_1_6668 \tab cancer \tab 2 \tab 30.8 \cr
   20230403_TIMSTOF_9_S1-D3_1_6627 \tab cancer \tab 2 \tab 64.1 \cr
}}

\item{skipping_MaxLFQ_median_norm}{logical - if median normalization after MaxLFQ calculation should be skipped; default = FALSE; applied only if MaxLFQ protein estimation is selected}

\item{HCPC_analysis}{boolean; should a HCPC be performed or not}

\item{costum_colors}{if you would like to use your own colors for condition coloring please provide a named color vector (e.g. c(condition1 = "black", condition2 = "grey")); names should have the same naming and length like the conditions set in Spectronaut}

\item{condition_comparisons}{condition comparisons for pairwise- comparison; e.g. condition_comparisons <- cbind(c("condition1","control"),c("condition3","control") )}

\item{number_of_cores_statistics}{number of processor cores to be used for the calculations default = 2;

\code{parallel::detectCores()-2} for faster processing (will detect the number of cores in the system and use nearly all cores)}

\item{build_HTML_report}{boolean; if a HTML report of the analysis should be generated or not}

\item{report_copy}{if TRUE --> copy Spectronaut input report to SpectroPipeR project folder 01_input_data}
}
\value{
SpectroPipeR list object containing tables and plots of the analysis in addition to the automatically saved tables and plots.
For the description of the generated figures and tables please read the manual & vignettes

\strong{The SpectroPipeR list element contains:}
\itemize{
\item SpectroPipeR_data
\item SpectroPipeR_data_quant
\item SpectroPipeR_data_MVA
\item SpectroPipeR_data_stats
}

\strong{SpectroPipeR_data:}\tabular{ll}{
   \if{html}{\out{<u>}} \strong{list element} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{description} \if{html}{\out{</u>}} \cr
   spectronaut_output \tab \emph{tibble:} Spectronaut report tibble provided for the analysis \cr
   SDRF_file \tab \emph{tibble:} intermediate SDRF table of the analysis \cr
   summary_distinct \tab \emph{tibble:} distinct ion, modified peptide, stripped peptides and \cr
    \tab protein group count per file filtered by provided Q-value \cr
   raw_file_names \tab \emph{tibble:} R.FileNames capped and uncapped version together with \cr
    \tab R.Condition and R.Replicate \cr
   ion_id_median \tab \emph{numerical value:} median of ion intensity \cr
   ion_id_cutoff \tab \emph{numerical value:} ion ID count threshold to classify sample as outlier \cr
   PG_2_peptides_ID_raw \tab \emph{tibble:} with protein groups with at least 2 peptides with peptide \cr
    \tab and replicate count \cr
   summary_distinct_outlier \tab \emph{tibble:} if outlier are detected they are listed in this tibble \cr
   ID_rate_plot \tab \emph{ggplot2 plot:} ID rate plot \cr
   ID_rate_plot_filter \tab \emph{ggplot2 plot:} ion ID rate plot with ion ID cutoff line \cr
   sample_length \tab \emph{numberical value:} number of samples in the provided Spectronaut report \cr
   parameter \tab \emph{list:} parameters provided by the user \cr
   time_stamp_log_file \tab \emph{string:} time stamp of the log file (format: \%Y_\%m_\%d__\%H_\%M) \cr
   log_file_name \tab \emph{string:} analysis log file name \cr
}


\strong{SpectroPipeR_data_quant:}\tabular{ll}{
   \if{html}{\out{<u>}} \strong{list element} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{description} \if{html}{\out{</u>}} \cr
   data_input_normalized \tab \emph{tibble:} Spectronaut report tibble provided for the analysis \cr
   MedianNormalizationFactor \tab \emph{tibble:} ion normalization factor table \cr
   MedianNormalizationFactor_outlier \tab \emph{tibble:} table containing detected norm. outliers on ion level \cr
   NormFactor_plot \tab \emph{ggplot2 plot:} norm. factor plots \cr
   iBAQ_intensities \tab \emph{tibble:} table containing the iBAQ int. \cr
   iBAQ_intensities_summary \tab \emph{tibble:} table containing the per condition summarized iBAQ int. \cr
   protein_data \tab \emph{tibble:} protein intensity table (e.g. Hi3 or MaxLFQ int.) \cr
   PG_2_peptides_ID_raw \tab \emph{tibble:} with protein groups with at least 2 peptides with peptide \cr
    \tab and replicate count \cr
   protein_data_normalization_factor \tab \emph{tibble:} normalization factor table for protein int. data \cr
   peptide_intensity_filtered_2pep_hi3 \tab \emph{tibble:} if Hi3 protein int. was selected a table containing \cr
    \tab the peptides and intensities used for Hi3 protein intensity calculation \cr
   peptide_intensity \tab \emph{tibble:} peptides intensity table based on norm. ion intensity \cr
   parameter \tab \emph{list:} parameters provided by the user updated with the \cr
    \tab norm_quant_module() parameters \cr
   CV_cumulative_frequency \tab \emph{tibble:} cumulative frequency table on peptide \cr
    \tab and protein intensity level \cr
   sample_length \tab \emph{numberical value:} number of samples in the provided Spectronaut report \cr
}


\strong{SpectroPipeR_data_MVA:}\tabular{ll}{
   \if{html}{\out{<u>}} \strong{list element} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{description} \if{html}{\out{</u>}} \cr
   PCA_peptide_intensity \tab \emph{PCA list element:} PCA list element of peptide int. \cr
   PCA_protein_intensity \tab \emph{PCA list element:} PCA list element of protein int. \cr
   UMAP_protein_intensity \tab \emph{umap element:} UMAP element of protein int. \cr
   peptide_intensity_correlation \tab \emph{matrix:} Spearman correlation scores of peptide int. \cr
   protein_intensity_correlation \tab \emph{matrix:} Spearman correlation scores of protein int. \cr
}


\strong{SpectroPipeR_data_stats:}\tabular{ll}{
   \if{html}{\out{<u>}} \strong{list element} \if{html}{\out{</u>}} \tab \if{html}{\out{<u>}} \strong{description} \if{html}{\out{</u>}} \cr
   stat_results \tab \emph{tibble:} statistical analysis results table \cr
   stat_column_description \tab \emph{tibble:} statistical analysis results table column description \cr
   stats_results_iBAQ_quantiles \tab \emph{tibble:} statistical analysis results table containing \cr
    \tab the iBAQ quantilies (Q1-Q10) of the protein per group for a better \cr
    \tab ratio judgement \cr
   stat_results_filtered \tab \emph{tibble:} filtered (user defined FC and adj. p-value) statistical \cr
    \tab analysis results table \cr
}
}
\description{
Function for performing the whole SpectroPipeR analysis workflow.
}
\details{
\strong{batch adjustment}

Batch effects refer to systematic differences between batches (groups) of samples in high-throughput experiments.
These differences can arise due to various factors, such as batch variations in sample preparation, handling, processing
procedures and measurement orders. Batch effects can obscure the true biological signal and lead to incorrect conclusions
if not properly accounted for.
In the SpectroPipeR pipeline, the ComBat tool was employed to adjust for batch effects in the datasets where the batch
covariate was known. ComBat utilizes the methodology described in \href{https://pubmed.ncbi.nlm.nih.gov/16632515/}{Johnson et al. 2007}.
It uses an empirical Bayes (EB) framework for adjusting data for batch effects that is robust to outliers in small sample sizes
and performs comparable to existing methods for large samples. \href{https://pubmed.ncbi.nlm.nih.gov/16632515/}{Johnson et al. 2007:}
This method incorporates systematic batch biases common across genes in making adjustments, assuming that phenomena resulting in
batch effects often affect many genes in similar ways (i.e. increased expression, higher variability, etc). Specifically, the
the L/S model parameters are estimated that represent the batch effects by pooling information across peptides in each
batch to shrink the batch effect parameter estimates toward the overall mean of the batch effect estimates (across genes).
These EB estimates are then used to adjust the data for batch effects, providing more robust adjustments for the batch effect on each peptide.
In SpectroPipeR a parametric ComBAT emperical Bayes adjustment is implemented by utilizing the sva-package.

\strong{covariate adjustment}

If a covariate adjustment of peptide intensity data was performed using the users input formula,
a linear mixed model (LMM) was calculated based on that formula per peptide and the outcoming
residuals were added to the mean peptide intensity over the samples. This means that the adjusted
peptide intensities retain their intensity level (low intense peptides keep their low intensity and
high intense ions keep their higher intensity).
}
\examples{
\donttest{
# load library
library(SpectroPipeR)

# use default parameters list
params <- list(output_folder = "../SpectroPipeR_test_folder")

# example input file
example_file_path <- system.file("extdata",
                                "SN_test_HYE_mix_file.tsv",
                                package="SpectroPipeR")
# perform the analysis
SpectroPipeR_analysis <- SpectroPipeR(file = example_file_path,
                                     parameter = params,
                                     condition_comparisons = cbind(c("HYE mix A","HYE mix B"))
                                     )
}
}
