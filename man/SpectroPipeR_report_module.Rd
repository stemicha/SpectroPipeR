% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_module.R
\name{SpectroPipeR_report_module}
\alias{SpectroPipeR_report_module}
\title{SpectroPipeR: report module}
\usage{
SpectroPipeR_report_module(
  SpectroPipeR_data = NULL,
  SpectroPipeR_data_quant = NULL,
  SpectroPipeR_data_stats = NULL,
  open_rendered_report = FALSE
)
}
\arguments{
\item{SpectroPipeR_data}{(mandatory) is the SpectroPipeR_data list object from read_spectronaut_module() object e.g. \code{SpectroPipeR_data} please see example below}

\item{SpectroPipeR_data_quant}{(mandatory) is the SpectroPipeR_data_quant list object from norm_quant_module() object e.g. \code{SpectroPipeR_data_quant} please see example below}

\item{SpectroPipeR_data_stats}{(optional) is the SpectroPipeR_data_quant list object from statistics_module() object e.g. \code{SpectroPipeR_data_stats} please see example below; if you have only 1 replicate data set a statistical analysis is not possible so leave this parameter to NULL a report without statitsics will be generated}

\item{open_rendered_report}{logical - if rendered report should be opened (TRUE or FALSE); default = FALSE}
}
\value{
generates and exports an interactive standalone html report to the output folder
}
\description{
Function for generating an interactive analysis report (HTML file). This is the final part of the SpectroPipeR pipeline where the analysis modules output serves as input.
}
\examples{
\donttest{
#load library
library(SpectroPipeR)

# use default parameters list
params <- list(output_folder = "../SpectroPipeR_test_folder")

# example input file
example_file_path <- system.file("extdata",
                                "SN_test_HYE_mix_file.tsv",
                                 package="SpectroPipeR")

# step 1: load Spectronaut data module
SpectroPipeR_data <- read_spectronaut_module(file = example_file_path,
                                            parameter = params,
                                            print.plot = FALSE)

# step 2: normalize & quantification module
SpectroPipeR_data_quant <- norm_quant_module(SpectroPipeR_data = SpectroPipeR_data)

# step 3: MVA module
SpectroPipeR_MVA <- MVA_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
          HCPC_analysis = FALSE)

# step 4: statistics module
SpectroPipeR_data_stats <- statistics_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
                                            condition_comparisons = cbind(c("HYE mix A",
                                                                            "HYE mix B")))

# step 5: report module
SpectroPipeR_report_module(SpectroPipeR_data = SpectroPipeR_data,
                          SpectroPipeR_data_quant = SpectroPipeR_data_quant,
                          SpectroPipeR_data_stats = SpectroPipeR_data_stats)
}
}
