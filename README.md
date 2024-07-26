
# SpectroPipeR <img src="vignettes/figures/SpectroPipeR_hexbin_logo.png" align="right" width="200" />

a streamlining post Spectronaut™ DIA-MS data analysis R package

The manual can be found under https://stemicha.github.io/SpectroPipeR/

<!-- badges: start -->
<!-- badges: end -->

## Introduction

During proteome studies, researchers frequently face various challenges, some of which can be mitigated while others cannot. A prevalent issue is the bottleneck in downstream data analysis, which arises due to a limited number of bioinformaticians, rapid generation of raw data, and variations in data analysis methods or workflows.

To tackle this problem SpectroPipeR was developed.  This pipeline is designed to simplify data analysis tasks, significantly reduce the workload for scientists, be easily expandable, user-friendly even for those with minimal bioinformatics knowledge, generate standardized analysis, outputs, and reports for each project, and produce publication-ready tables, figures, and reports.

SpectroPipeR comprises a set of R functions that facilitate a comprehensive, fully automated, and standardized data analysis of Spectronaut DIA-MS data. This includes ID rate summary, ON/OFF analysis, normalization, batch or covariate adjustment, iBAQ and maxLFQ quantification, multivariate analysis, peptide-centric statistical analysis (ROPECA or modified t-test), and interactive HTML report generation. The output is presented through a variety of clear graphs and tables in a well-structured folder system. The comprehensive standalone HTML report is extremely useful for existing Electronic Laboratory Notebooks (ELN) or Laboratory Information Management Systems (LIMS) to quickly obtain a project-specific overview.

SpectroPipeR consists of a global parameter setting and four analysis modules and one reporting module that are executed sequentially. This modular approach allows flexibility where specific analyses like ID- and intensity plots can be run independently or as part of the complete pipeline.

After each module execution, dynamic console feedback is provided and written to a log file to help identify errors early on. Upon completion of all modules, a comprehensive set of tables and plots categorized in different folders is generated to summarize the project from various perspectives.

SpectroPipeR includes also a module for XIC plotting, capable of generating protein-specific XIC plots for each ion associated with the protein. This is complemented by a range of significant metrics designed to aid in evaluating the accuracy of both identification and quantification.

## Installation

You can install the development version of SpectroPipeR like so:

- install R using https://cran.r-project.org.
- <u>optional:</u> install RStudio IDE using https://posit.co/download/rstudio-desktop/
- inside R execute the following code:

``` r
#install devtools
install.packages("devtools")
# install SpectroPipeR from github
devtools::install_github("stemicha/SpectroPipeR")

# quit and restart R or restart R session in Rstudio
```

## Requirements

For the interactive html report feature SpectroPipeR needs **Quarto CLI**.
Quarto is an open-source scientific and technical publishing system.
You can install the Quarto CLI using [Quarto get started installation](https://quarto.org/docs/get-started/).

## Spectronaut report requirements

SpectroPipeR requires certain columns from the Spectronaut output report that are not included by default.

_The following steps are advised:_

1. **download and installation of SpectroPipeR (if not already be done)**
2. **load the SpectroPipeR package and utilize the `Spectronaut_export_scheme()` function to create the necessary Spectronaut report scheme (SpectroPipeR_report.rs) in the output folder provided.**
``` r
# output_location: path to output folder for the SpectroPipeR_report.rs Spectronaut report scheme
Spectronaut_export_scheme(output_location = "../SpectroPipeR_test_folder")
```
3. **import the generated SpectroPipeR report scheme (SpectroPipeR_report.rs) into Spectronaut**
4. **conduct an analysis of your raw mass spectrometry data in Spectronaut and define conditions during analysis setup process**
5. **produce the output report (*.tsv) using the imported SpectroPipeR report scheme**
6. **open R and use/edit code below to perform your analysis**

``` r
# load library
library(SpectroPipeR)

# edit strings inside quotes

# use/edit default parameters list // define output folder (mandatory !)
params <- list(output_folder = "[your path to output folder]")

# add your path to your Spectronaut report (*.tsv)
file_path <- "[your path to Spectronaut report (*.tsv)]"

# add your condition comparisons as named in Spectronaut analysis setup
condition_comp = cbind(c("[condition_2]","[condition_1]"),
                       c("[condition_3]","[condition_1]")
                        )

# launch analysis
SpectroPipeR_analysis <- SpectroPipeR(file = file_path,
                                      parameter = params,
                                      condition_comparisons = condition_comp
                                      )

```


Spectronaut output report should contain the following columns to work in SpectroPipeR (these are included when using the scheme generated by `Spectronaut_export_scheme()`):



<u>mandatory Spectronaut report columns:</u>

R.FileName, R.Condition, R.Replicate, R.Instrument Name, R.Raw File Name, R.MS1 Mass Analyzer, R.MS2 Mass Analyzer, R.Run Date, PG.ProteinGroups, PG.Organisms, PG.IBAQ, PEP.StrippedSequence, EG.ModifiedPeptide, PEP.NrOfMissedCleavages, EG.UserGroup, EG.Qvalue, EG.PEP, EG.Cscore, EG.NormalizationFactor, EG.TotalQuantity (Settings), EG.SignalToNoise, EG.Identified, EG.ApexRT, EG.IntCorrScore, EG.DatapointsPerPeak, EG.DatapointsPerPeak (MS1), FG.Charge, FG.Id, FG.XICDBID, FG.LabeledSequence, FG.ShapeQualityScore, FG.MS1Quantity, FG.MS2Quantity, FG.MS1RawQuantity, FG.MS2RawQuantity

EG.TotalQuantity (Settings) is used for the quantification. 
Per default MS2 level should be selected in the quantification setting in Spectronaut™.

## code example

### all-in-one function example:

``` r
# load library
library(SpectroPipeR)

# use default parameters list
params <- list(output_folder = "../SpectroPipeR_test_folder")

# example input file // or path to your Spectronaut report (*.tsv)
example_file_path <- system.file("extdata", "SN_test_HYE_mix_file.tsv", package="SpectroPipeR")

# launch analysis
SpectroPipeR_analysis <- SpectroPipeR(file = example_file_path,
                                      parameter = params,
                                      condition_comparisons = cbind(c("HYE mix A",
                                                                      "HYE mix B"))
                                      )

```




### single function execution example:

``` r
# load library
library(SpectroPipeR)


# use default parameters list
params <- list(output_folder = "../SpectroPipeR_test_folder")

# example input file // or path to your Spectronaut report (*.tsv)
example_file_path <- system.file("extdata",
                                 "SN_test_HYE_mix_file.tsv",
                                 package="SpectroPipeR")

# step 1: load Spectronaut data module
SpectroPipeR_data <- read_spectronaut_module(file = example_file_path,
                                             parameter = params)

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

```

If all the necessary packages are correctly installed, you can also execute SpectroPipeR analysis from the terminal using bash. Sometimes, extensive analyses require more computational power and memory, so it might be beneficial to run them on a server.

Here is a basic example that you can use as a starting point, which you can save as "SpectroPipeR_terminal.sh"

``` bash
#!/bin/bash

# Command line arguments
input_file=$1
parameter=$2
condition_comparisons=$3


# usage
# bash SpectroPipeR_terminal.sh "SN_test_HYE_mix_file.tsv" "output_folder=SpectroPipeR_test_folder;stat_test=modt;ion_q_value_cutoff=0.001" "HYE mix A,HYE mix B;HYE mix B,HYE mix A"

# R script
Rscript -e "
library(SpectroPipeR)

print(\"$parameter\")
print(\"$input_file\")
print(\"$condition_comparisons\")


# split parameters
params <- strsplit(\"$parameter\",split = \";\")
params <- strsplit(unlist(params),\"=\")
names(params) <- lapply(params, function(x) x[1])
params <- lapply(params, function(x) x[2])

# convert parameters class
convert_params <- function(params) {
  numeric_params = c(\"ion_q_value_cutoff\",
                     \"id_drop_cutoff\",
                     \"normalization_factor_cutoff_outlier\",
                     \"fold_change\",
                     \"p_value_cutoff\")
  logical_params = c(\"filter_oxidized_peptides\",
                     \"paired\")
  # iterate over each parameter
  for (param in names(params)) {
    # check if parameter is in numeric_params
    if (param %in% numeric_params) {
      # convert to numeric
      params[[param]] <- as.numeric(params[[param]])
    }
    # check if parameter is in logical_params
    else if (param %in% logical_params) {
      # convert to logical
      params[[param]] <- as.logical(params[[param]])
    }
  }
  return(params)
}

params <- convert_params(params)

# split condition comparisons
cond_comp <- strsplit(\"$condition_comparisons\",split = \";\")
cond_comp <- sapply(unlist(cond_comp), function(x) strsplit(x,split = \",\"))
cond_comp <- do.call(cbind,cond_comp)

class(params)
print(params)

# perform the analysis
SpectroPipeR_analysis <- SpectroPipeR(file = \"$input_file\",
                                     parameter = params,
                                     condition_comparisons = cond_comp
                                     )
"


```

example terminal usage:

bash SpectroPipeR_terminal.sh [input_file] [named parameters separated by ;] [condition comparisons separated by ;]

``` terminal
bash SpectroPipeR_terminal.sh "SN_test_HYE_mix_file.tsv" "output_folder=SpectroPipeR_test_folder;stat_test=modt;ion_q_value_cutoff=0.001" "HYE mix A,HYE mix B;HYE mix B,HYE mix A"
```
