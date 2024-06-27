
#test ID list
test_that("read_spectronaut_module works", {
  # parameter list
  params <- list(output_folder = tempdir(), # character: output folder
                 ion_q_value_cutoff = 0.001,              # numeric: Biognosys default Q-value = 0.01 = 1% error rate
                 id_drop_cutoff = 0.3,                    # numeric: value 0-1 (1 = 100%)
                 normalization_method = "median",         # character: "median" or "spectronaut; default = auto-detection
                 normalization_factor_cutoff_outlier = 4, # numeric: abs. norm. factor cut-off -> outlier definition
                 filter_oxidized_peptides = TRUE,         # logical: TRUE or FALSE
                 protein_intensity_estimation = "MaxLFQ", # character: Hi3 or MaxLFQ
                 stat_test = "modt",                      # character: statistical test to be used "rots" or "modt" or "t"
                 type_slr = "median",                     # character: "median" or "tukey"
                 fold_change = 1.5,                       # numeric: fold-change cutoff
                 p_value_cutoff = 0.05,                   # numeric: p-value cutoff
                 paired = FALSE                           # logical: should a paired statistics be used?
  )

  # example input file
  example_file_path <- system.file("extdata", "SN_test_HYE_mix_file.tsv", package="SpectroPipeR")

  # step 1: load Spectronaut data module
  SpectroPipeR_data <- read_spectronaut_module(file = example_file_path,
                                         parameter = params,
                                         print.plot = FALSE)

  # step 2: normalize & quantification module
  SpectroPipeR_data_quant <- norm_quant_module(SpectroPipeR_data = SpectroPipeR_data)

  # step 3: MVA module
  SpectroPipeR_MVA <- MVA_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
             HCPC_analysis = F)

  # step 4: statistics module
  SpectroPipeR_data_stats <- statistics_module(SpectroPipeR_data_quant = SpectroPipeR_data_quant,
                                         condition_comparisons = cbind(c("HYE mix A","HYE mix B")))


  # remove temp. output folder
  unlink(params$output_folder)


})
