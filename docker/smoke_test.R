suppressPackageStartupMessages({
  library(SpectroPipeR)
})

params <- list(
  output_folder = tempdir(),
  ion_q_value_cutoff = 0.001,
  id_drop_cutoff = 0.3,
  normalization_method = "median",
  normalization_factor_cutoff_outlier = 4,
  filter_oxidized_peptides = TRUE,
  protein_intensity_estimation = "MaxLFQ",
  stat_test = "modt",
  type_slr = "median",
  fold_change = 1.5,
  p_value_cutoff = 0.05,
  paired = FALSE
)

example_file_path <- system.file("extdata", "SN_test_HYE_mix_file.tsv", package = "SpectroPipeR")
if (!nzchar(example_file_path)) stop("bundled HYE example file not found in installed package")

cat("\n=== STEP 1: read_spectronaut_module ===\n")
d <- read_spectronaut_module(file = example_file_path, parameter = params, print.plot = FALSE)
stopifnot(!is.null(d))

cat("\n=== STEP 2: norm_quant_module ===\n")
q <- norm_quant_module(SpectroPipeR_data = d)
stopifnot(!is.null(q))

cat("\n=== STEP 3: MVA_module ===\n")
m <- MVA_module(SpectroPipeR_data_quant = q, HCPC_analysis = FALSE)
stopifnot(!is.null(m))

cat("\n=== STEP 4: statistics_module ===\n")
s <- statistics_module(
  SpectroPipeR_data_quant = q,
  condition_comparisons = cbind(c("HYE mix A", "HYE mix B"))
)
stopifnot(!is.null(s))

cat("\n=== ALL MODULES PASSED ===\n")
