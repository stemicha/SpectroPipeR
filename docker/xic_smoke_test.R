suppressPackageStartupMessages({
  library(SpectroPipeR)
})

# Output dir override (env var) so we can run this on host or in container
out_dir <- Sys.getenv("XIC_OUT_DIR", unset = file.path(tempdir(), "xic_smoke_out"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Spectronaut_report_path <- system.file(
  "extdata/HYE_demo_data", "HYE_demo_data_Report_SpectroPipeR.tsv",
  package = "SpectroPipeR"
)
Spectronaut_xicDB_path  <- system.file(
  "extdata/HYE_demo_data/XIC_DBs",
  package = "SpectroPipeR"
)
stopifnot(nzchar(Spectronaut_report_path), file.exists(Spectronaut_report_path))
stopifnot(nzchar(Spectronaut_xicDB_path),  dir.exists(Spectronaut_xicDB_path))
stopifnot(length(list.files(Spectronaut_xicDB_path, pattern = "\\.db$")) >= 1)

cat("\n=== XIC_plot_module ===\n")
cat("report : ", Spectronaut_report_path, "\n")
cat("xic db : ", Spectronaut_xicDB_path,  "\n")
cat("out    : ", out_dir, "\n\n")

XIC_plot_module(
  Spectronaut_report_path = Spectronaut_report_path,
  Spectronaut_xicDB_path  = Spectronaut_xicDB_path,
  protein_groups          = c("P29311", "P38720"),
  output_path             = out_dir,
  export_csv_files        = TRUE,
  number_of_cores         = 2
)

pdfs <- list.files(out_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
csvs <- list.files(out_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

cat("\nproduced PDFs:\n"); for (f in pdfs) cat(" -", f, "(", file.info(f)$size, "bytes )\n")
cat("\nproduced CSVs:\n"); for (f in csvs) cat(" -", f, "(", file.info(f)$size, "bytes )\n")

if (length(pdfs) == 0) stop("XIC_plot_module produced no PDF files")
if (length(csvs) == 0) stop("XIC_plot_module produced no CSV files (export_csv_files = TRUE)")
if (any(file.info(pdfs)$size == 0)) stop("at least one PDF is 0 bytes")

cat("\n=== XIC PLOT TEST PASSED ===\n")
cat("PDFs:", length(pdfs), " CSVs:", length(csvs), "\n")
