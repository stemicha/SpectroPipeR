# usethis::use_build_ignore(files = "tests/util__generate_test_file.R")


library(tidyverse)

test_file<- read_delim("../input_data/20240517_060926_HYE_mix_demo/HYE_mix_demo_Report_SpectroPipeR (Normal).tsv",
           delim = "\t",
           col_names = T,
           col_types = cols(EG.Qvalue = col_character(),
                            PG.IBAQ = col_character())
)



test_file_filter<- test_file %>%
  distinct(R.FileName,
           R.Condition,
           R.Replicate,
           `R.Raw File Name`,
           `R.Instrument Name`,
           `R.MS1 Mass Analyzer`,
           `R.MS2 Mass Analyzer`,
           `R.Run Date`,
           PG.ProteinGroups,
           PG.Organisms,
           PG.IBAQ,
           PEP.StrippedSequence,
           EG.ModifiedPeptide,
           PEP.NrOfMissedCleavages,
           EG.UserGroup,
           EG.Qvalue,
           EG.PEP,
           EG.Cscore,
           EG.NormalizationFactor,
           `EG.TotalQuantity (Settings)`,
           EG.SignalToNoise,
           EG.Identified,
           EG.ApexRT,
           EG.IntCorrScore,
           EG.DatapointsPerPeak,
           `EG.DatapointsPerPeak (MS1)`,
           FG.Charge,
           FG.Id,
           FG.XICDBID,
           FG.LabeledSequence,
           FG.ShapeQualityScore,
           FG.MS1Quantity,
           FG.MS2Quantity,
           FG.MS1RawQuantity,
           FG.MS2RawQuantity)

#use only 500 proteins per species / demo
proteins <- unlist(test_file_filter %>%
  distinct(PG.Organisms,PG.ProteinGroups) %>%
  group_by(PG.Organisms) %>%
  slice_sample(n = 30) %>%
  distinct(PG.ProteinGroups))

#write output to extdata folder
write_delim(test_file_filter %>%
              filter(PG.ProteinGroups %in% proteins) %>%
              mutate(R.Condition = str_replace_all(R.Condition,"A_manual","HYE mix A"),
                     R.Condition = str_replace_all(R.Condition,"B_manual","HYE mix B")),
            file = "inst/extdata/SN_test_HYE_mix_file.tsv",delim = "\t")

