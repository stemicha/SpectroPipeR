# usethis::use_build_ignore(files = "tests/util__generate_test_file_XIC_test.R")


library(tidyverse)
library(RSQLite)

ions <- c("_YEEMVENMK_.2","_YGPSLMPGGSEEAWPHIK_.3")


XIC_folder <- "../input_data/20240517_060926_HYE_mix_demo/20240516_154329_HYE_mix_demo_XIC-DBs/"

files <- list.files(path = XIC_folder,full.names = T)
#create tmp dir
dir.create("tmp_XIC")
#copy files
file.copy(from = files,to = "tmp_XIC")

files_new <- list.files(path = "tmp_XIC",full.names = T)

# get report
tmp_report <- read_delim(file = "../input_data/20240517_060926_HYE_mix_demo/HYE_mix_demo_Report_SpectroPipeR (Normal).tsv",delim = "\t")

#fgid selection
fgid_selection <- unlist(tmp_report %>%
                           filter(FG.Id %in% ions) %>%
                           distinct(FG.XICDBID))

for(i in 1:length(files_new)){
  helfRlein::statusbar(run = i, max.run = length(files_new), width = 60L, info = "processing...")
  # Step 2: Connect to the SQLite database
  con <- dbConnect(RSQLite::SQLite(), dbname = files_new[i])

  # Step 3: Load the data from the SQLite file
  data <- dbGetQuery(con, "SELECT * FROM IonTraces")

  # Step 4: Filter the data based on your criteria
  # For example, let's filter rows where column 'a' is greater than 5
  filtered_data <- data %>%
    filter(FGID%in%fgid_selection)

  # Step 5: Delete the original table
  dbExecute(con, "DROP TABLE IF EXISTS IonTraces")

  # Step 6: Write the filtered data back to the SQLite file
  dbWriteTable(con, "IonTraces", filtered_data)

  # Step 7: reduce file size
  dbExecute(con, "VACUUM")
  # Step 8: Close the connection
  dbDisconnect(con)
  file.copy(from = files_new[i],to = paste0("inst/extdata/HYE_demo_data/XIC_DBs/file",i,".xic.db"))
}
unlink("tmp_XIC",recursive = T)

write_delim(x = tmp_report %>%
              filter(FG.Id %in% ions),
            file = "inst/extdata/HYE_demo_data/HYE_demo_data_Report_SpectroPipeR.tsv",delim = "\t")
