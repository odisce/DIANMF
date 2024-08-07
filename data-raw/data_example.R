## code to prepare `data_example` dataset goes here

library(dplyr)
require(data.table)
library(MSnbase)

file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");

temp_dt <- rawData.onDiskMSnExp %>%
  MSnbase::filterRt(., c(62, 78)) %>%
  MSnbase::filterMz(., c(95,115)) 

fdt_sub <- fData(temp_dt) %>% as.data.table()

event_index <- fdt_sub[, which(
  (msLevel == 1) |
    (msLevel == 2 & precursorMZ <= 160)
)]

temp_dt[event_index] %>%
  MSnbase::writeMSData(., "../test_data.mzml")

file <- "../test_data.mzml";
data_example <- MSnbase::readMSData(file, mode = "onDisk");

usethis::use_data(data_example, overwrite = TRUE)
