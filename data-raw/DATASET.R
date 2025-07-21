## code to prepare `DATASET` dataset goes here
library(dplyr)
require(data.table)
library(MSnbase)
library(tools)

input_dir <- "../mzml/"
input_files <- list.files(input_dir, pattern = ".mzml", full.names = TRUE, ignore.case = TRUE)

for(i in c(17,18)){
  file <- input_files[i]
  rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk")
  temp_dt <- rawData.onDiskMSnExp %>%
    MSnbase::filterRt(., c(280, 320)) %>%
    MSnbase::filterMz(., c(250, 320)) %>%
    MSnbase::filterEmptySpectra()

  fdt_sub <- fData(temp_dt) %>% as.data.table()
  event_index <- fdt_sub[, which(
    (msLevel == 1) |
      (msLevel == 2 )
  )]

  data_example <- temp_dt[event_index]
  data_example %>%
    MSnbase::writeMSData(., paste0("./inst/extdata/test_data_",file_path_sans_ext(basename(file)),".mzml"))
}

# create the msexp object on both mzml files
# rt_range = c(280, 320)
# data_example <- get_test_peaks()
# 
# usethis::use_data(data_example, overwrite = TRUE)
 