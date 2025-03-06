## code to prepare `DATASET` dataset goes here
library(dplyr)
require(data.table)
library(MSnbase)

file <- "//fouet/spi/scidospi/06_Data/BarbierSaintHilaire_ComparativeEvaluationData_2020/DIA/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML"
rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk")

temp_dt <- rawData.onDiskMSnExp %>%
  MSnbase::filterRt(., c(537, 560)) %>%
  MSnbase::filterEmptySpectra()

fdt_sub <- fData(temp_dt) %>% as.data.table()
event_index <- fdt_sub[, which(
  (msLevel == 1) |
    (msLevel == 2 )
)]

data_example <- temp_dt[event_index]
x <- fData(data_example)

data_example %>%
  MSnbase::writeMSData(., "./inst/extdata/test_data.mzml")


# file <- "./inst/extdata/test_data.mzML"
# file <- get_test_data()
data_example <- get_test_peaks()

usethis::use_data(data_example, overwrite = TRUE)
