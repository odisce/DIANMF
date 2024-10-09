## code to prepare `data_example` data-set goes here

library(dplyr)
require(data.table)
library(MSnbase)

file <- "C:/Users/DK273056/Documents/1workflow/DIA_NMF_workflow/mzML/20170803_FS-DIA-E2-10ng-rep3_pos_51.mzML";
rawData.onDiskMSnExp <- MSnbase::readMSData(file, mode = "onDisk");

# for atropine
# temp_dt <- rawData.onDiskMSnExp %>%
#   MSnbase::filterRt(., c(62, 78)) %>%
#   MSnbase::filterMz(., c(95,115)) 
#   # %>%  MSnbase::filterEmptySpectra() # delete the empty spectra
# 
# fdt_sub <- fData(temp_dt) %>% as.data.table()
# event_index <- fdt_sub[, which(
#   (msLevel == 1) |
#     (msLevel == 2 & precursorMZ <= 160)
# )]

# for dextromethorphan
temp_dt <- rawData.onDiskMSnExp %>%
  MSnbase::filterRt(., c(425, 440)) %>%
  MSnbase::filterMz(., c(100, 290)) %>%
  MSnbase::filterEmptySpectra() # delete the empty spectra

fdt_sub <- fData(temp_dt) %>% as.data.table()
event_index <- fdt_sub[, which(
  (msLevel == 1) |
    (msLevel == 2 & precursorMZ <= 280)
)]

data_example <- temp_dt[event_index]
x <- fData(data_example)

data_example %>%
  MSnbase::writeMSData(., "./inst/extdata/test_data.mzml")

file <- "./inst/extdata/test_data.mzML"
data_example <- MSnbase::readMSData(file, mode = "onDisk")

# # to order the spectrum indexes from 1:length(spectra)
# fData(data_example)$spectrum <- seq(1, nrow(fData(data_example)), 1)  

usethis::use_data(data_example, overwrite = TRUE)
