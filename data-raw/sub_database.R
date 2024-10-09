## code to prepare `sub_database` dataset goes here

# subset of the reference data base NewDB
NewDB <- readRDS("~/1workflow/DIA_NMF_workflow/Data/NewDB.rds")
class(NewDB)

sIndex.mz <- c(unname(unlist(NewDB$index[7975, 'mz'])), unname(unlist(NewDB$index[3924, 'mz'])),
               unname(unlist(NewDB$index[7678, 'mz'])), unname(unlist(NewDB$index[7963, 'mz'])))
sIndex.nce <- c(unname(unlist(NewDB$index[7975, 'NCE'])), unname(unlist(NewDB$index[3924, 'NCE'])),
                unname(unlist(NewDB$index[7678, 'NCE'])), unname(unlist(NewDB$index[7963, 'NCE'])))
sIndex.name <- c(unname(unlist(NewDB$index[7975, 'name'])), unname(unlist(NewDB$index[3924, 'name'])),
                 unname(unlist(NewDB$index[7678, 'name'])), unname(unlist(NewDB$index[7963, 'name'])))
sIndex.composition <- c(unname(unlist(NewDB$index[7975, 'composition'])), unname(unlist(NewDB$index[3924, 'composition'])),
                        unname(unlist(NewDB$index[7678, 'composition'])), unname(unlist(NewDB$index[7963, 'composition'])))
sIndex.id <- c(1,2,3,4)
sIndex.pol <- c("POS", "POS", "POS", "POS")

sIndex <- data.table(
  'id' = sIndex.id,
  'Polarity' = sIndex.pol,
  'mz' = sIndex.mz,
  'NCE' = sIndex.nce,
  'composition' = sIndex.composition,
  'name' = sIndex.name)

sSpectra <- list(NewDB$spectra[7975]$EXTRACT_MSMS1_d0ca6f846, NewDB$spectra[3924]$EXTRACT_MSMS1_9814bf7a1,
                 NewDB$spectra[7678]$EXTRACT_MSMS1_498f208d1, NewDB$spectra[7963]$EXTRACT_MSMS1_7017005c6)

sub_database <- list(
  'index' = sIndex,
  'spectra' = sSpectra )

usethis::use_data(sub_database, overwrite = TRUE)
