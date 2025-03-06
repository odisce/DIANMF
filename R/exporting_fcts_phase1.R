# this script contains exporting functions 
# features_list <- readRDS("~/DIA_NMF_R_package/Results_ms1_then_ms2/features_list.rds")

#' Export MS1 pure spectra
#'
#' @param features.l `list` created by `DIANMF::DIANM.f`
#'
#' @return `Spectra` object of MS1 pure spectra.
#' @export
#' 
#' @import dplyr
#' @importFrom Spectra Spectra
#' @importFrom BiocGenerics Reduce
#' @importFrom S4Vectors DataFrame
exportMS1Spectra <- function(features.l){
  
  ms1Spectra <- lapply(1:length(features.l), function(f){
    ms1PureSpectra <- features.l[[f]]$ms1_pure_spectra
    ms1PureSpectra <- ms1PureSpectra[, c("rank", "mz", "value")]
    ms1PureSpectra <- ms1PureSpectra %>%
      group_by(rank) %>%
      group_split()
    return(ms1PureSpectra)
  })
  
  ms1Spectra <- BiocGenerics::Reduce(c, ms1Spectra)
  
  spd <- S4Vectors::DataFrame(
    msLevel = rep(1L, length(ms1Spectra))
  )
  spd$mz <- lapply(ms1Spectra, function(df) df$mz)
  spd$intensity <- lapply(ms1Spectra, function(df) df$value)
  sps <- Spectra::Spectra(spd) # the default backend MsBackendMemory 

  return(sps)
}
# ms1Spect <- exportMS1Spectra(features.l)
# 
# # export the spectra in a specific format .mgf or .msp
# library(MsBackendMgf)
# fl_mgf <- "output.mgf"
# ss_mgf <- export(ms1Spect, backend = MsBackendMgf(), file = fl_mgf )
# my_mgf_file <- MsBackendMgf::readMgf(fl_mgf)
# my_mgf_file$mz@listData
# my_mgf_file$intensity@listData
# 
# library(MsBackendMsp)  
# fl_msp <- "output.msp"  # the output file name
# ss_msp <- export(ms1Spect, MsBackendMsp(), file = fl_msp)
# readLines(fl_msp, n = 10)

# spectraVariables(sps)
# msLevel(sps)
# print(object.size(sps), units = "Mb")
# mz(sps[10])
# intensity(sps[10])
# ------------------------------------------------------------------------------


#' Export MS2 pure spectra
#'
#' @inheritParams exportMS1Spectra
#'
#' @return `Spectra` object of MS2 pure spectra.
#' @export
#'
#' @import dplyr
#' @importFrom Spectra Spectra
#' @importFrom BiocGenerics Reduce
#' @importFrom S4Vectors DataFrame
exportMS2Spectra <- function(features.l){
 
  ms2Spectra <- lapply(1:length(features.l), function(f){
    ms2PureSpectra <- features.l[[f]]$ms2_pure_spectra
    ms2PureSpectra <- ms2PureSpectra %>%
      group_by(rank) %>%
      group_split()
    return(ms2PureSpectra)
  })
  
  ms2Spectra <- BiocGenerics::Reduce(c, ms2Spectra)
  
  spd <- S4Vectors::DataFrame(
    msLevel = rep(2L, length(ms2Spectra))
  )
  spd$mz <- lapply(ms2Spectra, function(df) df$mz)
  spd$intensity <- lapply(ms2Spectra, function(df) df$value)
  sps <- Spectra::Spectra(spd) 
  
  return(sps)
}
# ms2Spect <- exportMS2Spectra(features.l)

# ------------------------------------------------------------------------------

#' Export MS2 pure spectra from different isolation windows
#'
#' @inheritParams exportMS1Spectra
#'
#' @return `list` of `Spectra`
#' @export
#'
#' @import dplyr
#' @importFrom Spectra Spectra
#' @importFrom S4Vectors DataFrame
exportMS2Spectra_isolW <- function(features.l){

  ms2Spectra <- lapply(1:length(features.l), function(f){
    ms2PureSpectra <- features.l[[f]]$ms2_pure_spectra
    ms2PureSpectra <- ms2PureSpectra %>% 
      dplyr::arrange(isolationWindowTargetMz)
    nb_isolation_windows <- length(unique(ms2PureSpectra$isolationWindowTargetMz))
    
    ms2_spd_isoW <- S4Vectors::DataFrame(
      msLevel = rep(2L, nb_isolation_windows),
      id = seq(1, nb_isolation_windows, by = 1)
    )
    ms2_spectrum_isoW <- ms2PureSpectra %>%
      group_by(isolationWindowTargetMz) %>%
      group_split()
    
    ms2_spd_isoW$mz <- lapply(ms2_spectrum_isoW, function(df) df$mz)
    ms2_spd_isoW$intensity <- lapply(ms2_spectrum_isoW, function(df) df$value)
    
    ms2_sp_isoW <- Spectra::Spectra(ms2_spd_isoW)
    
    return(ms2_sp_isoW)
  })
  
  return(ms2Spectra)
}
# ms2Spect <- exportMS2Spectra_isolW(features.l)
# mz(ms2Spect[[1]])
