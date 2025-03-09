# this script contains exporting functions 

#' Export MS pure spectra
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
exportMSSpectra <- function(features.l){
  
  # msSpectra <- lapply(1:length(features.l), function(f){
  #   f <- 1
  #   
  #   ms1PureSpectra <- features.l[[f]]$ms1_pure_spectra
  #   ms1PureSpectra <- ms1PureSpectra[, c("rank", "mz", "value")]
  #   ms1PureSpectra <- ms1PureSpectra %>%
  #     arrange(rank) %>%
  #     group_by(rank) %>%
  #     group_split()
  #   
  #   info <- features.l[[f]]$xcms_assigned_sources
  #   info$source <- as.numeric(info$source)
  #   info <- info %>%
  #     arrange(source) %>%
  #     group_by(source) %>%
  #     group_split()
  #   
  #   ms2PureSpectra <- features.l[[f]]$ms2_pure_spectra
  #   ms2PureSpectra <- ms2PureSpectra %>% 
  #     dplyr::arrange(isolationWindowTargetMz) %>%
  #     group_by(rank) %>%
  #     group_split()
  #   
  #   res <- lapply(1:length(ms1PureSpectra), function(source_idx){
  #     source_idx <- 1
  #     
  #     spd <- S4Vectors::DataFrame(
  #       msLevel = c(1L, rep( 2L, length(ms2PureSpectra))),
  #       name = rep(f, length(ms2PureSpectra) +1),  # iteration index (first)
  #       id =  unique(ms2PureSpectra[[source_idx]]$rank)  # rank (source index)
  #     )
  #     
  #     spd$COMMENT <-  c(info[[source_idx]]$xic_label) # xcmspeaks
  #     spd_ms1$mz <- lapply(ms1PureSpectra, function(df) df$mz)
  #     spd_ms1$intensity <- lapply(ms1PureSpectra, function(df) df$value)
  #     sps_ms1 <- Spectra::Spectra(spd_ms1) # the default backend MsBackendMemory 
  #     
  #     
  #     
  #   })
  #   
  #   unlist(lapply(ms1PureSpectra, function(df) unique(df$rank)))
  #   
  #   
  # 
  #   
  #   # ms2 pure spectra from different isolation windows  
  #  
  #   
  #   nb_isolation_windows <- length(unique(ms2PureSpectra$isolationWindowTargetMz))
  #   
  #   ms2_spectrum_isoW <- ms2PureSpectra %>%
  #     group_by(rank) %>%
  #     group_split()
  #   
  #   ms2_spd_isoW <- S4Vectors::DataFrame(
  #     msLevel = rep(2L, length(ms2_spectrum_isoW)),
  #     id = unlist(lapply(ms2_spectrum_isoW, function(df) unique(df$rank)))  # rank
  #   )
  #  
  #   ms2_spd_isoW$mz <- lapply(ms2_spectrum_isoW, function(df) df$mz)
  #   ms2_spd_isoW$intesity <- lapply(ms2_spectrum_isoW, function(df) df$value)
  #   ms2_spd_isoW$isolationWindowTargetMz <- lapply(ms2_spectrum_isoW, function(df) df$isolationWindowTargetMz)
  #   ms2_spd_isoW$isolationWindowLowerMz <- lapply(ms2_spectrum_isoW, function(df) df$isolationWindowLowerMz)
  #   ms2_spd_isoW$isolationWindowUpperMz <- lapply(ms2_spectrum_isoW, function(df) df$isolationWindowUpperMz)
  #   isolationWindowUpperMz 
  #   ms2_sp_isoW <- Spectra::Spectra(ms2_spd_isoW)
  #   
  #   return(list(
  #     'sps_ms1' = sps_ms1,
  #     'sps_ms2' = ms2_spd_isoW
  #   ))
  # })
  # 
  # return(msSpectra)
}

# features.l <- readRDS("~/DIA_NMF_R_package/Results_ms1_then_ms2/features_list.rds")


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

# exportMS2Spectra <- function(features.l){
#  
#   ms2Spectra <- lapply(1:length(features.l), function(f){
#     ms2PureSpectra <- features.l[[f]]$ms2_pure_spectra
#     ms2PureSpectra <- ms2PureSpectra %>%
#       group_by(rank) %>%
#       group_split()
#     return(ms2PureSpectra)
#   })
#   
#   ms2Spectra <- BiocGenerics::Reduce(c, ms2Spectra)
#   
#   spd <- S4Vectors::DataFrame(
#     msLevel = rep(2L, length(ms2Spectra))
#   )
#   spd$mz <- lapply(ms2Spectra, function(df) df$mz)
#   spd$intensity <- lapply(ms2Spectra, function(df) df$value)
#   sps <- Spectra::Spectra(spd) 
#   
#   return(sps)
# }
# # ms2Spect <- exportMS2Spectra(features.l)

# ------------------------------------------------------------------------------

# exportMS2Spectra_isolW <- function(features.l){
# 
#   ms2Spectra <- lapply(1:length(features.l), function(f){
#     ms2PureSpectra <- features.l[[f]]$ms2_pure_spectra
#     ms2PureSpectra <- ms2PureSpectra %>% 
#       dplyr::arrange(isolationWindowTargetMz)
#     nb_isolation_windows <- length(unique(ms2PureSpectra$isolationWindowTargetMz))
#     
#     ms2_spd_isoW <- S4Vectors::DataFrame(
#       msLevel = rep(2L, nb_isolation_windows),
#       id = seq(1, nb_isolation_windows, by = 1)
#     )
#     ms2_spectrum_isoW <- ms2PureSpectra %>%
#       group_by(isolationWindowTargetMz) %>%
#       group_split()
#     
#     ms2_spd_isoW$mz <- lapply(ms2_spectrum_isoW, function(df) df$mz)
#     ms2_spd_isoW$intensity <- lapply(ms2_spectrum_isoW, function(df) df$value)
#     
#     ms2_sp_isoW <- Spectra::Spectra(ms2_spd_isoW)
#     
#     return(ms2_sp_isoW)
#   })
#   
#   return(ms2Spectra)
# }
# # ms2Spect <- exportMS2Spectra_isolW(features.l)
# # mz(ms2Spect[[1]])

# ------------------------------------------------------------------------------


# this script contains functions that export for every feature its spectra

#' Export Features MS spectra 
#'
#' @param ms1features.df `data.frame` ms1 features obtained by xcms or with `DIANMF::extract_xcms_features`. 
#' @inheritParams exportMSSpectra 
#'
#' @return `list` of `Spectra` objects for every feature.
#' @export
#'
#' @import dplyr
#' @importFrom Spectra Spectra
#' @importFrom S4Vectors DataFrame
exportMSFeaturesSpectra <- function(ms1features.df, features.l){
  
  res <- lapply(ms1features.df$featureid, function(f_id){
    print(f_id)
    # f_id <- "FT08993"
    iter <- ms1features.df[featureid == f_id, ]$iteration %>%
      strsplit(",") %>%               
      unlist() %>%                     
      as.numeric() %>%                
      .[1]
    
    if( iter != 0 ){
      sub_info <- features.l[[iter]]
      
      # ms1 spectrum -----------------------------
      ms1_spectrum <- sub_info$ms1_pure_spectra
      ms1_spd <- S4Vectors::DataFrame(
        msLevel = rep(1L, nrow(ms1_spectrum))
      )
      ms1_spd$mz <- list(ms1_spectrum$mz)
      ms1_spd$intensity <- list(ms1_spectrum$value)
      ms1_sp <- Spectra::Spectra(ms1_spd) 
      
      # ms2 spectrum -----------------------------
      ms2_spectrum <- sub_info$ms2_pure_spectra
      ms2_spd <- S4Vectors::DataFrame(
        msLevel = rep(2L, nrow(ms2_spectrum))
      )
      ms2_spd$mz <- list(ms2_spectrum$mz)
      ms2_spd$intensity <- list(ms2_spectrum$value)
      ms2_sp <- Spectra::Spectra(ms2_spd)
      
      # ms2 spectra (from different isolation window) ------------------
      ms2_spectrum <- ms2_spectrum %>%
        dplyr::arrange(isolationWindowTargetMz)
      
      nb_isolation_windows <- length(unique(ms2_spectrum$isolationWindowTargetMz))
      ms2_spd_isoW <- S4Vectors::DataFrame(
        msLevel = rep(2L, nb_isolation_windows),
        id = seq(1, nb_isolation_windows, by = 1)
      )
      
      ms2_spectrum_isoW <- ms2_spectrum %>%
        group_by(isolationWindowTargetMz) %>%
        group_split()
      
      ms2_spd_isoW$mz <- lapply(ms2_spectrum_isoW, function(df) df$mz)
      ms2_spd_isoW$intensity <- lapply(ms2_spectrum_isoW, function(df) df$value)
      
      ms2_sp_isoW <- Spectra::Spectra(ms2_spd_isoW)
      
      return(list(
        'ms1_spectrum' = ms1_sp,
        'ms2_spectrum' = ms2_sp,
        'ms2_spectra' = ms2_sp_isoW ))
    }
  })
  
  return(res)
}

# features.l <- features_list <- readRDS("~/DIA_NMF_R_package/Results_ms1_then_ms2/features_list.rds")
# ms1features.df <- ms1_features <- readRDS("~/DIA_NMF_R_package/Results3/ms1_features.rds")
# all_spct <- exportMSFeatures(ms1features.df[1:10, ], features.l)
