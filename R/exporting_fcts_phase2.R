# this script contains functions that export for every feature its spectra

#' Export Features MS spectra 
#'
#' @param ms1features.df `data.frame` ms1 features obtained by xcms or with `DIANMF::extract_xcms_features`. 
#' @inheritParams exportMS1Spectra 
#'
#' @return `list` of `Spectra` objects for every feature.
#' @export
#'
#' @import dplyr
#' @importFrom Spectra Spectra
#' @importFrom S4Vectors DataFrame
exportMSFeatures <- function(ms1features.df, features.l){
  
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
