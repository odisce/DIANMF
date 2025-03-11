# the following function, extract for every identified pure source in all iterations one Spectra object;
# this object contains 1 MS1 pure spectra and the 10 MS2 pure spectra from different isolation windows.

#' Export MS pure spectra objects
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
exportMSSpectra <- function(features.l) {
  
  MSSpectra <- lapply(seq_along(features.l), function(iteration) {
    
    ms1PureSpectra <- features.l[[iteration]]$ms1_pure_spectra %>%
      arrange(rank) %>%
      group_by(rank) %>%
      group_split()
    
    ms2PureSpectra <- features.l[[iteration]]$ms2_pure_spectra %>%
      arrange(isolationWindowTargetMz) %>%
      group_by(rank) %>%
      group_split()
    
    info <- features.l[[iteration]]$xcms_assigned_sources %>%
      arrange(source) %>%
      group_by(source) %>%
      group_split()
    
    res <- lapply(seq_along(ms1PureSpectra), function(source_idx) {
      # MS1
      ms1_df <- ms1PureSpectra[[source_idx]]
      ms1_spd <- S4Vectors::DataFrame(
        msLevel = 1L,
        name = iteration,
        id = unique(ms1_df$rank),
        rtime = mean(info[[source_idx]]$rt),
        COMMENT = paste(info[[source_idx]]$xic_label, collapse = ";")
      )
      
      ms1_spd$mz <- list(ms1_df$mz)
      ms1_spd$intensity <- list(ms1_df$value)
      ms1_spectra <- Spectra::Spectra(ms1_spd)
      
      # MS2
      ms2_df <- ms2PureSpectra[[source_idx]]
      isolation_spectra <- lapply(unique(ms2_df$isolationWindowTargetMz), function(targetMz) {
        subset_df <- ms2_df[ms2_df$isolationWindowTargetMz == targetMz, ]  
        
        ms2_spd <- S4Vectors::DataFrame(
          msLevel = 2L,
          name = iteration,
          id = unique(subset_df$rank), 
          isolationWindowTargetMz = as.numeric(targetMz),  
          isolationWindowLowerMz = as.numeric(unique(subset_df$isolationWindowLowerMz)),
          isolationWindowUpperMz = as.numeric(unique(subset_df$isolationWindowUpperMz))  )
        
        ms2_spd$mz <- list(subset_df$mz)
        ms2_spd$intensity <- list(subset_df$value)
        return(Spectra::Spectra(ms2_spd))
      })
      
      isolation_spectra <- do.call(c, isolation_spectra)
      msSpectra <- c(ms1_spectra, isolation_spectra)
      return(msSpectra) 
    })
    
    return(res)
  })
  
  MSSpectra <- unlist(MSSpectra, recursive = FALSE)
  return(MSSpectra)
}
# pure_spectra <- exportMSSpectra(features.l) 
# 
# # msSpectra <- pure_spectra[[1]]
# library(MsBackendMgf )
# export(msSpectra, backend = MsBackendMgf(), file = "./output_spectra.mgf")
# sps <- msSpectra 
# spectraVariables(sps) # all 
# msLevel(sps)
# rtime(sps)
# mz(sps)
# intensity(sps)
# 
# ms1_Spectra <- filterMsLevel(msSpectra, 1L)
# ms2_Spectra <- filterMsLevel(msSpectra, 2L)
# 
# ms2_spec <- filterIsolationWindow(ms2_Spectra, mz = 115)
# mz(ms2_spec)
# intensity(ms2_spec)



# for every feature (or peak), we want one Spectra object, alot of redundancy here
# exportMSFeaturesSpectra <- function(ms1_features.df, features.l){
#   
#   res <- lapply(ms1_features.df$featureid, function(f_id){
#     # f_id <- "FT00001"
#     print(f_id)
#     iter <- ms1_features.df[featureid == f_id, ]$iteration %>%
#       strsplit(",") %>%               
#       unlist() %>%                     
#       as.numeric() %>%                
#       .[1]
#     
#     if( iter != 0 ){
#       sub_info <- features.l[[iter]]
#       
#       feature_msSpectra <- exportMSSpectra(features.l = list(sub_info))
#     }
#   })
#   
#   return(res)
# }