# functions to extract MS1 and MS2 mixed matrix/data

#' Extract MS1 data
#'
#' @param idx.pg 
#' @param eics_peaks.mat 
#' @param rawData.onDiskMSnExp 
#' @param ppm 
#' @param rt_index 
#' @param mz_range 
#'
#' @return mixed MS1 matrix
#' @export
#' @import data.table
#' @import magrittr
#' @import MSnbase
#'
extract_ms1_matrix.f <- function(idx.pg, eics_peaks.mat, rawData.onDiskMSnExp, ppm = 7, rt_index = TRUE, mz_range = TRUE){

  scans_info <- MSnbase::fData(rawData.onDiskMSnExp);
  ms1Scans_info <- scans_info[scans_info$msLevel == '1', ];
  
  idx.ms1 <- which(ms1Scans_info$retentionTime >= eics_peaks.mat[idx.pg, 'rtmin'] &
                     ms1Scans_info$retentionTime <= eics_peaks.mat[idx.pg, 'rtmax'] );
  idx.ms1 <- ms1Scans_info[idx.ms1, 'spectrum'];
  
  mz <- eics_peaks.mat[idx.pg, 'mz'];  # mz of this peak

  peak_apex_rt <- eics_peaks.mat[idx.pg, 'rt'];
  ms1Scans_peak <- scans_info[ idx.ms1, ];
  
  idx.apex.ms1 <- which.min(abs(ms1Scans_peak$retentionTime - peak_apex_rt));
  
  # get the ms1 scans related to the ms1 peak
  spec.exp_ms1 <- MSnbase::spectra(rawData.onDiskMSnExp[idx.ms1]);
  
  spec.exp_rt <- sapply(spec.exp_ms1, rtime)
  ## get retention time of each spectra: 
  spec.exp_ms1 <- lapply(1:length(spec.exp_ms1), function(i){
    res <- data.table(
      "mz" = mz(spec.exp_ms1[[i]]),
      "intensity" = intensity(spec.exp_ms1[[i]])
      # "rt" = rtime(spec.exp_ms1[[i]])
    )
    return(res)
  });
  
  if(!is.null(mz_range)){ # the isolation window range of where this MS1 peak was fragmented
    info.swath <- isolationWindows.range(rawData.onDiskMSnExp)
    idx.swath <- which( info.swath[, 'lowerMz'] <= mz & info.swath[, 'upperMz'] >= mz )
    mz_range <- info.swath[idx.swath, ]
  };
  
  ms1_mixed_matrix <- extract_eics(spectra_list = spec.exp_ms1, ppm = ppm,
                                       apex_index = idx.apex.ms1,
                                       rt_index = spec.exp_rt, mz_range = mz_range);
  
  return(ms1_mixed_matrix)
}
# spectra_list <- spec.exp_ms1
# ppm <- 7
# apex_index <- idx.apex.ms1
# mz_range <- isol_window_mz_range
# mz_range <- NULL
# rt_index <- spec.exp_rt

#' Extract MS2 data
#'
#' @param idx.pg 
#' @param eics_peaks.mat 
#' @param rawData.onDiskMSnExp 
#' @param ppm 
#' @param rt_index 
#' @param mz_range 
#'
#' @return mixed MS2 matrix
#' @export
#' @import data.table
#' @import magrittr
#' @import MSnbase
#'
extract_ms2_matrix.f <- function(idx.pg, eics_peaks.mat, rawData.onDiskMSnExp, ppm = 7, rt_index = TRUE, mz_range = TRUE){

  scans_info <- MSnbase::fData(rawData.onDiskMSnExp);
  ms1Scans_info <- scans_info[scans_info$msLevel == '1', ];

  idx.ms1 <- which(ms1Scans_info$retentionTime >= eics_peaks.mat[idx.pg, 'rtmin'] &
                     ms1Scans_info$retentionTime <= eics_peaks.mat[idx.pg, 'rtmax'] );
  idx.ms1.ext <- c( idx.ms1[1]-1, idx.ms1, idx.ms1[length(idx.ms1)] +1 )
  idx.ms1 <- ms1Scans_info[idx.ms1, 'spectrum'];
  idx.ms1.ext <- ms1Scans_info[idx.ms1.ext, 'spectrum'];
  
  mz <- eics_peaks.mat[idx.pg, 'mz'];  # mz of this peak
  
  info.swath <- isolationWindows.range(rawData.onDiskMSnExp)
  idx.swath <- which( info.swath[, 'lowerMz'] <= mz & info.swath[, 'upperMz'] >= mz )
  mz_range <- info.swath[idx.swath, ]
  
  # extract the MS2 scans related to these ms1 scans and from the isolation window: idx.swath
  idx.ms2 <- idx.ms1 + idx.swath
  idx.ms2.ext <- idx.ms1.ext + idx.swath
  
  ms1Scans_peak <- scans_info[ idx.ms1, ];
  peak_apex_rt <- eics_peaks.mat[idx.pg, 'rt'];
  
  if( idx.swath <=5 ){  # as MS1 data extracted, nothing should be changed
    
    idx.apex.ms1 <- which.min(abs(ms1Scans_peak$retentionTime - peak_apex_rt));  
    
    # extract the ms2 spectra
    spec.exp_ms2 <- MSnbase::spectra(rawData.onDiskMSnExp[idx.ms2]);
    spec.exp_rt <- sapply(spec.exp_ms2, rtime)

    spec.exp_ms2 <- lapply(1:length(spec.exp_ms2), function(i){
      res <- data.table(
        "mz" = mz(spec.exp_ms2[[i]]),
        "intensity" = intensity(spec.exp_ms2[[i]])
        # "rt" = rtime(spec.exp_ms2[[i]])
      )
      return(res)
    });
    
    ms2_mixed_matrix <- extract_eics(spectra_list = spec.exp_ms2, ppm = ppm,
                                         apex_index = idx.apex.ms1,
                                         rt_index = spec.exp_rt, mz_range = mz_range);
  } else {  # for every MS1 scan, its nearest MS2 scan is the MS2 scan from the previous MS1 scan  
    
    idx.apex.ms1 <- which( idx.ms1.ext == idx.ms1[which.min(abs(ms1Scans_peak$retentionTime - peak_apex_rt))]);
    
    spec.exp_ms2 <- MSnbase::spectra(rawData.onDiskMSnExp[idx.ms2.ext]);
    spec.exp_rt <- sapply(spec.exp_ms2, rtime)
    
    spec.exp_ms2 <- lapply(1:length(spec.exp_ms2), function(i){
      res <- data.table(
        "mz" = mz(spec.exp_ms2[[i]]),
        "intensity" = intensity(spec.exp_ms2[[i]])
        # "rt" = rtime(spec.exp_ms2[[i]])
      )
      return(res)
    });
    
    ms2_mixed_matrix <- extract_eics(spectra_list = spec.exp_ms2, ppm = ppm,
                                         apex_index = idx.apex.ms1,
                                         rt_index = spec.exp_rt, mz_range = mz_range);
    
    
    ms2_mixed_matrix <- ms2_mixed_matrix[, 1:(ncol(ms2_mixed_matrix) - 2)];
  }
  
  return(ms2_mixed_matrix)
}



