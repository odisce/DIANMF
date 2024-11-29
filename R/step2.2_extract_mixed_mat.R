#' Check if the vector has at least 4 consecutive non-zero values.
#'
#' @param row vector of `numeric` to evaluate.
#'
#' @return `Logical` `TRUE` if the vector contains at least 4  consecutive non-zero values, `FALSE` otherwise.
has_four_consecutive_non_zero <- function(row) {
  # Find consecutive non-zero values
  non_zero_streaks <- rle(row != 0)
  any(non_zero_streaks$lengths[non_zero_streaks$values] >= 4)
}
# mat <- matrix(c(
#   7, 1, 2, 3, 0,
#   2, 0, 1, 0, 0,
#   4, 3, 4, 7, 0,
#   1, 0, 2, 0, 3
# ), nrow = 4, byrow = TRUE)
# row_filter <- apply(mat, 1, has_four_consecutive_non_zero)
# filtered_mat <- mat[row_filter, , drop = FALSE]


#' Extract MS1 or MS2 mixed matrix of a peak.
#'
#' @param peak.idx `numeric(1)` MS1 peak index.
#' @param ms1_peaks.df MS1 peaks `data.frame`.  
#' @inheritParams detect_peaks_by_xcms
#' @param ppm.n `numeric(1)` defining the m/z tolerated deviation in parts per million (ppm).
#' @inheritParams extract_eics
#' @param iso_win_index `numeric(1)` SWATH isolation window index; just for MS2 level.
#'
#' @return Mixed `matrix` related to peak.idx.
#' 
#' @export
#' 
#' @importFrom data.table as.data.table 
#' @importFrom MSnbase filterMsLevel filterRt fData filterIsolationWindow spectra rtime
#' @import magrittr
#' @importFrom stats spectrum
extract_ms_matrix.f <- function(peak.idx, ms1_peaks.df, rawData.onDiskMSnExp, ppm.n = 7, rt_index = TRUE, mz_range = TRUE, iso_win_index = NULL){
  
  # raw_fdata <- fData(rawData.onDiskMSnExp) %>% as.data.table()
  rawData_ms1 <- MSnbase::filterMsLevel(rawData.onDiskMSnExp, 1L) %>%
    MSnbase::filterRt(., c(ms1_peaks.df[peak.idx, 'rtmin'], ms1_peaks.df[peak.idx, 'rtmax'])) %>%
    MSnbase::fData(.) %>%
    as.data.table(.)
  
  if (!is.null(iso_win_index)) {
    iso_win_mat <- isolationWindows.range(rawData.onDiskMSnExp) %>%
      as.data.table()
    iso_mass <- mean(unlist(iso_win_mat[iso_win_index,]))
    rawData_sub <- MSnbase::filterIsolationWindow(rawData.onDiskMSnExp, iso_mass) %>%
      fData() %>%
      as.data.table()
    ms1_L <- FALSE
    ## Select the closest MS2 scans from MS1 inside the peak
    ## Extract MS1 scans belonging to the peak
    scans_info <- rawData_ms1 %>% {
      .[msLevel == 1, {
        ms1_rt <- retentionTime
        rawData_sub[which.min(abs(retentionTime - ms1_rt)),]
      }, by = .(spectrum_ms1 = spectrum)]
    }
    
  } else {
    scans_info <- rawData_ms1
    ms1_L <- TRUE
  }
  
  # idx.ms <- which(raw_fdata$spectrum %in% scans_info[, spectrum]);
  idx.ms <- scans_info[, spectrum]

  mz <- ms1_peaks.df[peak.idx, 'mz'];  # mz of this peak
  peak_apex_rt <- ms1_peaks.df[peak.idx, 'rt']

  idx.apex.ms <- which.min(abs(rawData_ms1$retentionTime - peak_apex_rt))
  
  # get the ms scans related to the ms peak
  spec.exp_ms <- MSnbase::spectra(rawData.onDiskMSnExp[idx.ms])
  
  spec.exp_rt <- sapply(spec.exp_ms, rtime)
  ## get retention time of each spectra: 
  spec.exp_ms <- lapply(1:length(spec.exp_ms), function(i){
    res <- data.table(
      "mz" = mz(spec.exp_ms[[i]]),
      "intensity" = intensity(spec.exp_ms[[i]])
      # "rt" = rtime(spec.exp_ms1[[i]])
    )
    return(res)
  });
  
  if(!is.null(mz_range) && isTRUE(ms1_L)){ # the isolation window range of where this MS1 peak was fragmented
    info.swath <- isolationWindows.range(rawData.onDiskMSnExp)
    idx.swath <- which( info.swath[, 'lowerMz'] <= mz & info.swath[, 'upperMz'] >= mz )
    mz_range <- info.swath[idx.swath, ]
  };
  if (isFALSE(ms1_L) & !is.null(mz_range)) {
    warning("Selected MS2 data, the mz_range was set to NULL")
    mz_range <- NULL
  }
  
  ms_mixed_matrix <- extract_eics(spectra_list = spec.exp_ms, ppm.n = ppm.n,
                                   apex_index = idx.apex.ms,
                                   rt_index = spec.exp_rt, mz_range = mz_range);
  
  # filter the rows; Keep only the ones that have at least 4 consecutive non-zero values
  row_filter <- apply(ms_mixed_matrix, 1, has_four_consecutive_non_zero)
  ms_mixed_matrix <- ms_mixed_matrix[row_filter, , drop = FALSE]
  
  return(ms_mixed_matrix)
}


#' Extract MS2 matrices of a peak.
#'
#' @inheritParams extract_ms_matrix.f
#' @param ms1_pure_spectrum `data.frame` MS1 pure spectrum of the peak. 
#' @inheritParams detect_peaks_by_xcms
#' @param info.swath `data.frame` SWATH isolation windows.
#'
#' @return `list` of ms2 matrices, every matrix is from one isolation window.
#' 
#' @export
#' 
#' @importFrom data.table %between%
extract_ms2_matrices <- function(peak.idx, ms1_peaks.df, ppm.n, ms1_pure_spectrum, rawData.onDiskMSnExp, info.swath){
  
  min_mz <- min(ms1_pure_spectrum$mz_value);
  max_mz <- max(ms1_pure_spectrum$mz_value);
  idx.swath <- which(info.swath$lowerMz < max_mz | info.swath$upperMz > min_mz); 
  
  # so from every isolation window in idx.swath we will extract the ms2 data
  res_ms2 <- lapply(idx.swath, function(i){
    ms2_mat <- extract_ms_matrix.f(peak.idx = peak.idx, ms1_peaks.df = ms1_peaks.df, rawData.onDiskMSnExp = rawData.onDiskMSnExp,
                                   ppm.n = ppm.n, rt_index = TRUE, mz_range = NULL, iso_win_index = i);
    return(ms2_mat)
  });

  return(res_ms2)
}


#-------------------------------------------------------------------------------------------------------------------------------------
# old functions
# extract_ms1_matrix.f <- function(idx.pg, eics_peaks.mat, rawData.onDiskMSnExp, ppm = 7, rt_index = TRUE, mz_range = TRUE){
#   
#   scans_info <- MSnbase::fData(rawData.onDiskMSnExp);
#   ms1Scans_info <- scans_info[scans_info$msLevel == '1', ];
#   
#   idx.ms1 <- which(ms1Scans_info$retentionTime >= eics_peaks.mat[idx.pg, 'rtmin'] &
#                      ms1Scans_info$retentionTime <= eics_peaks.mat[idx.pg, 'rtmax'] );
#   idx.ms1 <- ms1Scans_info[idx.ms1, 'spectrum'];
#   
#   mz <- eics_peaks.mat[idx.pg, 'mz'];  # mz of this peak
#   
#   peak_apex_rt <- eics_peaks.mat[idx.pg, 'rt'];
#   ms1Scans_peak <- scans_info[ idx.ms1, ];
#   
#   idx.apex.ms1 <- which.min(abs(ms1Scans_peak$retentionTime - peak_apex_rt));
#   
#   # get the ms1 scans related to the ms1 peak
#   spec.exp_ms1 <- MSnbase::spectra(rawData.onDiskMSnExp[idx.ms1]);
#   
#   spec.exp_rt <- sapply(spec.exp_ms1, rtime)
#   ## get retention time of each spectra: 
#   spec.exp_ms1 <- lapply(1:length(spec.exp_ms1), function(i){
#     res <- data.table(
#       "mz" = mz(spec.exp_ms1[[i]]),
#       "intensity" = intensity(spec.exp_ms1[[i]])
#       # "rt" = rtime(spec.exp_ms1[[i]])
#     )
#     return(res)
#   });
#   
#   if(!is.null(mz_range)){ # the isolation window range of where this MS1 peak was fragmented
#     info.swath <- isolationWindows.range(rawData.onDiskMSnExp)
#     idx.swath <- which( info.swath[, 'lowerMz'] <= mz & info.swath[, 'upperMz'] >= mz )
#     mz_range <- info.swath[idx.swath, ]
#   };
#   
#   ms1_mixed_matrix <- extract_eics(spectra_list = spec.exp_ms1, ppm = ppm,
#                                    apex_index = idx.apex.ms1,
#                                    rt_index = spec.exp_rt, mz_range = mz_range);
#   
#   return(ms1_mixed_matrix)
# }
# spectra_list <- spec.exp_ms1
# ppm <- 7
# apex_index <- idx.apex.ms1
# mz_range <- isol_window_mz_range
# mz_range <- NULL
# rt_index <- spec.exp_rt


# extract_ms2_matrix.f <- function(idx.pg, eics_peaks.mat, rawData.onDiskMSnExp, ppm = 7, rt_index = TRUE, mz_range = TRUE){
# 
#   scans_info <- MSnbase::fData(rawData.onDiskMSnExp);
#   ms1Scans_info <- scans_info[scans_info$msLevel == '1', ];
# 
#   idx.ms1 <- which(ms1Scans_info$retentionTime >= eics_peaks.mat[idx.pg, 'rtmin'] &
#                      ms1Scans_info$retentionTime <= eics_peaks.mat[idx.pg, 'rtmax'] );
#   idx.ms1.ext <- c( idx.ms1[1]-1, idx.ms1, idx.ms1[length(idx.ms1)] +1 )
#   idx.ms1 <- ms1Scans_info[idx.ms1, 'spectrum'];
#   idx.ms1.ext <- ms1Scans_info[idx.ms1.ext, 'spectrum'];
#   
#   mz <- eics_peaks.mat[idx.pg, 'mz'];  # mz of this peak
#   
#   info.swath <- isolationWindows.range(rawData.onDiskMSnExp)
#   idx.swath <- which( info.swath[, 'lowerMz'] <= mz & info.swath[, 'upperMz'] >= mz )
#   mz_range <- info.swath[idx.swath, ]
#   
#   # extract the MS2 scans related to these ms1 scans and from the isolation window: idx.swath
#   idx.ms2 <- idx.ms1 + idx.swath
#   idx.ms2.ext <- idx.ms1.ext + idx.swath
#   
#   ms1Scans_peak <- scans_info[ idx.ms1, ];
#   peak_apex_rt <- eics_peaks.mat[idx.pg, 'rt'];
#   
#   if( idx.swath <=5 ){  # as MS1 data extracted, nothing should be changed
#     
#     idx.apex.ms1 <- which.min(abs(ms1Scans_peak$retentionTime - peak_apex_rt));  
#     
#     # extract the ms2 spectra
#     spec.exp_ms2 <- MSnbase::spectra(rawData.onDiskMSnExp[idx.ms2]);
#     spec.exp_rt <- sapply(spec.exp_ms2, rtime)
# 
#     spec.exp_ms2 <- lapply(1:length(spec.exp_ms2), function(i){
#       res <- data.table(
#         "mz" = mz(spec.exp_ms2[[i]]),
#         "intensity" = intensity(spec.exp_ms2[[i]])
#         # "rt" = rtime(spec.exp_ms2[[i]])
#       )
#       return(res)
#     });
#     
#     ms2_mixed_matrix <- extract_eics(spectra_list = spec.exp_ms2, ppm = ppm,
#                                          apex_index = idx.apex.ms1,
#                                          rt_index = spec.exp_rt, mz_range = mz_range);
#   } else {  # for every MS1 scan, its nearest MS2 scan is the MS2 scan from the previous MS1 scan  
#     
#     idx.apex.ms1 <- which( idx.ms1.ext == idx.ms1[which.min(abs(ms1Scans_peak$retentionTime - peak_apex_rt))]);
#     
#     spec.exp_ms2 <- MSnbase::spectra(rawData.onDiskMSnExp[idx.ms2.ext]);
#     spec.exp_rt <- sapply(spec.exp_ms2, rtime)
#     
#     spec.exp_ms2 <- lapply(1:length(spec.exp_ms2), function(i){
#       res <- data.table(
#         "mz" = mz(spec.exp_ms2[[i]]),
#         "intensity" = intensity(spec.exp_ms2[[i]])
#         # "rt" = rtime(spec.exp_ms2[[i]])
#       )
#       return(res)
#     });
#     
#     ms2_mixed_matrix <- extract_eics(spectra_list = spec.exp_ms2, ppm = ppm,
#                                          apex_index = idx.apex.ms1,
#                                          rt_index = spec.exp_rt, mz_range = mz_range);
#     
#     
#     ms2_mixed_matrix <- ms2_mixed_matrix[, 1:(ncol(ms2_mixed_matrix) - 2)];
#   }
#   
#   return(ms2_mixed_matrix)
# }

