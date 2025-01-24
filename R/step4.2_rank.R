# Rank of factorization method:
# At every iteration p_j (peak processing), we will search for the already processed peaks that co-elute (and not just overlap with rtmin or rtmax) of p_j.
# rank = nb of "processed co-eluting" peaks with p_j before p_j (unique peaks and not several peaks related to the same compound)
# we do not need to filter the mz (think why?)
# The key to a succession of this strategy is the decreasing order of the peak intensity at the apexes

#' find rank of a mixed data
#'
#' @inheritParams extract_ms_matrix.f
#' @param max_r max rank allowed
#' @param min_rt `numeric`
#' @param max_rt `numeric`
#'
#' @return rank of factorization (number of pure compounds)
find_rank <- function(ms1_peaks.df, peak.idx, min_rt, max_rt, max_r){
  if( peak.idx == 1 ){
    return( min(10, max_r) )
  } else {
    ms1_peaks_sub <- ms1_peaks.df[ 1:(peak.idx-1), ] # processed peaks
    
    # 1st method
    # peak_rtmin <- ms1_peaks.df[peak.idx, "rtmin"]
    # peak_rtmax <- ms1_peaks.df[peak.idx, "rtmax"]
    
    left_coelute <- which( ms1_peaks_sub[, "rtmax"] > min_rt & ms1_peaks_sub[, "rtmin"] < min_rt )
    right_coelute <- which( ms1_peaks_sub[, "rtmin"] < max_rt & ms1_peaks_sub[, "rtmax"] > max_rt )
    
    # # 2nd method
    # left_coelute <- which( ms1_peaks_sub[, "rtmax"] > window_rt_min & ms1_peaks_sub[, "rtmin"] < window_rt_min )
    # right_coelute <- which( ms1_peaks_sub[, "rtmin"] < window_rt_max & ms1_peaks_sub[, "rtmax"] > window_rt_max )
    
    # count the co-eluting peaks
    coelute <- unique(c(left_coelute, right_coelute))   
    coelute_peaks <- ms1_peaks_sub[coelute, "is_ion"]
    r <- length(unique(coelute_peaks))
    
    if(r == 0){
      return( min(10, max_r) )  # I should find a way
    } else {
      r <- r + 5
      r <- min(r, 10) 
    }
  }
  
  if( r > max_r ){
    r <- max_r  }
  
  return(r)
}


# check the peak from the raw data
is_true_peak <- function(chromatogram = chrom) {
  if (!is.data.frame(chromatogram)) {
    stop("empty data, no chromatogram.")  }
  
  apex_index <- which.max(chromatogram$intensity)
  if (apex_index == 1 || apex_index == nrow(chromatogram)) {
    return(FALSE)  }
  
  intensities <- chromatogram$intensity
  apex_intensity <- intensities[apex_index]
  
  left_valid <- all(intensities[(apex_index - 2):(apex_index - 1)] < apex_intensity) &&  # no need for the 1st condition !!
    (intensities[apex_index - 2] < intensities[apex_index - 1])
  right_valid <- all(intensities[(apex_index + 1):(apex_index + 2)] < apex_intensity) &&
    (intensities[apex_index + 1] > intensities[apex_index + 2])
  
  return(left_valid && right_valid)
}
