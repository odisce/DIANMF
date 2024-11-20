# Rank of factorization method:
# At every iteration p_j (peak processing), we will search for the already processed peaks that co-elute (and not just overlap with rtmin or rtmax) of p_j.
# rank = nb of "processed co-eluting" peaks with p_j before p_j (unique peaks and not several peaks related to the same compound)
# we do not need to filter the mz (think why?)
# The key to a succession of this strategy is the decreasing order of the peak intensity at the apexes

#' find rank of a mixed data
#'
#' @inheritParams extract_ms_matrix.f
#' @inheritParams check_ms1_ions 
#' @param max_r max rank allowed 
#'
#' @return rank of factorization (number of pure compounds)
find_rank <- function(ms1_peaks.df, peak.idx, rt_prec, rt_tol = 5, max_r){
  if( peak.idx == 1 ){
    r <- 3
  } else {
    ms1_peaks_sub <- ms1_peaks.df[ 1:(peak.idx-1), ] # processed peaks
    peak_rtmin <- ms1_peaks.df[peak.idx, "rtmin"]
    peak_rtmax <- ms1_peaks.df[peak.idx, "rtmax"]
    
    left_coelute <- which( ms1_peaks_sub[, "rtmax"] > peak_rtmin & ms1_peaks_sub[, "rtmin"] < peak_rtmin & (rt_prec-ms1_peaks_sub[, "rt"]) < rt_tol )
    right_coelute <- which( ms1_peaks_sub[, "rtmin"] < peak_rtmax & ms1_peaks_sub[, "rtmax"] > peak_rtmax & (ms1_peaks_sub[, "rt"] - rt_prec) < rt_tol )
    
    coelute <- unique(c(left_coelute, right_coelute))   
    x <- ms1_peaks_sub[coelute, ]
    coelute_peaks <- ms1_peaks_sub[coelute, "is_ion"]
    r <- length(unique(coelute_peaks))
    
    if(r == 0){
      return(3)  # this can be improved (I have something in mined)
    } else {
      r <- r + 3
    }
    if( r > 10 | r > max_r ){
      return(0)  }
  }
  
  return(r)
}
