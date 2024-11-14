# Rank of factorization methods
# At every iteration p_j (peak processing), we will search for the already processed peaks that overlap with rtmin or rtmax of p_j.
# rank = nb of "processed" peaks before p_j
# we do not need to filter the mz (think why?)
# The key to a succession of this strategy is the decreasing order of the peak intensity at the apexes

find_rank <- function(ms1_peaks.df, peak.idx, rt_tol = 3, max_r){
  if( peak.idx == 1 ){
    r <- 3
  } else {
    ms1_peaks_sub <- ms1_peaks.df[ 1:peak.idx, ] # processed peaks
    peak_rtmin <- ms1_peaks.df[peak.idx, "rtmin"]
    peak_rtmax <- ms1_peaks.df[peak.idx, "rtmax"]
    
    left_coelute <- which( ms1_peaks_sub[, "rtmax"] > peak_rtmin + rt_tol & ms1_peaks_sub[, "rtmin"] < peak_rtmin - rt_tol )
    right_coelute <- which( ms1_peaks_sub[, "rtmin"] < peak_rtmax - rt_tol & ms1_peaks_sub[, "rtmax"] > peak_rtmax + rt_tol )
    
    coelute <- unique(c(left_coelute, right_coelute))      
    coelute_peaks <- ms1_peaks_sub[coelute, "is_ion"]
    r <- length(unique(coelute_peaks))
    
    if( r > 10 | r > max_r ){
      r <- 0  }
    if(r == 0){
      # print("no co-eluting peaks")
      r <- 3  # this can be improved (peak7)
    } else {
      r <- r + 2
      print(paste(r, "co-eluting peaks"))
    }
  }
  
  return(r)
}