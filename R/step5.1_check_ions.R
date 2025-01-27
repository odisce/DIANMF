#' Check which ions are also detected as MS1 peaks.
#'
#' @inheritParams extract_ms1_pure_spectrum
#' @param comp_ms1 index of the corresponding peak compound.
#' @inheritParams extract_ms_matrix.f
#' @param rt_prec `numeric(1)` apex retention time.
#' @param rt_tol_ions `numeric(1)` retention time tolerance, to check if an ion is also a detected peak from the peaks data frame.
#'
#' @return `numeric` MS1 peaks indexes.
check_ms1_ions <- function(W_ms1, comp_ms1, ms1_peaks.df, rt_prec, rt_tol_ions = 2){
  
  # ions_maybe_peaks <- which( W_ms1[,comp_ms1] >= rowSums(W_ms1) )  # 1st method
  # ions_maybe_peaks <- which( W_ms1[,comp_ms1] >= rowSums(W_ms1[, -comp_ms1]) ) # 2nd method, bad
  
  # ions_maybe_peaks <- apply(W_ms1, 1, function(row) {  # 3rd method: very bad
  #   row[comp_ms1] >= ( max(row[-comp_ms1]) + 5000 )
  # })
  # ions_maybe_peaks_indices <- which(ions_maybe_peaks)
  
  # ions_maybe_peaks <- as.numeric(names(ions_maybe_peaks))
  
  # ions_maybe_peaks <- as.numeric(rownames(W_ms1)) #4rt method (use the retention time info)
  
  # 5th method, use the frequency of each ion; i.e. how many times it is part from other components
  # if it is involved in many other sources 
  non_zero_counts <- rowSums(W_ms1[, -comp_ms1] != 0)
  ions_maybe_peaks <- as.numeric(rownames(W_ms1)[non_zero_counts < max( floor( (ncol(W_ms1) -1) /2 )-1, 2 ) ])

  ions_peaks <- c()
  for (i in 1:length(ions_maybe_peaks)) {
    idx <- which( ions_maybe_peaks[i] >= ms1_peaks.df[, 'mzmin'] & 
                    ions_maybe_peaks[i] <= ms1_peaks.df[, 'mzmax'] &
                    abs(rt_prec - ms1_peaks.df[, 'rt']) < rt_tol_ions  )
    
    if( length(unname(idx)) > 0 ){
      ions_peaks <- c(ions_peaks,  idx)
    }
  }
  
  return(ions_peaks)
}

# &
#   abs(rt_min - ms1_peaks.df[, 'rtmin']) < rt_tol_ions &
#   abs(rt_max - ms1_peaks.df[, 'rtmax']) < rt_tol_ions